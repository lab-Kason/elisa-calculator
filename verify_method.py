"""Independent cross-check of calc.py, using a second implementation.

Why this file exists
--------------------
`test_calc.py` checks calc.py's behaviour (does it raise, fall back, flag
ranges).  This file instead re-derives every number from scratch with the
closed-form least-squares formulas and compares against calc.py, so a bug in
calc.py cannot hide behind the same code path.

It also proves the two sample datasets in `sample_data/` are self-consistent:
* A is built on the exact line corrected OD = 0.0008 * conc + 0.05, so the
  expected answer is known analytically (wastewater = 600 pg/mL at both
  dilutions).
* B is a realistic, noisy plate used to exercise out-of-range flags.

Run: python verify_method.py
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from calc import add_corrected, fit_linear, match_columns, quantify_samples, to_numeric

STD_COLS = ["Standard", "Concentration", "OD450", "OD600"]
SMP_COLS = ["Sample", "Dilution factor", "OD450", "OD600"]

failures: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    print(("PASS " if condition else "FAIL ") + name + (f"  ({detail})" if detail else ""))
    if not condition:
        failures.append(name)


def load(path: str, columns: list[str]) -> pd.DataFrame:
    return match_columns(pd.read_csv(path), columns)


# ---------------------------------------------------------------------------
# 1. independent implementation: ordinary least squares, closed form
#    slope = (n*Sxy - Sx*Sy) / (n*Sxx - Sx^2)
#    intercept = (Sy - slope*Sx) / n
# ---------------------------------------------------------------------------
def ols_closed_form(x: list[float], y: list[float]) -> tuple[float, float, float]:
    n = len(x)
    sx, sy = sum(x), sum(y)
    sxx = sum(v * v for v in x)
    sxy = sum(a * b for a, b in zip(x, y))
    denom = n * sxx - sx * sx
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    resid = [yi - (slope * xi + intercept) for xi, yi in zip(x, y)]
    ss_res = sum(r * r for r in resid)
    my = sy / n
    ss_tot = sum((yi - my) ** 2 for yi in y)
    r2 = 1.0 - ss_res / ss_tot if ss_tot else float("nan")
    return slope, intercept, r2


def expected_samples(smp: pd.DataFrame, slope: float, intercept: float) -> list[dict]:
    """Second implementation of the app's sample calculation."""
    out = []
    for _, row in smp.iterrows():
        y = row["corrected"]
        dilution = row["dilution"]
        if pd.isna(dilution) or float(dilution) <= 0:
            dilution = 1.0
        measured = (float(y) - intercept) / slope
        out.append({"name": row["name"], "dilution": float(dilution),
                    "measured": measured, "reported": measured * float(dilution)})
    return out


def prepare(standards_csv: str, samples_csv: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    std = load(standards_csv, STD_COLS).rename(
        columns={"Standard": "name", "Concentration": "conc", "OD450": "od450", "OD600": "od600"}
    )
    std = add_corrected(std)
    smp = load(samples_csv, SMP_COLS).rename(
        columns={"Sample": "name", "Dilution factor": "dilution", "OD450": "od450", "OD600": "od600"}
    )
    smp["dilution"] = to_numeric(smp["dilution"])
    smp = add_corrected(smp)
    return std, smp


def verify_dataset(paths: tuple[str, str], label: str) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    std, smp = prepare(*paths)
    std = std.loc[to_numeric(std["conc"]).notna() & std["corrected"].notna()].copy()
    std["conc"] = to_numeric(std["conc"])

    # calc.py's answer ...
    fit = fit_linear(std["conc"], std["corrected"])
    results = quantify_samples(smp, fit)

    # ... vs the independent closed-form answer
    slope_i, intercept_i, r2_i = ols_closed_form(
        std["conc"].tolist(), std["corrected"].tolist()
    )
    check(f"{label}: slope matches independent OLS",
          math.isclose(fit["slope"], slope_i, rel_tol=1e-12),
          f"calc={fit['slope']:.10g} independent={slope_i:.10g}")
    check(f"{label}: intercept matches independent OLS",
          math.isclose(fit["intercept"], intercept_i, rel_tol=1e-12),
          f"calc={fit['intercept']:.10g} independent={intercept_i:.10g}")
    check(f"{label}: R2 matches independent OLS",
          math.isclose(fit["r2"], r2_i, rel_tol=1e-10),
          f"calc={fit['r2']:.10g} independent={r2_i:.10g}")

    second = expected_samples(smp, slope_i, intercept_i)
    for got, want in zip(results.to_dict("records"), second):
        check(f"{label}: {want['name']} measured conc matches 2nd implementation",
              math.isclose(got["measured_conc"], want["measured"], rel_tol=1e-12, abs_tol=1e-9),
              f"calc={got['measured_conc']:.6f} independent={want['measured']:.6f}")
        check(f"{label}: {want['name']} reported conc matches 2nd implementation",
              math.isclose(got["reported_conc"], want["reported"], rel_tol=1e-12, abs_tol=1e-9),
              f"calc={got['reported_conc']:.6f} independent={want['reported']:.6f}")

    # round trip: reported conc -> absorbances -> should reproduce the raw OD
    for row in results.to_dict("records"):
        back = intercept_i + slope_i * (row["reported_conc"] / row["dilution"])
        check(f"{label}: round trip {row['name']}",
              math.isclose(back, row["corrected"], rel_tol=1e-9, abs_tol=1e-12),
              f"reconstructed corrected OD={back:.6f} raw={row['corrected']:.6f}")
    return fit, std, results


# ---------------------------------------------------------------------------
# 2. dataset A: analytic ground truth
# ---------------------------------------------------------------------------
print("=" * 72)
print("Dataset A — clean plate (exact line: corrected OD = 0.0008*conc + 0.05)")
print("=" * 72)
fit_a, std_a, res_a = verify_dataset(
    ("sample_data/example_A_clean.csv", "sample_data/example_A_samples_clean.csv"), "A"
)
check("A: slope is exactly 0.0008", math.isclose(fit_a["slope"], 0.0008, rel_tol=1e-9),
      f"{fit_a['slope']:.12g}")
check("A: intercept is exactly 0.05", math.isclose(fit_a["intercept"], 0.05, abs_tol=1e-12),
      f"{fit_a['intercept']:.12g}")
check("A: R2 is 1", math.isclose(fit_a["r2"], 1.0, abs_tol=1e-12), f"{fit_a['r2']:.12f}")
check("A: original sample = 600 pg/mL (analytic)",
      math.isclose(res_a.loc[0, "reported_conc"], 600.0, rel_tol=1e-6),
      f"{res_a.loc[0, 'reported_conc']:.6f}")
check("A: 1:2 diluted sample still reports 600 pg/mL after x2 dilution",
      math.isclose(res_a.loc[1, "reported_conc"], 600.0, rel_tol=1e-6),
      f"{res_a.loc[1, 'reported_conc']:.6f}")
print("\nTable A (computed by calc.py):")
print(res_a[["name", "dilution", "corrected", "measured_conc", "reported_conc", "status"]]
      .to_string(index=False, float_format=lambda v: f"{v:.4f}"))

# ---------------------------------------------------------------------------
# 3. dataset B: noisy plate + out-of-range handling
# ---------------------------------------------------------------------------
print()
print("=" * 72)
print("Dataset B — noisy plate (out-of-range / missing OD600 / replicates)")
print("=" * 72)
fit_b, std_b, res_b = verify_dataset(
    ("sample_data/example_B_messy.csv", "sample_data/example_B_samples_messy.csv"), "B"
)
check("B: R2 is realistic (0.99-1.0)", 0.99 < fit_b["r2"] <= 1.0, f"{fit_b['r2']:.6f}")
statuses = res_b["status"].tolist()
check("B: one sample above the curve", any("高于" in s for s in statuses), str(statuses))
check("B: one sample below the curve", any("低于" in s for s in statuses), str(statuses))
missing_od600 = res_b.loc[res_b["name"] == "OD600 missing demo", "corrected"].iloc[0]
check("B: blank OD600 falls back to OD450 (0.501)", math.isclose(missing_od600, 0.501, abs_tol=1e-12),
      f"{missing_od600:.6f}")

# the same wastewater measured at 1x / 1:5 / 1:50 must report the same
# concentration (true value 2200 pg/mL); the 1:50 well is out of range, so it
# is allowed to drift a little
reported_b = res_b.set_index("name")["reported_conc"]
check("B: original dilution reports ~2200 pg/mL",
      math.isclose(reported_b["Wastewater original"], 2200, rel_tol=0.01), f"{reported_b['Wastewater original']:.1f}")
check("B: 1:5 dilution reports the same ~2200 pg/mL",
      math.isclose(reported_b["Wastewater 1:5"], 2200, rel_tol=0.01), f"{reported_b['Wastewater 1:5']:.1f}")
check("B: 1:50 out-of-range dilution still ~2200 pg/mL (within 2%)",
      math.isclose(reported_b["Wastewater 1:50"], 2200, rel_tol=0.02), f"{reported_b['Wastewater 1:50']:.1f}")
check("B: replicate A ~730 pg/mL",
      math.isclose(reported_b["Replicate A"], 730, rel_tol=0.02), f"{reported_b['Replicate A']:.1f}")
check("B: replicate B ~745 pg/mL",
      math.isclose(reported_b["Replicate B"], 745, rel_tol=0.02), f"{reported_b['Replicate B']:.1f}")
check("B: missing-OD600 sample ~600 pg/mL",
      math.isclose(reported_b["OD600 missing demo"], 600, rel_tol=0.01), f"{reported_b['OD600 missing demo']:.1f}")
print("\nTable B (computed by calc.py):")
print(res_b[["name", "dilution", "corrected", "measured_conc", "reported_conc", "status"]]
      .to_string(index=False, float_format=lambda v: f"{v:.4f}"))

# ---------------------------------------------------------------------------
print()
print("=" * 72)
if failures:
    print(f"FAILED: {len(failures)} check(s): {failures}")
    raise SystemExit(1)
print("METHOD VERIFIED — every calc.py number matches an independent implementation.")
print("=" * 72)
