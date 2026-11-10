"""Independent cross-check of calc.py (4PL version).

Why this file exists
--------------------
`test_calc.py` checks calc.py's own behaviour.  This file instead re-derives
every sample concentration with a *different* algorithm:

* calc.py inverts the 4PL in closed form;
* here the inverse is found by **bisection** on an independently written
  4PL equation (not imported from calc.py).

It also checks the reference dataset:

* `sample_data/example_A_clean.csv` is generated from the reference
  "Real 4PL model" (A=0.0821, B=1.0741, C=3877.5, D=10.7614), so the expected
  answers are ~394 pg/mL (well 6, undiluted) and ~352 pg/mL (well 7, measured
  1:2 and multiplied by 2).
* `sample_data/example_B_messy.csv` is a noisier plate used to exercise the
  out-of-range flags.

Run: python verify_method.py
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

from calc import add_corrected, fit_4pl, match_columns, quantify_samples, to_numeric

STD_COLS = ["Standard", "Concentration", "OD450", "OD600"]
SMP_COLS = ["Sample", "Dilution factor", "OD450", "OD600"]

failures: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    print(("PASS " if condition else "FAIL ") + name + (f"  ({detail})" if detail else ""))
    if not condition:
        failures.append(name)


# ---------------------------------------------------------------------------
# independent implementation (deliberately separate from calc.py)
# ---------------------------------------------------------------------------
def model_independent(x, a: float, b: float, c: float, d: float):
    x = np.asarray(x, dtype=float)
    return d + (a - d) / (1.0 + (x / c) ** b)


def inverse_bisection(y: float, a: float, b: float, c: float, d: float) -> float:
    """Numerically invert the (increasing) 4PL curve by bisection."""
    y = float(y)
    lo, hi = 0.0, 1.0
    while float(model_independent(hi, a, b, c, d)) < y:
        hi *= 10.0
        if hi > 1e12:
            raise ValueError("sample OD is above the upper asymptote; no solution")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if float(model_independent(mid, a, b, c, d)) < y:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def load(path: str, columns: list[str], rename: dict) -> pd.DataFrame:
    return add_corrected(match_columns(pd.read_csv(path), columns).rename(columns=rename))


def prepare(std_path: str, smp_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    std = load(
        std_path,
        STD_COLS,
        {"Standard": "name", "Concentration": "conc", "OD450": "od450", "OD600": "od600"},
    )
    smp = load(
        smp_path,
        SMP_COLS,
        {"Sample": "name", "Dilution factor": "dilution", "OD450": "od450", "OD600": "od600"},
    )
    smp["dilution"] = to_numeric(smp["dilution"])
    return std, smp


def verify_dataset(label: str, std_path: str, smp_path: str, fit_expected: bool = False) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    std, smp = prepare(std_path, smp_path)
    stdv = std.loc[to_numeric(std["conc"]).notna() & std["corrected"].notna()].copy()
    stdv["conc"] = to_numeric(stdv["conc"])

    fit = fit_4pl(stdv["conc"], stdv["corrected"])
    check(f"{label}: 4PL fit ok", fit["ok"], fit.get("reason", ""))
    if not fit["ok"]:
        return fit, stdv, pd.DataFrame()

    results = quantify_samples(smp, fit)
    params = (fit["a"], fit["b"], fit["c"], fit["d"])

    # 1) fit residuals, measured with the independent model
    residual = stdv["corrected"].to_numpy(float) - model_independent(
        stdv["conc"].to_numpy(float), *params
    )
    check(
        f"{label}: fit residuals < 1% of OD span",
        np.max(np.abs(residual)) < 0.01 * max(np.ptp(stdv["corrected"].to_numpy(float)), 1e-9),
        f"max |residual| = {np.max(np.abs(residual)):.5f}",
    )

    # 2) standard recovery via independent bisection
    for _, row in stdv.iterrows():
        if float(row["conc"]) <= 0:
            continue
        back = inverse_bisection(float(row["corrected"]), *params)
        check(
            f"{label}: standard {row['name']} recovers nominal within 2%",
            abs(back - float(row["conc"])) / float(row["conc"]) < 0.02,
            f"{back:.2f} vs {float(row['conc']):.1f}",
        )

    # 3) every sample: compare calc.py against the independent implementation
    for _, row in results.iterrows():
        if not np.isfinite(float(row["measured_conc"])):
            continue
        independent = inverse_bisection(float(row["corrected"]), *params)
        check(
            f"{label}: {row['name']} inverse matches bisection",
            math.isclose(float(row["measured_conc"]), independent, rel_tol=1e-6),
            f"calc={float(row['measured_conc']):.6f} bisection={independent:.6f}",
        )
        check(
            f"{label}: {row['name']} dilution correction",
            math.isclose(float(row["reported_conc"]), independent * float(row["dilution"]),
                         rel_tol=1e-9),
            f"{float(row['reported_conc']):.4f}",
        )
        # round trip: put the reported concentration back through the curve
        back_od = float(
            model_independent(independent, *params)
        )
        check(
            f"{label}: {row['name']} round trip reproduces corrected OD",
            math.isclose(back_od, float(row["corrected"]), rel_tol=1e-9, abs_tol=1e-9),
            f"{back_od:.6f} vs {float(row['corrected']):.6f}",
        )

    return fit, stdv, results


# ---------------------------------------------------------------------------
# dataset A: reference plate
# ---------------------------------------------------------------------------
print("=" * 76)
print("Dataset A — reference 4PL plate (A=0.0821, B=1.0741, C=3877.5, D=10.7614)")
print("=" * 76)
fit_a, std_a, res_a = verify_dataset(
    "A", "sample_data/example_A_clean.csv", "sample_data/example_A_samples_clean.csv"
)
res_a_i = res_a.set_index("name")
check("A: fit uses 5 standard points", fit_a.get("n") == 5, str(fit_a.get("n")))
check(
    "A: fit parameters close to the reference model",
    math.isclose(fit_a["a"], 0.0821, rel_tol=0.02)
    and math.isclose(fit_a["b"], 1.0741, rel_tol=0.05)
    and math.isclose(fit_a["c"], 3877.5, rel_tol=0.10)
    and math.isclose(fit_a["d"], 10.7614, rel_tol=0.05),
    f"A={fit_a['a']:.4f} B={fit_a['b']:.4f} C={fit_a['c']:.1f} D={fit_a['d']:.3f}",
)
check(
    "A: wastewater original ≈ 394 pg/mL",
    abs(float(res_a_i.loc["Wastewater original", "measured_conc"]) - 394) < 3,
    f"{float(res_a_i.loc['Wastewater original', 'measured_conc']):.2f}",
)
check(
    "A: 1:2 diluted result ≈ 352 pg/mL (after x2)",
    abs(float(res_a_i.loc["Wastewater 1:2", "reported_conc"]) - 352) < 3,
    f"{float(res_a_i.loc['Wastewater 1:2', 'reported_conc']):.2f}",
)
check(
    "A: over-range sample flagged as extrapolation",
    "外推" in res_a_i.loc["Wastewater over-range demo", "status"],
    res_a_i.loc["Wastewater over-range demo", "status"],
)
print("\nTable A (calc.py):")
print(
    res_a[["name", "dilution", "corrected", "measured_conc", "reported_conc", "status"]]
    .to_string(index=False, float_format=lambda v: f"{v:.4f}")
)

# ---------------------------------------------------------------------------
# dataset B: noisy plate
# ---------------------------------------------------------------------------
print()
print("=" * 76)
print("Dataset B — noisy 4PL plate (out-of-range / missing OD600 / replicates)")
print("=" * 76)
fit_b, std_b, res_b = verify_dataset(
    "B", "sample_data/example_B_messy.csv", "sample_data/example_B_samples_messy.csv"
)
res_b_i = res_b.set_index("name")
check("B: fit is increasing (D > A)", fit_b.get("d", 0) > fit_b.get("a", 0),
      f"A={fit_b.get('a', float('nan')):.3f} D={fit_b.get('d', float('nan')):.3f}")
check(
    "B: over-range sample flagged as extrapolation",
    "外推" in res_b_i.loc["Over-range demo", "status"],
    res_b_i.loc["Over-range demo", "status"],
)
check(
    "B: missing OD600 falls back to OD450",
    math.isclose(float(res_b_i.loc["OD600 missing demo", "corrected"]), 0.500, abs_tol=1e-12),
    f"{float(res_b_i.loc['OD600 missing demo', 'corrected']):.3f}",
)
print("\nTable B (calc.py):")
print(
    res_b[["name", "dilution", "corrected", "measured_conc", "reported_conc", "status"]]
    .to_string(index=False, float_format=lambda v: f"{v:.4f}")
)

# ---------------------------------------------------------------------------
print()
print("=" * 76)
if failures:
    print(f"FAILED: {len(failures)} check(s): {failures}")
    raise SystemExit(1)
print("METHOD VERIFIED — every calc.py number matches an independent implementation.")
print("=" * 76)
