"""Plain-python checks for calc.py.  Run: python test_calc.py"""

import math

import numpy as np
import pandas as pd

from calc import (
    add_corrected,
    fit_4pl,
    invert_4pl,
    match_columns,
    model_4pl,
    quantify_samples,
    to_numeric,
)

STD_COLS = ["Standard", "Concentration", "OD450", "OD600"]
SMP_COLS = ["Sample", "Dilution factor", "OD450", "OD600"]

REFERENCE = {"A": 0.0821, "B": 1.0741, "C": 3877.5, "D": 10.7614}


def check(name: str, condition: bool, detail: str = "") -> None:
    print(("PASS " if condition else "FAIL ") + name + (f"  ({detail})" if detail else ""))
    if not condition:
        raise SystemExit(1)


def load_standards(path: str) -> pd.DataFrame:
    df = match_columns(pd.read_csv(path), STD_COLS).rename(
        columns={"Standard": "name", "Concentration": "conc", "OD450": "od450", "OD600": "od600"}
    )
    return add_corrected(df)


def load_samples(path: str) -> pd.DataFrame:
    df = match_columns(pd.read_csv(path), SMP_COLS).rename(
        columns={"Sample": "name", "Dilution factor": "dilution", "OD450": "od450", "OD600": "od600"}
    )
    df["dilution"] = to_numeric(df["dilution"])
    return add_corrected(df)


# --- corrected absorbance ---------------------------------------------------
std = add_corrected(
    pd.DataFrame(
        {
            "name": ["S1", "S2", "S3", "S4"],
            "conc": [2000.0, 1000.0, 500.0, 250.0],
            "od450": [1.912, 1.024, 0.541, 0.301],
            "od600": [0.012, 0.024, 0.021, 0.001],
        }
    )
)
check(
    "corrected = OD450 - OD600",
    np.allclose(std["corrected"], [1.900, 1.000, 0.520, 0.300]),
    str(std["corrected"].tolist()),
)

no_ref = add_corrected(pd.DataFrame({"name": ["x"], "conc": [1.0], "od450": [0.5], "od600": [np.nan]}))
check("blank OD600 falls back to OD450", no_ref["corrected"].iloc[0] == 0.5)

# --- 4PL model sanity -------------------------------------------------------
check(
    "4PL at X=0 equals A",
    math.isclose(float(model_4pl(np.array([0.0]), **{k.lower(): v for k, v in REFERENCE.items()})[0]),
                 REFERENCE["A"], rel_tol=1e-12),
)
check(
    "4PL at X=C equals (A+D)/2",
    math.isclose(
        float(model_4pl(np.array([REFERENCE["C"]]), **{k.lower(): v for k, v in REFERENCE.items()})[0]),
        (REFERENCE["A"] + REFERENCE["D"]) / 2,
        rel_tol=1e-12,
    ),
)

# --- inverse 4PL against the reference slide ---------------------------------
inv = invert_4pl(
    np.array([0.925, 0.454]),
    REFERENCE["A"], REFERENCE["B"], REFERENCE["C"], REFERENCE["D"],
)
check("OD 0.925 -> ~393.5 pg/mL", abs(float(inv[0]) - 393.5) < 2.0, f"{float(inv[0]):.2f}")
check("OD 0.454 -> ~175.9 pg/mL", abs(float(inv[1]) - 175.9) < 2.0, f"{float(inv[1]):.2f}")
check("OD 0.454 at dilution 2 -> ~352 pg/mL", abs(float(inv[1]) * 2 - 352.0) < 3.0)

# --- dataset A: reference plate, end to end ---------------------------------
std_a = load_standards("sample_data/example_A_clean.csv")
smp_a = load_samples("sample_data/example_A_samples_clean.csv")
stdv_a = std_a.loc[to_numeric(std_a["conc"]).notna() & std_a["corrected"].notna()].copy()
stdv_a["conc"] = to_numeric(stdv_a["conc"])
fit_a = fit_4pl(stdv_a["conc"], stdv_a["corrected"])
check("A: 4PL fit succeeds", fit_a["ok"], fit_a.get("reason", ""))
check("A: uses 5 standard points", fit_a["n"] == 5, str(fit_a.get("n")))
check("A: fit is increasing (D > A)", fit_a["d"] > fit_a["a"], f"A={fit_a['a']:.3f} D={fit_a['d']:.2f}")
check(
    "A: one-degree-of-freedom warning is shown",
    any("自由度" in w for w in fit_a["warnings"]),
    str(fit_a["warnings"]),
)

res_a = quantify_samples(smp_a, fit_a).set_index("name")
check("A: wastewater original ~394 pg/mL",
      abs(float(res_a.loc["Wastewater original", "measured_conc"]) - 394) < 3,
      f"{float(res_a.loc['Wastewater original', 'measured_conc']):.1f}")
check("A: 1:2 dilution reports ~352 pg/mL after x2",
      abs(float(res_a.loc["Wastewater 1:2", "reported_conc"]) - 352) < 3,
      f"{float(res_a.loc['Wastewater 1:2', 'reported_conc']):.1f}")
check("A: over-range sample is flagged",
      "外推" in res_a.loc["Wastewater over-range demo", "status"],
      res_a.loc["Wastewater over-range demo", "status"])

# standard recovery through the fitted curve (blank excluded: its nominal is 0)
back = invert_4pl(stdv_a["corrected"].to_numpy(), fit_a["a"], fit_a["b"], fit_a["c"], fit_a["d"])
pos = stdv_a["conc"].to_numpy() > 0
recovery = back[pos] / stdv_a["conc"].to_numpy()[pos]
check("A: standard recovery 99-101%", np.all(np.abs(recovery - 1) < 0.01), str(np.round(recovery, 4)))

# --- dataset B: noisy plate -------------------------------------------------
std_b = load_standards("sample_data/example_B_messy.csv")
smp_b = load_samples("sample_data/example_B_samples_messy.csv")
stdv_b = std_b.loc[to_numeric(std_b["conc"]).notna() & std_b["corrected"].notna()].copy()
stdv_b["conc"] = to_numeric(stdv_b["conc"])
fit_b = fit_4pl(stdv_b["conc"], stdv_b["corrected"])
check("B: noisy 4PL fit succeeds", fit_b["ok"], fit_b.get("reason", ""))
res_b = quantify_samples(smp_b, fit_b).set_index("name")
check("B: OD600 missing falls back to OD450 (0.500)",
      math.isclose(float(res_b.loc["OD600 missing demo", "corrected"]), 0.500, abs_tol=1e-12),
      f"{float(res_b.loc['OD600 missing demo', 'corrected']):.3f}")
check("B: over-range sample is flagged",
      "外推" in res_b.loc["Over-range demo", "status"],
      res_b.loc["Over-range demo", "status"])
check("B: all other samples calculated",
      res_b.loc[["Wastewater original", "Wastewater 1:2", "Replicate A", "Replicate B"], "measured_conc"].notna().all())
check("B: replicates agree within 10%",
      abs(float(res_b.loc["Replicate A", "measured_conc"]) - float(res_b.loc["Replicate B", "measured_conc"]))
      / float(res_b.loc["Replicate A", "measured_conc"]) < 0.10)

# --- exactly 4 points (no blank) can still be fitted ------------------------
x4 = np.array([2000.0, 1000.0, 500.0, 250.0])
y4 = model_4pl(x4, REFERENCE["A"], REFERENCE["B"], REFERENCE["C"], REFERENCE["D"])
fit4 = fit_4pl(x4, y4)
check("4-point fit succeeds", fit4["ok"], fit4.get("reason", ""))
check("4-point fit has 0 residual degrees of freedom", fit4["dof"] == 0)
inv4 = invert_4pl(np.array([0.925, 0.454]), fit4["a"], fit4["b"], fit4["c"], fit4["d"])
check("4-point fit inverts reference ODs",
      abs(float(inv4[0]) - 393.5) < 2 and abs(float(inv4[1]) - 175.9) < 2,
      str(np.round(inv4, 1)))

# --- dilution factor and permissive edge cases ------------------------------
samples = add_corrected(
    pd.DataFrame(
        {
            "name": ["orig", "dil2"],
            "dilution": [1.0, 2.0],
            "od450": [0.925, 0.454],
            "od600": [0.0, 0.0],
        }
    )
)
res = quantify_samples(samples, fit_a).set_index("name")
check("dilution factor is multiplied back",
      math.isclose(float(res.loc["dil2", "reported_conc"]),
                   float(res.loc["dil2", "measured_conc"]) * 2.0, rel_tol=1e-12))

missing = add_corrected(
    pd.DataFrame({"name": ["a"], "dilution": [np.nan], "od450": [np.nan], "od600": [np.nan]})
)
res_missing = quantify_samples(missing, fit_a)
check("missing OD does not raise", res_missing["status"].iloc[0] == "无法计算")
check("missing dilution defaults to 1", res_missing["dilution"].iloc[0] == 1.0)

invalid = add_corrected(
    pd.DataFrame({"name": ["neg"], "dilution": [0.0], "od450": [0.5], "od600": [0.0]})
)
res_invalid = quantify_samples(invalid, fit_a)
check("invalid dilution defaults to 1 with a note",
      res_invalid["dilution"].iloc[0] == 1.0 and "稀释倍数" in res_invalid["note"].iloc[0])

# --- rejected fits ----------------------------------------------------------
too_few = fit_4pl([2000.0, 1000.0, 500.0], [3.0, 2.0, 1.0])
check("three points are not enough", not too_few["ok"], too_few.get("reason", ""))
same_x = fit_4pl([1000.0, 1000.0, 1000.0, 1000.0], [0.1, 0.2, 0.3, 0.4])
check("identical concentrations rejected", not same_x["ok"], same_x.get("reason", ""))
flat_y = fit_4pl([0.0, 250.0, 500.0, 1000.0, 2000.0], [0.1, 0.1, 0.1, 0.1, 0.1])
check("identical OD values rejected", not flat_y["ok"], flat_y.get("reason", ""))

# --- CSV header handling ----------------------------------------------------
templ = match_columns(
    pd.DataFrame(columns=["Standard", "Concentration", "OD450", "OD600"]),
    STD_COLS,
)
check("template CSV loads", list(templ.columns) == STD_COLS)
lower = match_columns(
    pd.DataFrame(columns=["standard", "concentration", "od450", "od600"]),
    STD_COLS,
)
check("headers are case-insensitive", list(lower.columns) == STD_COLS)
weird = match_columns(pd.DataFrame(columns=["a", "b", "c", "d"]), STD_COLS)
check("same column count is renamed positionally", list(weird.columns) == STD_COLS)
short = match_columns(pd.DataFrame(columns=["Standard", "OD450"]), STD_COLS)
check("missing columns are added", list(short.columns) == STD_COLS)

print("\nAll checks passed.")
