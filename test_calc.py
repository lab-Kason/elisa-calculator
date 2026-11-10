"""Plain-python checks for calc.py.  Run: python test_calc.py"""

import numpy as np
import pandas as pd

from calc import add_corrected, fit_linear, quantify_samples


def check(name: str, condition: bool, detail: str = "") -> None:
    print(("PASS " if condition else "FAIL ") + name + (f"  ({detail})" if detail else ""))
    if not condition:
        raise SystemExit(1)


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

# --- exact synthetic curve: fit must recover slope/intercept ----------------
exact = add_corrected(
    pd.DataFrame(
        {
            "name": ["S1", "S2", "S3", "S4"],
            "conc": [2000.0, 1000.0, 500.0, 250.0],
            "od450": [0.0008 * c + 0.05 + 0.01 for c in (2000, 1000, 500, 250)],
            "od600": [0.01, 0.01, 0.01, 0.01],
        }
    )
)
fit_exact = fit_linear(exact["conc"], exact["corrected"])
check("exact linear data fits", fit_exact["ok"], str(fit_exact))
check("slope recovered", np.isclose(fit_exact["slope"], 0.0008, atol=1e-12), f"{fit_exact['slope']:.6g}")
check("intercept recovered", np.isclose(fit_exact["intercept"], 0.05, atol=1e-12), f"{fit_exact['intercept']:.6g}")
check("R2 = 1 for exact data", np.isclose(fit_exact["r2"], 1.0), f"{fit_exact['r2']:.6f}")

# --- sample inversion + dilution factor ------------------------------------
samples = add_corrected(
    pd.DataFrame(
        {
            "name": ["Sample 1 (1x)", "Sample 2 (1x)", "Sample 2 (2x)"],
            "dilution": [1.0, 1.0, 2.0],
            "od450": [0.45, 0.25, 0.33],
            "od600": [0.0, 0.0, 0.0],
        }
    )
)
res = quantify_samples(samples, fit_exact)
meas = res["measured_conc"].to_numpy(dtype=float)
reported = res["reported_conc"].to_numpy(dtype=float)
# (0.45 - 0.05) / 0.0008 = 500 ; (0.25 - 0.05) / 0.0008 = 250 ; (0.33 - 0.05) / 0.0008 = 350
check("measured concentration", np.allclose(meas, [500.0, 250.0, 350.0]), str(meas.round(2).tolist()))
check("dilution factor is multiplied back", np.allclose(reported, [500.0, 250.0, 700.0]), str(reported.round(2).tolist()))
check("in-range samples are flagged OK", (res["status"] == "在标准曲线范围内").all(), str(res["status"].tolist()))

# --- out-of-range flag ------------------------------------------------------
far = add_corrected(
    pd.DataFrame({"name": ["too high"], "dilution": [1.0], "od450": [2.5], "od600": [0.0]})
)
res_far = quantify_samples(far, fit_exact)
check("above-top absorbance flagged", "高于标准曲线范围" in res_far["status"].iloc[0], res_far["status"].iloc[0])

# --- permissive edge cases --------------------------------------------------
one_point = fit_linear([100.0], [0.4])
check("one point is not enough", not one_point["ok"], one_point["reason"])

same_x = fit_linear([100.0, 100.0], [0.4, 0.5])
check("identical concentrations are rejected", not same_x["ok"], same_x["reason"])

missing = add_corrected(
    pd.DataFrame({"name": ["a"], "dilution": [np.nan], "od450": [np.nan], "od600": [np.nan]})
)
res_missing = quantify_samples(missing, fit_exact)
check("missing OD does not raise", res_missing["status"].iloc[0] == "无法计算")
check("missing dilution defaults to 1", res_missing["dilution"].iloc[0] == 1.0)

empty = quantify_samples(pd.DataFrame(columns=["name", "dilution", "od450", "od600"]), fit_exact)
check("empty sample table returns empty results", empty.empty)




# --- CSV header handling ----------------------------------------------------
from calc import match_columns

templ = match_columns(
    pd.DataFrame(columns=["Standard", "Concentration", "OD450", "OD600"]),
    ["Standard", "Concentration", "OD450", "OD600"],
)
check("template CSV loads", list(templ.columns) == ["Standard", "Concentration", "OD450", "OD600"])

lower = match_columns(
    pd.DataFrame(columns=["standard", "concentration", "od450", "od600"]),
    ["Standard", "Concentration", "OD450", "OD600"],
)
check("headers are case-insensitive", list(lower.columns) == ["Standard", "Concentration", "OD450", "OD600"])

weird = match_columns(pd.DataFrame(columns=["a", "b", "c", "d"]), ["Standard", "Concentration", "OD450", "OD600"])
check("same column count is renamed positionally", list(weird.columns) == ["Standard", "Concentration", "OD450", "OD600"])

short = match_columns(pd.DataFrame(columns=["Standard", "OD450"]), ["Standard", "Concentration", "OD450", "OD600"])
check("missing columns are added", list(short.columns) == ["Standard", "Concentration", "OD450", "OD600"])

print("\nAll checks passed (incl. CSV headers).")
