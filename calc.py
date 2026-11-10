"""Core ELISA calculations for the Streamlit app.

This module deliberately has no Streamlit imports, so the maths can be tested
with plain Python (see test_calc.py).

Conventions
-----------
* Canonical column names used internally:
    standards: name, conc, od450, od600
    samples:   name, dilution, od450, od600
* Corrected absorbance = OD450 - OD600.  If OD600 is blank it is treated
  as 0, so the corrected value falls back to OD450.
* The standard curve is a straight line y = slope * x + intercept, where
  x = standard concentration and y = corrected absorbance.
* Nothing is hard-coded: the standard concentrations are whatever the user
  types.  Sample concentration = (absorbance - intercept) / slope, then
  multiplied by the sample's dilution factor.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

RESULT_COLUMNS = [
    "name",
    "dilution",
    "od450",
    "od600",
    "corrected",
    "measured_conc",
    "reported_conc",
    "status",
    "note",
]


def to_numeric(series: pd.Series) -> pd.Series:
    """Best-effort numeric conversion; unparsable entries become NaN."""
    return pd.to_numeric(series, errors="coerce")


def add_corrected(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with numeric OD columns and a `corrected` column."""
    out = df.copy()
    od450 = to_numeric(out["od450"]) if "od450" in out else pd.Series(dtype=float)
    od600 = to_numeric(out["od600"]) if "od600" in out else pd.Series(dtype=float)
    out["od450"] = od450
    out["od600"] = od600
    out["corrected"] = od450 - od600.fillna(0.0)
    return out


def valid_points(df: pd.DataFrame, x_col: str, y_col: str) -> pd.DataFrame:
    """Rows where both x and y are numeric."""
    x = to_numeric(df[x_col])
    y = to_numeric(df[y_col])
    return df.loc[x.notna() & y.notna()].copy()


def fit_linear(x: pd.Series, y: pd.Series) -> dict:
    """Least-squares straight line through the standard points.

    Returns a dict with either {"ok": False, "reason": ...} or the fit
    parameters and summary statistics.
    """
    x = to_numeric(pd.Series(x)).to_numpy(dtype=float)
    y = to_numeric(pd.Series(y)).to_numpy(dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = int(len(x))
    if n < 2:
        return {"ok": False, "reason": f"需要至少 2 个有效标准点（目前 {n} 个）", "n": n}
    if float(np.ptp(x)) == 0.0:
        return {"ok": False, "reason": "标准品浓度全部相同，无法拟合", "n": n}

    slope, intercept = np.polyfit(x, y, 1)
    slope = float(slope)
    intercept = float(intercept)
    fitted = slope * x + intercept
    ss_res = float(np.sum((y - fitted) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return {
        "ok": True,
        "n": n,
        "slope": slope,
        "intercept": intercept,
        "r2": r2,
        "x_min": float(x.min()),
        "x_max": float(x.max()),
        "y_min": float(y.min()),
        "y_max": float(y.max()),
    }


def sample_concentration(fit: dict, corrected: float) -> float:
    """Inverse of the standard curve: absorbance -> concentration."""
    if not fit.get("ok"):
        raise ValueError(fit.get("reason", "no fit"))
    slope = fit["slope"]
    if not math.isfinite(slope) or slope == 0:
        raise ValueError("斜率为 0，无法反算浓度")
    return (float(corrected) - fit["intercept"]) / slope


def quantify_samples(samples: pd.DataFrame, fit: dict) -> pd.DataFrame:
    """Turn sample absorbances into concentrations using the fitted curve.

    Deliberately permissive: missing values or out-of-range readings produce
    a note instead of an exception, because students may enter imperfect data.
    """
    rows: list[dict] = []
    for _, row in samples.iterrows():
        dilution_raw = to_numeric(pd.Series([row.get("dilution")])).iloc[0]
        note_parts: list[str] = []
        if pd.isna(dilution_raw) or float(dilution_raw) <= 0:
            dilution = 1.0
            note_parts.append("稀释倍数未填或无效，按 1 处理")
        else:
            dilution = float(dilution_raw)

        corrected = row.get("corrected")
        base = {
            "name": row.get("name"),
            "dilution": dilution,
            "od450": row.get("od450"),
            "od600": row.get("od600"),
            "corrected": corrected,
            "measured_conc": None,
            "reported_conc": None,
            "status": "无法计算",
            "note": "",
        }

        if not fit.get("ok"):
            base["note"] = fit.get("reason", "标准曲线尚未拟合")
            rows.append(base)
            continue
        if pd.isna(corrected):
            base["note"] = "缺少 OD450/OD600 数据"
            rows.append(base)
            continue

        try:
            measured = sample_concentration(fit, corrected)
        except ValueError as exc:  # slope is 0 / not finite
            base["note"] = str(exc)
            rows.append(base)
            continue

        reported = measured * dilution
        if corrected < fit["y_min"]:
            status = "低于标准曲线范围（外推）"
            note_parts.append("低于标准曲线最低点，结果仅供参考")
        elif corrected > fit["y_max"]:
            status = "高于标准曲线范围（外推）"
            note_parts.append("高于标准曲线最高点，建议加大稀释后重测")
        else:
            status = "在标准曲线范围内"
        if measured < 0:
            note_parts.append("反算浓度为负值，请检查 OD 或空白")

        base.update(
            measured_conc=measured,
            reported_conc=reported,
            status=status,
            note="；".join(note_parts),
        )
        rows.append(base)

    if not rows:
        return pd.DataFrame(columns=RESULT_COLUMNS)
    return pd.DataFrame(rows, columns=RESULT_COLUMNS)


def match_columns(df: pd.DataFrame, expected: list[str]) -> pd.DataFrame:
    """Make an uploaded CSV fit the expected headers, as forgiving as possible.

    Matching is case/whitespace insensitive; if the header count matches it
    falls back to positional renaming; otherwise missing columns are added
    as empty so the table still loads.
    """
    lookup = {str(c).strip().lower(): c for c in df.columns}
    if all(e.lower() in lookup for e in expected):
        return df.rename(columns={lookup[e.lower()]: e for e in expected})[expected]
    if df.shape[1] == len(expected):
        out = df.copy()
        out.columns = expected
        return out
    return df.reindex(columns=expected)
