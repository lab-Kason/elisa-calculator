"""Core ELISA calculations for the Streamlit app.

This module deliberately has no Streamlit imports, so the maths can be tested
with plain Python (see test_calc.py).

Model
-----
The standard curve is a 4-parameter logistic (4PL) function

    Y = D + (A - D) / (1 + (X / C) ** B)

    A = lower asymptote ("zero" response)
    B = Hill slope
    C = inflection point (EC50), in the same concentration unit as X
    D = upper asymptote ("max" response)

Conventions
-----------
* Canonical column names used internally:
    standards: name, conc, od450, od600
    samples:   name, dilution, od450, od600
* Corrected absorbance = OD450 - OD600.  If OD600 is blank it is treated
  as 0, so the corrected value falls back to OD450.
* The blank well (concentration 0), when entered, is used as the zero
  standard and included in the fit, as in the reference 4PL analysis.
* At least 4 valid standard points are required to fit 4 parameters; 5 points
  (blank + S1-S4) are recommended.
* Sample concentration = inverse 4PL of the corrected absorbance, then
  multiplied by the sample's dilution factor.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

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


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# 4PL model
# ---------------------------------------------------------------------------
def model_4pl(x, a: float, b: float, c: float, d: float):
    """Evaluate Y = D + (A - D) / (1 + (X / C) ** B)."""
    x = np.asarray(x, dtype=float)
    a, b, c, d = float(a), float(b), float(c), float(d)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        base = np.power(x / c, b)
        y = d + (a - d) / (1.0 + base)
    return y


def invert_4pl(y, a: float, b: float, c: float, d: float) -> np.ndarray:
    """Back-calculate X from Y.

    Rearranging Y = D + (A - D) / (1 + (X/C)^B) gives

        X = C * ((A - Y) / (Y - D)) ** (1 / B)

    Only defined when Y is strictly between the asymptotes A and D;
    returns NaN for out-of-range Y.
    """
    y = np.asarray(y, dtype=float)
    a, b, c, d = float(a), float(b), float(c), float(d)
    num = a - y
    den = y - d
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        ratio = np.where(den != 0.0, num / den, np.nan)
        ratio = np.where(ratio >= 0.0, ratio, np.nan)
        x = c * np.power(ratio, 1.0 / b)
    return x


def _unpack(q):
    """Unconstrained parameters -> (a, b, c, d) with d > a, b > 0, c > 0."""
    a, k, log_c, log_b = q
    d = a + np.exp(k)
    c = np.exp(log_c)
    b = np.exp(log_b)
    return a, b, c, d


def _residual(q, x, y, weighted):
    a, b, c, d = _unpack(q)
    y_hat = model_4pl(x, a, b, c, d)
    r = y - y_hat
    if weighted:
        # weight = 1 / y^2 (the convention used on the reference slide)
        r = r / np.maximum(np.abs(y_hat), 1e-6)
    return r


def fit_4pl(x, y, weighted: bool = False) -> dict:
    """Fit the 4PL model by (weighted) non-linear least squares.

    Returns a dict.  On success: ok=True plus the parameters A/B/C/D, R2/RMSE,
    n, dof, the standard x/y range, the equation string and warnings.  On
    failure: ok=False with a human-readable `reason`.
    """
    x = to_numeric(pd.Series(x)).to_numpy(dtype=float)
    y = to_numeric(pd.Series(y)).to_numpy(dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = int(x.size)

    if n < 4:
        return {
            "ok": False,
            "reason": f"需要至少 4 个有效标准点（目前 {n} 个）。"
                      "建议 5 个：空白孔（浓度 0）+ S1–S4。",
            "n": n,
        }
    if np.any(x < 0):
        return {"ok": False, "reason": "标准品浓度不能为负数，请检查浓度列。", "n": n}
    if float(np.ptp(x)) == 0.0:
        return {"ok": False, "reason": "标准品浓度全部相同，无法拟合 4PL。", "n": n}
    x_pos = x[x > 0]
    if x_pos.size < 2 or np.unique(x_pos).size < 2:
        return {
            "ok": False,
            "reason": "至少需要两个不同的非零标准品浓度（S1–S4）才能拟合 4PL。",
            "n": n,
        }
    if np.unique(y).size < 3:
        return {
            "ok": False,
            "reason": "标准品 OD 读数几乎相同，无法拟合曲线。请检查标准品孔与读数。",
            "n": n,
        }

    span = float(max(np.ptp(y), 1e-6))
    y_min = float(y.min())
    y_max = float(y.max())

    # Starting points
    a0 = y_min - 0.05 * span
    d0 = y_max + 0.05 * span
    c_guesses = {
        float(np.exp(np.mean(np.log(x_pos)))),
        float(np.median(x_pos)),
        float(x_pos.min() * 1.5),
        float(x_pos.max() / 1.5),
    }
    c_guesses = sorted(v for v in c_guesses if np.isfinite(v) and v > 0)
    b_guesses = (0.5, 1.0, 2.0, 4.0)

    lower = np.array(
        [a0 - 20.0 * span, -30.0, np.log(max(1e-6, x_pos.min() / 1e4)), np.log(1e-3)]
    )
    upper = np.array([d0 + 20.0 * span, 30.0, np.log(1e10), np.log(1e2)])

    best_q = None
    best_sse = np.inf
    for c0 in c_guesses:
        for b0 in b_guesses:
            q0 = np.array([a0, np.log(max(d0 - a0, 1e-9)), np.log(c0), np.log(b0)])
            try:
                sol = least_squares(
                    _residual,
                    q0,
                    args=(x, y, weighted),
                    bounds=(lower, upper),
                    method="trf",
                    x_scale="jac",
                    max_nfev=20000,
                    ftol=1e-12,
                    xtol=1e-12,
                    gtol=1e-12,
                )
            except (ValueError, FloatingPointError):
                continue
            if not np.all(np.isfinite(sol.x)):
                continue
            a, b, c, d = _unpack(sol.x)
            y_hat = model_4pl(x, a, b, c, d)
            sse = float(np.sum((y - y_hat) ** 2))
            if np.isfinite(sse) and sse < best_sse:
                best_sse = sse
                best_q = sol.x

    if best_q is None:
        return {
            "ok": False,
            "reason": "4PL 拟合不收敛。请检查标准品浓度与 OD 读数是否一一对应，"
                      "并确认 OD 随浓度递增。",
            "n": n,
        }

    a, b, c, d = _unpack(best_q)
    y_hat = model_4pl(x, a, b, c, d)
    sse = float(np.sum((y - y_hat) ** 2))
    sst = float(np.sum((y - y.mean()) ** 2))
    r2 = float(1.0 - sse / sst) if sst > 0 else float("nan")
    rmse = float(np.sqrt(sse / n))
    dof = int(n - 4)

    warnings: list[str] = []
    order = np.argsort(x)
    if np.any(np.diff(y[order]) <= 0):
        warnings.append(
            "标准品 OD 没有随浓度稳定上升（应大致 S4 < S3 < S2 < S1）。"
            "请检查标准品顺序和读数是否填错。"
        )
    if dof <= 1:
        warnings.append(
            f"只有 {n} 个标准点拟合 4 个参数（{dof} 个自由度）：曲线几乎/完全穿过所有点，"
            "R² 与置信区间没有参考价值。建议加入空白孔（浓度 0）或更多标准点。"
        )
    if not (a < d):
        warnings.append("拟合得到的下平台 A 不低于上平台 D，曲线可能不稳定，请检查数据。")

    equation = (
        f"Y = {d:.4f} + ({a:.4f} - {d:.4f}) / (1 + (X / {c:,.1f}) ^ {b:.4f})"
    )

    return {
        "ok": True,
        "n": n,
        "a": float(a),
        "b": float(b),
        "c": float(c),
        "d": float(d),
        "sse": sse,
        "sst": sst,
        "r2": r2,
        "rmse": rmse,
        "dof": dof,
        "x_min": float(np.min(x)),
        "x_max": float(np.max(x)),
        "y_min": float(np.min(y)),
        "y_max": float(np.max(y)),
        "equation": equation,
        "warnings": warnings,
        "weighted": bool(weighted),
    }


# ---------------------------------------------------------------------------
# Samples
# ---------------------------------------------------------------------------
def sample_concentration(fit: dict, corrected: float) -> float:
    """Inverse of the 4PL standard curve: corrected absorbance -> concentration."""
    if not fit.get("ok"):
        raise ValueError(fit.get("reason", "标准曲线尚未拟合"))
    if not math.isfinite(float(corrected)):
        raise ValueError("缺少校正吸光度")
    x = float(invert_4pl(np.array([corrected]), fit["a"], fit["b"], fit["c"], fit["d"])[0])
    if not math.isfinite(x):
        raise ValueError("校正吸光度落在 4PL 曲线之外（无解），请检查 OD 或空白")
    return x


def quantify_samples(samples: pd.DataFrame, fit: dict) -> pd.DataFrame:
    """Turn sample absorbances into concentrations using the fitted 4PL curve.

    Deliberately permissive: missing values or out-of-range readings produce a
    note instead of an exception, because students may enter imperfect data.
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
            measured = sample_concentration(fit, float(corrected))
        except ValueError as exc:
            base["note"] = str(exc)
            rows.append(base)
            continue

        reported = measured * dilution
        if measured > fit["x_max"]:
            status = "高于标准曲线范围（外推）"
            note_parts.append("高于最高标准点，建议加大稀释后重测")
        elif measured < fit["x_min"]:
            status = "低于标准曲线范围（外推）"
            note_parts.append("低于最低标准点，结果仅供参考")
        else:
            status = "在标准曲线范围内"
        if reported < 0:
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


# ---------------------------------------------------------------------------
# CSV uploads
# ---------------------------------------------------------------------------
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
