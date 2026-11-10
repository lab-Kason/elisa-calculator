"""ELISA Result Generator — Streamlit front-end.

Students type or paste their own plate-reader readings; the app builds the
standard curve and reports the IFN-γ-like concentration in the wastewater
sample.  No data, concentrations or results are hard-coded.

Run locally:
    streamlit run app.py

Deploy:
    push this folder to GitHub, then create the app on
    https://share.streamlit.io (main file: app.py).
"""

from __future__ import annotations

from datetime import datetime

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st

from calc import add_corrected, fit_linear, match_columns, quantify_samples, to_numeric

# Column headers used by the two tables and by the CSV templates.
STANDARD_COLUMNS = ["Standard", "Concentration", "OD450", "OD600"]
SAMPLE_COLUMNS = ["Sample", "Dilution factor", "OD450", "OD600"]

st.set_page_config(page_title="ELISA Result Generator", page_icon="🧪", layout="wide")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def default_standards() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Standard": ["S1", "S2", "S3", "S4"],
            "Concentration": [np.nan] * 4,
            "OD450": [np.nan] * 4,
            "OD600": [np.nan] * 4,
        }
    )


def default_samples() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Sample": ["Sample 1", "Sample 2"],
            "Dilution factor": [np.nan, np.nan],
            "OD450": [np.nan, np.nan],
            "OD600": [np.nan, np.nan],
        }
    )


def load_uploaded(uploaded, prefix: str, expected: list[str]) -> None:
    """Adopt a CSV once per new file, keeping table edits from other files."""
    if uploaded is None:
        return
    signature = (uploaded.name, uploaded.size)
    if st.session_state.get(f"{prefix}_sig") == signature:
        return
    try:
        df = pd.read_csv(uploaded)
    except Exception as exc:  # noqa: BLE001 - show the real error to the user
        st.error(f"无法读取 {uploaded.name}：{exc}")
        return
    st.session_state[f"{prefix}_data"] = match_columns(df, expected)
    st.session_state[f"{prefix}_sig"] = signature
    st.session_state[f"{prefix}_ver"] = st.session_state.get(f"{prefix}_ver", 0) + 1
    st.success(f"已载入 {uploaded.name} 中的 {len(df)} 行，可在下表继续修改。")


def clear_all() -> None:
    for prefix in ("std", "smp"):
        st.session_state.pop(f"{prefix}_data", None)
        st.session_state.pop(f"{prefix}_sig", None)
        st.session_state[f"{prefix}_ver"] = st.session_state.get(f"{prefix}_ver", 0) + 1


def editor_base(prefix: str, default_factory) -> pd.DataFrame:
    base = st.session_state.get(f"{prefix}_data")
    if base is None:
        base = default_factory()
    return base


def canonical_standards(df: pd.DataFrame) -> pd.DataFrame:
    out = df.rename(
        columns={
            "Standard": "name",
            "Concentration": "conc",
            "OD450": "od450",
            "OD600": "od600",
        }
    ).copy()
    for col in ("name", "conc", "od450", "od600"):
        if col not in out:
            out[col] = np.nan
    out["name"] = out["name"].astype("string").fillna("")
    return add_corrected(out)


def canonical_samples(df: pd.DataFrame) -> pd.DataFrame:
    out = df.rename(
        columns={
            "Sample": "name",
            "Dilution factor": "dilution",
            "OD450": "od450",
            "OD600": "od600",
        }
    ).copy()
    for col in ("name", "dilution", "od450", "od600"):
        if col not in out:
            out[col] = np.nan
    out["name"] = out["name"].astype("string").fillna("")
    out["dilution"] = to_numeric(out["dilution"])
    return add_corrected(out)


def fmt(value, digits: int = 3) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    if pd.isna(value):
        return ""
    return f"{float(value):.{digits}f}"


# ---------------------------------------------------------------------------
# sidebar
# ---------------------------------------------------------------------------
with st.sidebar:
    st.header("设置 Settings")
    unit = st.selectbox(
        "浓度单位 Concentration unit",
        ["pg/mL", "ng/mL", "mg/mL"],
        index=0,
        help="只影响显示与下载文件里的单位标注，不做换算：直接在「标准品」里用同一单位输入即可。",
    )
    st.divider()
    st.subheader("CSV 模板 Templates")
    st.download_button(
        "⬇️ Standards template",
        pd.DataFrame(columns=STANDARD_COLUMNS).to_csv(index=False),
        file_name="standards_template.csv",
        mime="text/csv",
        width="stretch",
    )
    st.download_button(
        "⬇️ Samples template",
        pd.DataFrame(columns=SAMPLE_COLUMNS).to_csv(index=False),
        file_name="samples_template.csv",
        mime="text/csv",
        width="stretch",
    )
    st.divider()
    if st.button("🧹 清空所有数据 Clear all", width="stretch"):
        clear_all()
        st.rerun()

# ---------------------------------------------------------------------------
# page header
# ---------------------------------------------------------------------------
st.title("🧪 ELISA Result Generator")
st.caption(
    "SCI 3102SEF DNA & Biochemical Analysis Laboratory Training II — Lab 4: "
    "detection of IFN-γ-like contaminants in wastewater by ELISA.  "
    "Fill in your own readings; the standard curve and sample concentrations "
    "are calculated automatically."
)

# ---------------------------------------------------------------------------
# input tabs
# ---------------------------------------------------------------------------
tab_std, tab_smp, tab_result = st.tabs(
    ["① 标准品 Standards", "② 样品 Samples", "③ 结果 Results & plot"]
)

with tab_std:
    st.subheader("Standard wells")
    st.caption(
        "每行一个标准品孔：填入浓度和 OD 读数。浓度不需要是 2000/1000/500/250，"
        "按你实际配制的浓度输入即可（4 个点以上也可以）。若量了空白孔，加一行、浓度填 0。"
    )
    up_std = st.file_uploader(
        "或上传标准的 CSV（表头：Standard, Concentration, OD450, OD600）",
        type=["csv"],
        key="std_upload",
    )
    load_uploaded(up_std, "std", STANDARD_COLUMNS)
    std_ver = st.session_state.get("std_ver", 0)
    std_raw = st.data_editor(
        editor_base("std", default_standards),
        key=f"std_editor_{std_ver}",
        num_rows="dynamic",
        hide_index=True,
        width="stretch",
        column_config={
            "Standard": st.column_config.TextColumn("Standard", help="例如 S1、S2、Blank"),
            "Concentration": st.column_config.NumberColumn(
                f"Concentration ({unit})", format="%.2f", min_value=0.0
            ),
            "OD450": st.column_config.NumberColumn("OD450", format="%.3f"),
            "OD600": st.column_config.NumberColumn("OD600", format="%.3f"),
        },
    )

with tab_smp:
    st.subheader("Wastewater sample wells")
    st.caption(
        "每行一个样品孔：填样品名、稀释倍数（原样填 1）和 OD 读数。"
        "同一浓度测了多孔就多填几行，例如 Sample 1-1 / Sample 1-2。"
    )
    up_smp = st.file_uploader(
        "或上传样品的 CSV（表头：Sample, Dilution factor, OD450, OD600）",
        type=["csv"],
        key="smp_upload",
    )
    load_uploaded(up_smp, "smp", SAMPLE_COLUMNS)
    smp_ver = st.session_state.get("smp_ver", 0)
    smp_raw = st.data_editor(
        editor_base("smp", default_samples),
        key=f"smp_editor_{smp_ver}",
        num_rows="dynamic",
        hide_index=True,
        width="stretch",
        column_config={
            "Sample": st.column_config.TextColumn("Sample", help="例如 Sample 1 (1x)"),
            "Dilution factor": st.column_config.NumberColumn(
                "Dilution factor", format="%.2f", min_value=0.0,
                help="原样 = 1，2 倍稀释 = 2，10 倍稀释 = 10，以此类推",
            ),
            "OD450": st.column_config.NumberColumn("OD450", format="%.3f"),
            "OD600": st.column_config.NumberColumn("OD600", format="%.3f"),
        },
    )

# ---------------------------------------------------------------------------
# calculations
# ---------------------------------------------------------------------------
std = canonical_standards(std_raw)
smp = canonical_samples(smp_raw)

std_valid = std.loc[to_numeric(std["conc"]).notna() & std["corrected"].notna()].copy()
std_valid["conc"] = to_numeric(std_valid["conc"])
fit = fit_linear(std_valid["conc"], std_valid["corrected"])
results = quantify_samples(smp, fit)

# ---------------------------------------------------------------------------
# results tab
# ---------------------------------------------------------------------------
with tab_result:
    if not fit["ok"]:
        st.info(
            "👈 先在「① 标准品」页填入至少 2 个标准品孔（浓度 + OD450 或 OD450/OD600），"
            "这里会自动出现标准曲线和样品浓度。"
        )
        if fit.get("reason"):
            st.caption(f"当前有效标准点：{fit.get('n', 0)} 个。{fit['reason']}")
    else:
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("R²", f"{fit['r2']:.5f}" if np.isfinite(fit["r2"]) else "—")
        m2.metric("斜率 Slope", f"{fit['slope']:.4f}")
        m3.metric("截距 Intercept", f"{fit['intercept']:.4f}")
        m4.metric("标准点数量", f"{fit['n']}")
        st.caption(
            f"拟合方程：校正吸光度 = {fit['slope']:.4f} × 浓度 ({unit}) "
            f"{fit['intercept']:+.4f}　|　"
            f"标准曲线 OD 范围：{fit['y_min']:.3f} – {fit['y_max']:.3f}"
        )
        if fit["slope"] <= 0:
            st.warning("拟合斜率为 0 或负值，请检查 OD 读数与浓度是否配对正确。")
        elif np.isfinite(fit["r2"]) and fit["r2"] < 0.95:
            st.warning("R² 偏低，建议检查标准品数据是否有录入错误。")

        grid = np.linspace(fit["x_min"], fit["x_max"], 100)
        line_df = pd.DataFrame(
            {
                "Concentration": grid,
                "Corrected absorbance": fit["slope"] * grid + fit["intercept"],
            }
        )
        points_df = pd.DataFrame(
            {
                "Concentration": std_valid["conc"].to_numpy(dtype=float),
                "Corrected absorbance": std_valid["corrected"].to_numpy(dtype=float),
                "Standard": std_valid["name"].astype(str).to_numpy(),
            }
        )
        x_enc = alt.X("Concentration:Q", title=f"Standard concentration ({unit})")
        y_enc = alt.Y(
            "Corrected absorbance:Q",
            title="Corrected absorbance (OD450 − OD600)",
            scale=alt.Scale(zero=False),
        )
        chart = (
            alt.Chart(line_df)
            .mark_line(color="#c62828", strokeWidth=2.5)
            .encode(x=x_enc, y=y_enc)
            + alt.Chart(points_df)
            .mark_circle(size=120, color="#1565c0")
            .encode(
                x=x_enc,
                y=y_enc,
                tooltip=[
                    alt.Tooltip("Standard:N", title="Standard"),
                    alt.Tooltip("Concentration:Q", title="Concentration"),
                    alt.Tooltip("Corrected absorbance:Q", title="Corrected OD", format=".3f"),
                ],
            )
        )
        st.altair_chart(chart.properties(height=420).interactive(), width="stretch")
        st.caption("🔵 标准品实测点　🔴 线性拟合线　（Corrected absorbance = OD450 − OD600；OD600 留空则等于 OD450）")

        st.divider()
        st.subheader("样品浓度 Sample results")
        if results.empty:
            st.info("还没有样品数据，去「② 样品」页填入 OD 读数。")
        else:
            display = pd.DataFrame(
                {
                    "Sample": results["name"],
                    "Dilution factor": results["dilution"].map(lambda v: f"{v:g}"),
                    "Corrected absorbance": results["corrected"].map(fmt),
                    f"Measured conc ({unit})": results["measured_conc"].map(fmt),
                    f"Reported conc ({unit})": results["reported_conc"].map(fmt),
                    "Status": results["status"],
                    "Note": results["note"],
                }
            )
            st.dataframe(display, hide_index=True, width="stretch")
            st.caption(
                "Measured conc 是把该孔吸光度直接代回标准曲线得到的浓度；"
                "Reported conc = Measured conc × 稀释倍数。"
            )
            if (results["status"].str.contains("外推")).any():
                st.warning("部分样品的吸光度落在标准曲线范围之外，结果属外推，建议调整稀释倍数重测。")

            export = results.rename(
                columns={
                    "name": "Sample",
                    "dilution": "Dilution factor",
                    "corrected": "Corrected absorbance",
                    "measured_conc": f"Measured conc ({unit})",
                    "reported_conc": f"Reported conc ({unit})",
                    "status": "Status",
                    "note": "Note",
                }
            )
            st.download_button(
                "⬇️ 下载结果 CSV",
                export.to_csv(index=False),
                file_name=f"elisa_results_{datetime.now():%Y%m%d_%H%M}.csv",
                mime="text/csv",
            )

        with st.expander("计算说明 How this is calculated"):
            st.markdown(
                f"""
1. **校正吸光度**：Corrected absorbance = OD450 − OD600（OD600 留空时按 OD450 处理）。
2. **标准曲线**：对标准品孔做最小二乘直线拟合，
   校正吸光度 = 斜率 × 浓度 + 截距，浓度用你输入的值（不预设 2000/1000/500/250）。
3. **样品浓度**：Measured conc = (校正吸光度 − 截距) ÷ 斜率；
   Reported conc = Measured conc × 稀释倍数。
4. 只有落进标准曲线 OD 范围（{fit['y_min']:.3f} – {fit['y_max']:.3f}）内的读数才算可靠；
   范围外会标注「外推」，报告里应说明。
                """
            )

st.divider()
st.caption(
    "本工具只做算术与作图，不对实验数据的正确性负责；最终解释以实验手册与课程要求为准。"
)
