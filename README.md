# ELISA Result Generator (SCI 3102SEF · Lab 4)

A tiny Streamlit app where students type or paste **their own** ELISA readings
and get the standard curve, R², and the IFN-γ-like concentration in the
wastewater sample — no data or answers are hard-coded in the app.

Built for the worksheet line:

> “With use of AI tools, plot the Standard concentration vs Corrected
> absorbance curve, and determine the concentration of IFN-γ-like contaminants
> in wastewater sample.”

## Files

| File | Purpose |
|---|---|
| `app.py` | Streamlit page (2 input tabs + results tab) |
| `calc.py` | Pure calculation: corrected OD, linear fit, inverse + dilution |
| `test_calc.py` | Plain-Python checks for `calc.py` (no Streamlit needed) |
| `requirements.txt` | Dependencies |

## Run it locally

```bash
cd elisa-calculator
python3.13 -m venv .venv          # any Python 3.10–3.13 works
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

The browser opens at <http://localhost:8501>. First run only needs the
`pip install` step; afterwards just `source .venv/bin/activate && streamlit run app.py`.

### Share with a whole class on the lab Wi-Fi (no deployment)

```bash
streamlit run app.py --server.address 0.0.0.0 --server.port 8501
```

Then students open `http://<your-lan-ip>:8501` (find the IP with `ipconfig getifaddr en0`).
Each browser session gets its own data; nothing is stored on the server.

### Share with a permanent public link (no terminal needed)

This repo is already published at
**<https://github.com/lab-Kason/elisa-calculator>** (public, branch `main`).
The steps below are for reference / re-publishing.

**Step 1 — publish to GitHub (GitHub Desktop, ~2 min)**

1. Open GitHub Desktop → *Sign in to GitHub.com* if you never have.
2. `File` → `Add Local Repository…` → choose this `elisa-calculator`
   folder (the one containing `.git`) → **Add Repository**.
3. Click **Publish repository** in the top bar → Name:
   `elisa-calculator`, description
   “Streamlit app: ELISA standard curve + IFN-γ quantification (SCI 3102SEF Lab 4)”,
   **untick “Keep this code private”** → **Publish repository**.

**Step 2 — deploy (Streamlit Community Cloud, ~2 min)**

1. Open <https://share.streamlit.io> → *Sign in with GitHub*.
2. **Create app** → *Deploy a public app from GitHub*.
3. Repository: `lab-Kason/elisa-calculator` · Branch: `main` ·
   Main file path: `app.py`.
4. **Deploy** — wait 1–2 minutes; the public URL
   `https://<something>.streamlit.app` is ready to share on OLE.

Health check for the deployed app: `https://<something>.streamlit.app/~/+/healthz`
should answer `ok`. No secrets, no extra config needed.

## Student workflow

1. **① 标准品** — one row per standard well: name, concentration, OD450, OD600.
   Concentrations are free text: use whatever the group actually prepared
   (2000/1000/500/250 or something else). A blank can be added as a row with
   concentration `0`.
2. **② 样品** — one row per sample well: name, dilution factor (undiluted = 1),
   OD450, OD600.
3. **③ 结果** — standard curve + R², per-sample corrected absorbance,
   measured concentration, and concentration × dilution factor. Download CSV
   for the report.

Rows can be added/deleted in the tables (bottom row) and pasted straight from
Excel. Alternatively upload a CSV made from `Standards template` /
`Samples template` in the sidebar.

## Formulas used

```
Corrected absorbance = OD450 − OD600          (OD600 blank → 0)
Standard curve       : y = slope·x + intercept      (least squares)
Measured conc        = (corrected − intercept) / slope
Reported conc        = measured conc × dilution factor
```

Readings outside the standard OD range are flagged
**低于 / 高于标准曲线范围（外推）** — the app still calculates them, but the
report should mention the extrapolation.

## Deliberate limitations

- It is a **teaching aid**, not a validated assay tool: no blank subtraction
  beyond what students enter, no 4PL/5PL, no LOD/LOQ, no replicate averaging
  (enter each well as its own row if you want that).
- Units are labels only: type standards and results in the same unit
  (the unit selector is in the sidebar). No pg/mL ↔ ng/mL conversion.
- Anything the student enters is treated as-is; impossible data produces a
  warning, not an error.

## Verify / extend

```bash
python test_calc.py     # calculation + CSV header checks
```

To add a 4-parameter logistic fit later, add a function in `calc.py`, then a
fit-method radio in `app.py`; the rest of the pipeline stays the same.

---

## 中文速查（给 TA）

- 这是个**单文件 Streamlit 网页**：学生自己输 OD 读数 → 自动算校正吸光度 →
  线性拟合标准曲线 → 反算废水样品 IFN-γ 浓度并乘回稀释倍数。App 里没有任何
  预置的标准品浓度或预设答案。
- 本地跑：`cd elisa-calculator && python3.13 -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt && streamlit run app.py`
- 课堂同一 Wi-Fi 共用：`streamlit run app.py --server.address 0.0.0.0`，把
  `http://<你的IP>:8501` 发给学生即可，每台电脑数据独立。
- 要公开链接就推 GitHub，再在 share.streamlit.io 新建 app，Main file 选 `app.py`。
- 侧边栏有 CSV 模板下载；表格里可以直接从 Excel 粘贴、按底部加行。
- 已验证：`test_calc.py` 全通过；用 Streamlit 1.67 的 AppTest 跑过空数据、
  有数据、清空数据三个场景，均无报错。
- 想升级成 4PL 拟合或加重复孔平均，告诉我即可，改 `calc.py` 一个函数即可。
