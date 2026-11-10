# ELISA 4PL Result Generator (SCI 3102SEF · Lab 4)

A tiny Streamlit app where students type or paste **their own** ELISA readings
and get the **4-parameter logistic (4PL)** standard curve, the fitted
parameters, and the IFN-γ-like concentration in the wastewater sample —
no data or answers are hard-coded in the app.

> **What changed vs the first version:** the earlier release fitted a straight
> line (`y = slope·x + intercept`). The course expects a **4PL** curve
> (`Y = D + (A − D) / (1 + (X/C)^B)`), so `calc.py` and the results tab were
> rewritten. The reference plate in `sample_data/` now reproduces the slide
> answers: **≈394 pg/mL** (well 6, undiluted) and **≈352 pg/mL** (well 7,
> measured at 1:2 and multiplied by 2).

Built for the worksheet line:

> “With use of AI tools, plot the Standard concentration vs Corrected
> absorbance curve, and determine the concentration of IFN-γ-like contaminants
> in wastewater sample.”

## Files

| File | Purpose |
|---|---|
| `app.py` | Streamlit page (2 input tabs + results tab) |
| `calc.py` | Pure calculation: corrected OD, 4PL fit, inverse + dilution |
| `test_calc.py` | Plain-Python checks for `calc.py` (no Streamlit needed) |
| `verify_method.py` | Re-derives every number with an independent bisection implementation |
| `test_app.py` | Headless UI smoke test (Streamlit `AppTest`) |
| `sample_data/` | Two example plates for testing; expected answers in `sample_data/README.md` |
| `requirements.txt` | Dependencies (note: **scipy** is now required for the 4PL fit) |

## Run it locally

```bash
cd elisa-calculator
python3.13 -m venv .venv          # any Python 3.10–3.13 works
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py
```

On Windows (PowerShell) without a venv:

```powershell
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

The browser opens at <http://localhost:8501>. First run only needs the
`pip install` step; afterwards just activate the venv and run streamlit again.

### Share with a whole class on the lab Wi-Fi (no deployment)

```bash
streamlit run app.py --server.address 0.0.0.0 --server.port 8501
```

Then students open `http://<your-lan-ip>:8501` (find the IP with `ipconfig` on
Windows). Each browser session gets its own data; nothing is stored on the
server.

### Share with a permanent public link (no terminal needed)

This repo is published at **<https://github.com/lab-Kason/elisa-calculator>**
(public, branch `main`), and the Streamlit app is already deployed on
Streamlit Community Cloud. Re-deploy after a push:

1. Push the updated files to the repo (Streamlit Cloud auto-redeploys; if the
   app was not created from this repo, go to <https://share.streamlit.io> →
   **New app** → repo `lab-Kason/elisa-calculator` → branch `main` →
   main file `app.py`).
2. Wait 1–2 minutes for the rebuild (adding **scipy** changes the environment,
   so the first rebuild can take a little longer).
3. Hard-refresh the app page; the title should read
   **“ELISA 4PL Result Generator”**.

Health check for the deployed app: `https://<something>.streamlit.app/~/+/healthz`
should answer `ok`. No secrets, no extra config needed.

## Student workflow

1. **① 标准品 Standards** — one row per standard well: name, concentration,
   OD450, OD600. At least **4 valid points** are required; the recommended
   layout is **Blank (concentration 0) + S1–S4**. Concentrations are free
   text: use whatever the group actually prepared, not necessarily
   2000/1000/500/250. A blank row is pre-filled as “Blank / 0”.
2. **② 样品 Samples** — one row per sample well: name, dilution factor
   (undiluted = 1), OD450, OD600.
3. **③ 结果 Results** — 4PL curve + fitted A/B/C/D, R²/RMSE, per-sample
   corrected absorbance, measured concentration (4PL inverse) and
   concentration × dilution factor. Download the CSV for the report.

Rows can be added/deleted in the tables (bottom row) and pasted straight from
Excel. Alternatively upload a CSV made from the *Standards template* /
*Samples template* in the sidebar.

## Formulas used

```
Corrected absorbance = OD450 − OD600      (OD600 blank → 0)

4PL standard curve   : Y = D + (A − D) / (1 + (X / C)^B)

    A = lower asymptote ("zero" response)
    B = Hill slope
    C = inflection point / EC50
    D = upper asymptote ("max" response)
    Y = corrected absorbance, X = concentration

Measured conc        = C × ((A − Y) / (Y − D))^(1/B)      (inverse 4PL)
Reported conc        = measured conc × dilution factor
```

Notes:

- The blank well (concentration 0), when entered, is used as the **zero
  standard** and included in the fit — this is the convention used in the
  reference 4PL analysis.
- The fit is **unweighted** non-linear least squares. `calc.fit_4pl()` also
  supports `weighted=True` (1/Y²) if the instructor prefers that convention;
  the UI currently uses the default.
- 4 PL has 4 parameters. With only 5 calibration points (blank + S1–S4) there
  is just **1 degree of freedom**, so R² and confidence intervals are weak
  indicators — the app warns about this.
- Readings with a concentration above the highest standard (or below the
  lowest, when no blank row is entered) are flagged
  **低于 / 高于标准曲线范围（外推）** — the app still calculates them, but the
  report should mention the extrapolation.

## Sample data and method verification

Two example plates live in `sample_data/` — see `sample_data/README.md` for the
expected answers and what each one is designed to exercise:

| Dataset | Purpose | Expected result |
|---|---|---|
| A — clean | plate generated from the reference 4PL model (`A=0.0821, B=1.0741, C=3877.5, D=10.7614`) | ≈394 pg/mL (original), ≈352 pg/mL (1:2 × 2); over-range demo flagged |
| B — messy | noisy curve, out-of-range well, blank OD600, replicates | results consistent with the fitted curve; over-range flagged |

`python verify_method.py` is the proof: it recomputes each inverse
concentration by **bisection** on an independently written 4PL equation
(different algorithm from `calc.py`), checks the dilution correction, does a
round-trip (concentration → OD), and verifies the reference answers above. It
must finish with **METHOD VERIFIED**.

## Deliberate limitations

- It is a **teaching aid**, not a validated assay tool: no blank subtraction
  beyond using the blank as the zero standard, no LOD/LOQ, no replicate
  averaging (enter each well as its own row if you want that), no weighted fit
  in the UI.
- Units are labels only: type standards and results in the same unit
  (the unit selector is in the sidebar). No pg/mL ↔ ng/mL conversion.
- Anything the student enters is treated as-is; impossible data produces a
  warning or a “无法计算” note, not a crash.

## Verify / extend

```bash
python test_calc.py      # calculation + CSV header checks
python verify_method.py  # independent re-derivation of every number
python test_app.py       # headless UI smoke test (needs streamlit)
```

To change the fit (e.g. weighted 1/Y²), add a checkbox in `app.py` that passes
`weighted=True` to `fit_4pl()`. To use a different plate layout, edit
`default_standards()` / `default_samples()`.

---

## 中文速查（给 TA）

- 这是一个 Streamlit 网页：学生自己输 OD 读数 → 自动算校正吸光度 →
  **拟合 4PL 标准曲线**（不是直线！）→ 反算废水样品 IFN-γ 浓度并乘回稀释倍数。
  App 里没有任何预置答案。
- **改动**：`calc.py` 从线性回归换成 4PL 非线性拟合；结果页显示 A/B/C/D、
  R²/RMSE、4PL 曲线；`requirements.txt` 新增 `scipy`。
- 本地跑：`pip install -r requirements.txt` 然后 `streamlit run app.py`。
- 推 GitHub 后 Streamlit Cloud 会自动重新部署；标题应变成
  “ELISA 4PL Result Generator”。
- 验证：`python test_calc.py` 和 `python verify_method.py` 应全部通过，
  示例数据应得到 ≈394（原液）和 ≈352 pg/mL（1:2 稀释 ×2）。
- 学生填写：标准品至少 4 个点（推荐 Blank + S1–S4）；样品填稀释倍数
  （原样 = 1）。超出标准曲线范围会提示「外推」。
