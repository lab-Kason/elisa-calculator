# Sample data — 4PL example plates

These CSVs can be uploaded through the two file uploaders in the app
(tab **① 标准品** and tab **② 样品**), or read with
`pandas.read_csv()` in tests.

## Dataset A — `example_A_clean.csv` + `example_A_samples_clean.csv`

A plate generated from the reference **"Real 4PL model"**:

```
A = 0.0821, B = 1.0741, C = 3877.5, D = 10.7614
```

Standard rows: Blank (0) + S1–S4 (2000 / 1000 / 500 / 250 pg/mL).

Expected output with `calc.py` (unweighted 4PL fit):

| Sample | Corrected OD | Dilution | Measured conc | Reported conc | Status |
|---|---|---|---|---|---|
| Wastewater original | 0.925 | 1 | ≈393.8 pg/mL | ≈393.8 pg/mL | in range |
| Wastewater 1:2 | 0.454 | 2 | ≈175.8 pg/mL | ≈351.5 pg/mL | in range |
| Wastewater over-range demo | 4.200 | 1 | ≈2 509.9 pg/mL | ≈2 509.9 pg/mL | flagged (above range) |

Fitted parameters on this rounded plate are close to the reference values
(A≈0.0820, B≈1.0709, C≈3961, D≈10.913; R²≈0.9999998). This is the same
answer as the lecture slide: **≈394 pg/mL** for well 6 and **≈352 pg/mL** for
well 7 after the ×2 dilution.

## Dataset B — `example_B_messy.csv` + `example_B_samples_messy.csv`

A noisier plate used to exercise the flags and permissive behaviour:

- noisy standard points (fit is still fine, parameters shift a little);
- an over-range sample (OD 4.500) → flagged **高于标准曲线范围（外推）**;
- a row with **OD600 blank** → corrected OD falls back to OD450 (0.500);
- two replicate wells (A and B) that should agree closely.

Approximate results: original ≈403.9 pg/mL, 1:2 diluted ≈375.2 pg/mL,
1:10 ≈265.4 pg/mL, replicates ≈217.2 / 219.4 pg/mL, missing-OD600 ≈208.1
pg/mL, over-range ≈2 978.0 pg/mL (flagged). Exact numbers depend only on the
CSV contents and the script; `verify_method.py` cross-checks them with an
independent bisection implementation.
