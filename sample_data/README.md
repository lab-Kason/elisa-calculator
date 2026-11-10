# Sample data (synthetic — for testing the app only, not real measurements)

Load a standards file in **① 标准品** and the matching samples file in
**② 样品** (upload, or copy-paste into the table). Use the *Standards template*
/ *Samples template* buttons in the sidebar if you need the exact header names.

## A — clean plate

`example_A_clean.csv` + `example_A_samples_clean.csv`

The standards sit exactly on `corrected OD = 0.0008 × concentration + 0.05`,
so every expected value can be derived by hand:

| Standard | Concentration (pg/mL) | OD450 | OD600 | corrected OD |
|---|---|---|---|---|
| S1 | 2000 | 1.69 | 0.04 | 1.65 |
| S2 | 1000 | 0.89 | 0.04 | 0.85 |
| S3 | 500 | 0.49 | 0.04 | 0.45 |
| S4 | 250 | 0.29 | 0.04 | 0.25 |

The wastewater is really **600 pg/mL**, measured twice (undiluted and 1:2):

| Sample | Dilution | corrected OD | measured conc | reported conc |
|---|---|---|---|---|
| Wastewater original | 1 | 0.530 | 600.0 | 600.0 |
| Wastewater 1:2 | 2 | 0.290 | 300.0 | 600.0 |

Expected fit: **slope 0.0008, intercept 0.05, R² = 1.00000**. This is the
textbook case — if the app shows anything else, something is broken.

## B — messy plate

`example_B_messy.csv` + `example_B_samples_messy.csv`

A realistic plate with scatter, out-of-range readings, replicate wells and a
missing OD600. The standards give **slope 0.00082737, intercept 0.00408696,
R² = 0.99997**, standard corrected-OD range 0.215 – 1.660.

| Sample | Dilution | corrected OD | expected status | expected reported conc (pg/mL) |
|---|---|---|---|---|
| Wastewater original | 1 | 1.824 | 高于标准曲线范围（外推） | 2199.6 |
| Wastewater 1:5 | 5 | 0.368 | 在标准曲线范围内 | 2199.2 |
| Wastewater 1:50 | 50 | 0.040 | 低于标准曲线范围（外推） | 2170.3 |
| Replicate A | 1 | 0.608 | 在标准曲线范围内 | 729.9 |
| Replicate B | 1 | 0.620 | 在标准曲线范围内 | 744.4 |
| OD600 missing demo | 1 | 0.501 (OD600 blank → OD450) | 在标准曲线范围内 | 600.6 |

Teaching points this plate demonstrates:

- the same wastewater measured at 1x and 1:5 agrees (**~2200 pg/mL**), so the
  dilution correction is right;
- the 1:50 well falls below the lowest standard, is flagged as extrapolation,
  and is the least accurate of the three — exactly why samples should be
  diluted into the working range;
- a blank OD600 is tolerated (corrected OD = OD450);
- replicate wells are entered as separate rows and can be averaged by hand.

All numbers above are reproduced by `python verify_method.py` in the repo root.
