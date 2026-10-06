# Breast Cancer Survival Analysis — METABRIC

MSc dissertation project: Evaluation of Risk Factors and Survival Rates of Patients with Breast Cancer Using Machine Learning and Traditional Methods.

## Overview

Predicting which breast cancer patients are at higher risk of death, and understanding what drives that risk, can help clinicians prioritise treatment and follow-up. This project compares a standard statistical model (Cox Proportional Hazards) against three machine learning survival models (Random Survival Forest, XGB-Cox, DeepSurv) on the METABRIC dataset.

## Data

[METABRIC](https://www.cbioportal.org/study/summary?id=brca_metabric) — 1,904 breast cancer patients with clinical and genomic data.

Features used: age at diagnosis, tumour size, number of positive lymph nodes, and ER/PR/HER2 receptor status.

## About this repo

The original project was built in R (`survival`, `survminer`, `randomForestSRC`, `xgboost`) and Python (`pycox`). This repo has two parts:

- `survival_analysis.py` — a Python version of the Kaplan-Meier, log-rank, and Cox PH analysis, built to run standalone without R. I cross-checked its output against my original results to make sure the rewrite was correct (see results below).
- `original_scripts/` — the original R scripts and the DeepSurv notebook, unchanged, including the ones for Random Survival Forest and XGB-Cox, which aren't re-implemented in Python here.

## Results

### Overall survival

![KM Overall](images/km_overall.png)

5-year survival: 96.5%. 10-year survival: 81.2%. Median survival: 197 months.

### Cox regression (multivariable)

![Cox Hazard Ratios](images/cox_hazard_ratios.png)

| Variable | Hazard ratio | p-value |
|---|---|---|
| Age at diagnosis | 0.991 | 0.004 |
| Tumour size | 0.999 | 0.673 |
| Lymph nodes positive | 1.023 | 0.052 |
| ER status (positive) | 0.970 | 0.777 |
| PR status (positive) | 0.902 | 0.244 |
| HER2 status (positive) | 1.208 | 0.115 |

Age was the only variable that stayed statistically significant once everything else was in the model. Lymph node involvement was borderline. ER, PR and HER2 status all looked important on their own (see the log-rank tests below) but lost significance once age, tumour size and lymph nodes were accounted for.

### Survival by receptor status

![KM by HER2 status](images/km_by_her2_status.png)

ER, PR, and HER2 status were each significantly associated with survival (log-rank p < 0.05 for all three). HER2-positive patients had noticeably worse survival than HER2-negative patients.

### Model comparison

| Model | C-index | AUC @12m | AUC @36m | AUC @60m |
|---|---|---|---|---|
| Cox PH | 0.448 | 0.568 | 0.600 | 0.565 |
| Random Survival Forest | 0.534 | 0.488 | 0.529 | 0.520 |
| XGB-Cox | 0.244 | 0.739 | 0.757 | 0.758 |
| DeepSurv | 0.240 (unstable) | — | — | — |

C-index measures how well a model ranks patients by risk over the whole follow-up period. Time-dependent AUC measures how well it separates patients at a specific point in time (e.g. who dies within 3 years). RSF had the best overall ranking; XGB-Cox was better at specific time windows. The two metrics don't always agree, which is normal.

Across all four models, age, lymph node involvement, and tumour size were consistently the strongest predictors. Receptor status added less once those were accounted for.

## Takeaway

No model was best across the board. Random Survival Forest ranked patients best overall, XGB-Cox was better at short/medium-term predictions, and Cox PH — the simplest and most interpretable model — stayed competitive with both. DeepSurv performed worst, which isn't unusual on a dataset this size (~1,900 patients), since deep learning models typically need a lot more data.

## Limitations

- Single cohort (METABRIC) — results may not generalise to other populations.
- Receptor-status subgroups are uneven in size (HER2-positive is only 12% of the cohort), so estimates for that group are less precise.
- This is a research comparison, not a validated clinical tool.

## Repo structure

```
breast-cancer-survival-analysis/
├── survival_analysis.py       # KM, log-rank, Cox PH (Python)
├── original_scripts/          # original R scripts + DeepSurv notebook
├── data/                      # dataset + exported result tables
├── images/                    # charts
└── README.md
```

## Running it

```bash
pip install pandas numpy scipy matplotlib
python survival_analysis.py
```

Random Survival Forest, XGB-Cox, and DeepSurv need R or the Python packages listed in `original_scripts/`.

## Next steps

- Run `original_scripts/` in an R environment to regenerate the ML model results directly here
- Add a published risk score (e.g. Nottingham Prognostic Index) as another baseline
