# Data Setup

1. Download the raw data from [CancerRxGene](https://www.cancerrxgene.org/).
2. Place the downloaded files into the `data/raw/` directory.

### Required Files
- `GDSC2_fitted_dose_response_27Oct23.xlsx`
- `model_list_20260814.csv`
- `mutations_summary_20260724.csv`

### Data Definitions
- **EGFR_STATUS** (4 classes): `Exon19del`, `L858R`, `Other`, `Wildtype`
  - `Other`: atypical EGFR variants (G719S, L861Q, S768I, ...). Kept as its own class and used for training.
  - `Wildtype`: defined as **EGFR driver-negative**, i.e. no EGFR variant in the Mutations Summary
    (which lists driver-annotated variants only).
  - One-hot columns: `EGFR_WT`, `EGFR_EX19DEL`, `EGFR_L858R`, `EGFR_OTHER`; T790M flagged in `EGFR_T790M`.
- **NSCLC cohort** (`IS_NSCLC`, Test Set): `TISSUE == 'Lung'` and `CANCER_TYPE` in
  {`Non-Small Cell Lung Carcinoma`, `Squamous Cell Lung Carcinoma`}.
  SCLC, pleural Mesothelioma and lung `Other Solid Cancers` (carcinoid, mucoepidermoid) are excluded.
