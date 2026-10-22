# Data Setup

1. Download the raw data from [CancerRxGene](https://www.cancerrxgene.org/).
2. Place the downloaded files into the `data/raw/` directory.

### Required Files
- `GDSC2_fitted_dose_response_27Oct23.xlsx`
- `model_list_20260814.csv`
- `mutations_summary_20260724.csv`
- `rnaseq_merged_rsem_tpm_20260323.csv` (Extracted from `rnaseq_merged_20260323.zip`; only this TPM matrix file is required for gene expression features)

### Data Definitions
- **EGFR_STATUS** (4 classes): `Exon19del`, `L858R`, `Other`, `Wildtype`
  - `Other`: atypical EGFR variants (G719S, L861Q, S768I, ...). Kept as its own class and used for training.
  - `Wildtype`: defined as **EGFR driver-negative**, i.e. no EGFR variant in the Mutations Summary
    (which lists driver-annotated variants only).
  - One-hot columns: `EGFR_WT`, `EGFR_EX19DEL`, `EGFR_L858R`, `EGFR_OTHER`; T790M flagged in `EGFR_T790M`.
- **NSCLC cohort** (`IS_NSCLC`, Test Set): `TISSUE == 'Lung'` and `CANCER_TYPE` in
  {`Non-Small Cell Lung Carcinoma`, `Squamous Cell Lung Carcinoma`}.
  SCLC, pleural Mesothelioma and lung `Other Solid Cancers` (carcinoid, mucoepidermoid) are excluded.
  (A blacklist filter `Lung & != SCLC` would return 132 rows incl. 21 Mesothelioma + 8 Other.)
- **EGFR_LOG2TPM**: EGFR expression from the RNA-Seq matrix. Values in the file are **already
  log2(TPM + 1)** despite the file name, so no transform is applied.
  - Duplicate profiles per model (Sanger + Broad): Sanger preferred.
  - Coverage: 932 / 957. The 25 missing values are imputed with the overall median (4.158).
- **EGFR_EXPR_MISSING**: 1 = no RNA-Seq profile (imputed), 0 = measured.
- **EGFR_RNA_SOURCE**: `Sanger` / `Broad` / empty (no profile).

### Pipeline
Run from the repository root, in order:

| Step | Script | Output (`data/processed/`) |
|---|---|---|
| 1-1 | `scripts/extract_osimertinib.py` | `gdsc2_osimertinib_fitted_dose_response.csv` (957) |
| 1-2 | `scripts/join_model_metadata.py` | `osimertinib_with_metadata.csv` |
| 1-3 | `scripts/label_egfr_status.py` | `osimertinib_with_egfr.csv` |
| 1-4 | `scripts/add_egfr_expression.py` | `osimertinib_with_egfr_expr.csv` |
| 1-5 | `scripts/split_train_test.py` | `osimertinib_pancancer_master.csv` (957), `train_pancancer.csv` (854), `test_nsclc.csv` (103) |

### Split Summary
| EGFR_STATUS | Train (Pan-Cancer) | Test (NSCLC) |
|---|---|---|
| Exon19del | 0 | 5 |
| L858R | 0 | 2 |
| Other | 14 | 0 |
| Wildtype | 840 | 96 |
| **Total** | **854** | **103** |
| EGFR_EXPR_MISSING = 1 | 20 | 5 |

Test sensitizing mutants (n=7): NCI-H1650, PC-14, HCC-827, PC-3_[JPC-3], LOU-NH91 (Exon19del);
NCI-H1975 (L858R + T790M), H3255 (L858R).

### Known Caveats
- **No Exon19del/L858R in Train**: `EGFR_EX19DEL` / `EGFR_L858R` are constant (0) in Train, so the
  model cannot learn them directly; sensitivity signal must come from `EGFR_LOG2TPM` and `Other`.
- **Imputation leakage**: the median is computed on all 957 rows (incl. Test). For strict evaluation,
  recompute it on Train only.
- **Non-responding mutants**: NCI-H1650 and LOU-NH91 (Exon19del) show high LN_IC50 (~0.3–0.4).
