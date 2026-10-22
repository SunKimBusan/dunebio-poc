"""
Step 1-5. Build Pan-Cancer master table and Train/Test split.

NSCLC cohort (Test Set, zero-shot transfer target):
  TISSUE == 'Lung' AND CANCER_TYPE in {Non-Small Cell Lung Carcinoma, Squamous Cell Lung Carcinoma}
  (whitelist; the blacklist 'Lung & != SCLC' would also pull in 21 Mesothelioma
   + 8 'Other Solid Cancers' = 132 rows, which are not NSCLC histologies)
Train Set: all remaining Pan-Cancer cell lines.

Outputs:
  data/processed/osimertinib_pancancer_master.csv (957)
  data/processed/train_pancancer.csv              (854)
  data/processed/test_nsclc.csv                   (103)
"""

import pandas as pd
from pathlib import Path

INPUT_FILE = Path("data/processed/osimertinib_with_egfr_expr.csv")
OUT_DIR = Path("data/processed")
MASTER_FILE = OUT_DIR / "osimertinib_pancancer_master.csv"
TRAIN_FILE = OUT_DIR / "train_pancancer.csv"
TEST_FILE = OUT_DIR / "test_nsclc.csv"

NSCLC_TYPES = {"Non-Small Cell Lung Carcinoma", "Squamous Cell Lung Carcinoma"}
EXPECTED = {"master": 957, "train": 854, "test": 103}


def main():
    df = pd.read_csv(INPUT_FILE, encoding="utf-8-sig")
    df["IS_NSCLC"] = (df["TISSUE"] == "Lung") & df["CANCER_TYPE"].isin(NSCLC_TYPES)

    train = df[~df["IS_NSCLC"]].reset_index(drop=True)
    test = df[df["IS_NSCLC"]].reset_index(drop=True)

    assert train["SANGER_MODEL_ID"].isin(test["SANGER_MODEL_ID"]).sum() == 0, "Train/Test leakage"
    for name, n in {"master": len(df), "train": len(train), "test": len(test)}.items():
        assert n == EXPECTED[name], f"{name}: {n} != {EXPECTED[name]}"

    df.to_csv(MASTER_FILE, index=False)
    train.to_csv(TRAIN_FILE, index=False)
    test.to_csv(TEST_FILE, index=False)

    print(f"Master: {len(df)} | Train: {len(train)} | Test: {len(test)}")
    summary = pd.DataFrame({
        "Train": train["EGFR_STATUS"].value_counts(),
        "Test": test["EGFR_STATUS"].value_counts(),
    }).fillna(0).astype(int)
    print("\nEGFR_STATUS:\n" + summary.to_string())
    print("\nTest CANCER_TYPE:\n" + test["CANCER_TYPE"].value_counts().to_string())
    print(f"\nEGFR_EXPR_MISSING  Train: {train['EGFR_EXPR_MISSING'].sum()} | Test: {test['EGFR_EXPR_MISSING'].sum()}")
    sens = test[test["EGFR_STATUS"].isin(["Exon19del", "L858R"])]
    print(f"\nTest sensitizing mutants (n={len(sens)}):")
    print(sens[["CELL_LINE_NAME", "EGFR_STATUS", "EGFR_T790M", "LN_IC50", "EGFR_LOG2TPM", "EGFR_EXPR_MISSING"]].to_string(index=False))
    print(f"\nSaved -> {MASTER_FILE}, {TRAIN_FILE}, {TEST_FILE}")


if __name__ == "__main__":
    main()
