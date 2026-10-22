"""
Step 1-2. Join Cell Model Passports Model List onto Osimertinib GDSC2 data.

Key: GDSC2.SANGER_MODEL_ID == model_list.model_id
Note: GDSC2 already has a 'CANCER_TYPE' column (TCGA-style codes, e.g. LUAD, BRCA).
      Model List 'cancer_type' (e.g. 'Non-Small Cell Lung Carcinoma') is renamed to
      MODEL_CANCER_TYPE to avoid collision; GDSC2 code is kept as GDSC_TCGA_TYPE.
"""

import pandas as pd
from pathlib import Path

DRUG_FILE = Path("data/processed/gdsc2_osimertinib_fitted_dose_response.csv")
MODEL_FILE = Path("data/raw/model_list_20260814.csv")
OUTPUT_FILE = Path("data/processed/osimertinib_with_metadata.csv")

MODEL_COLS = {
    "model_id": "SANGER_MODEL_ID",
    "model_name": "MODEL_NAME",
    "tissue": "TISSUE",
    "cancer_type": "MODEL_CANCER_TYPE",
    "cancer_type_detail": "CANCER_TYPE_DETAIL",
    "sample_site": "SAMPLE_SITE",
    "msi_status": "MSI_STATUS",
    "COSMIC_ID": "COSMIC_ID",
    "CCLE_ID": "CCLE_ID",
}


def main():
    drug = pd.read_csv(DRUG_FILE, encoding="utf-8-sig")
    models = pd.read_csv(MODEL_FILE, low_memory=False)

    models = models[[c for c in MODEL_COLS if c in models.columns]].rename(columns=MODEL_COLS)
    models = models.drop_duplicates("SANGER_MODEL_ID")

    drug = drug.rename(columns={"CANCER_TYPE": "GDSC_TCGA_TYPE"})
    merged = drug.merge(models, on="SANGER_MODEL_ID", how="left", validate="many_to_one")

    unmatched = merged["MODEL_CANCER_TYPE"].isna()
    # Final CANCER_TYPE: Model List value, fallback to GDSC TCGA code
    merged["CANCER_TYPE"] = merged["MODEL_CANCER_TYPE"].fillna(merged["GDSC_TCGA_TYPE"])

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(OUTPUT_FILE, index=False, encoding="utf-8-sig")

    print(f"Rows: {len(merged)} (input {len(drug)})")
    print(f"Unique cell lines: {merged['SANGER_MODEL_ID'].nunique()}")
    print(f"Unmatched in Model List: {unmatched.sum()}")
    if unmatched.any():
        print(merged.loc[unmatched, ["CELL_LINE_NAME", "SANGER_MODEL_ID"]].to_string(index=False))
    print("\nTop 15 CANCER_TYPE:")
    print(merged["CANCER_TYPE"].value_counts().head(15).to_string())
    nsclc = merged["CANCER_TYPE"].str.contains("Non-Small Cell Lung", na=False)
    print(f"\nNSCLC cell lines: {nsclc.sum()}")
    print(f"Saved -> {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
