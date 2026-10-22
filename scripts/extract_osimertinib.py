"""
Extract Osimertinib dose-response data from GDSC2 Fitted Dose Response dataset.

Notes on GDSC2 release 27Oct23:
- DRUG_NAME: 'Osimertinib' (AZD9291)
- DRUG_ID: 1919 (Note: In GDSC2, DRUG_ID 1084 corresponds to 'Rapamycin', while Osimertinib is assigned DRUG_ID 1919)
- Target: EGFR (EGFR signaling pathway)
"""

import sys
import openpyxl
import pandas as pd
from pathlib import Path

def extract_osimertinib(
    input_file: Path = Path("data/raw/GDSC2_fitted_dose_response_27Oct23.xlsx"),
    output_file: Path = Path("data/processed/gdsc2_osimertinib_fitted_dose_response.csv")
):
    print(f"Reading {input_file}...")
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Use openpyxl read_only mode for memory efficiency with large ~22MB excel file
    wb = openpyxl.load_workbook(input_file, read_only=True, data_only=True)
    sheet = wb.active
    rows = sheet.iter_rows(values_only=True)
    header = list(next(rows))

    drug_name_idx = header.index("DRUG_NAME")
    drug_id_idx = header.index("DRUG_ID")

    osimertinib_rows = []
    rapamycin_1084_rows = []

    for r in rows:
        d_name = str(r[drug_name_idx]) if r[drug_name_idx] is not None else ""
        d_id = r[drug_id_idx]

        # Check Osimertinib by name or DRUG_ID 1919
        if d_name.lower() in ["osimertinib", "azd9291", "azd-9291"] or d_id == 1919:
            osimertinib_rows.append(r)
        elif d_id == 1084:
            rapamycin_1084_rows.append(r)

    df_osim = pd.DataFrame(osimertinib_rows, columns=header)
    
    print(f"\n[Filtering Results]")
    print(f"- Total Osimertinib rows extracted (DRUG_NAME='Osimertinib' / DRUG_ID=1919): {len(df_osim)}")
    print(f"- Rows with DRUG_ID=1084 (Rapamycin, MTORC1 inhibitor): {len(rapamycin_1084_rows)}")
    
    # Save Osimertinib data
    df_osim.to_csv(output_file, index=False, encoding="utf-8-sig")
    print(f"\nSuccessfully saved filtered Osimertinib data to: {output_file}")
    
    # Print summary of the extracted data
    print("\n[Summary of Osimertinib Data]")
    print(f"- Number of Cell Lines: {df_osim['CELL_LINE_NAME'].nunique()}")
    print(f"- Cancer Types represented: {df_osim['CANCER_TYPE'].nunique()}")
    print("\nTop 5 Cancer Types:")
    print(df_osim['CANCER_TYPE'].value_counts().head(5).to_string())
    print("\nMetrics Summary (LN_IC50, AUC, RMSE, Z_SCORE):")
    print(df_osim[['LN_IC50', 'AUC', 'RMSE', 'Z_SCORE']].describe().to_string())
    
    return df_osim

if __name__ == "__main__":
    extract_osimertinib()
