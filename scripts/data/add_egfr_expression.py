"""
Step 1-4. Add EGFR expression (RNA-Seq TPM) feature.

Source: Cell Model Passports rnaseq_merged_rsem_tpm (wide format)
  row 1: model_id, row 2: model_name, row 3: data_source (Sanger/Broad),
  row 4: column header (gene_symbol, ensembl_gene_id, gene_id), rows 5+: genes.
Only the EGFR row is streamed (file is ~350MB), so the full matrix is never loaded.

Duplicate profiles for one model (Sanger + Broad) -> Sanger preferred, else mean of rest.
NOTE: despite the file name, values are already log2(TPM + 1)
  (EGFR max ~11, HCC-827 amplified ~9.7), so no further transform is applied.
Output columns:
  EGFR_LOG2TPM      : log2(TPM + 1) as provided (modeling feature);
                      missing values imputed with the overall median (all 957 rows)
  EGFR_EXPR_MISSING : 1 = no RNA-Seq profile (imputed), 0 = measured
  EGFR_RNA_SOURCE : Sanger / Broad / NaN (no RNA-Seq profile)
"""

import csv
import pandas as pd
from pathlib import Path

INPUT_FILE = Path("data/processed/osimertinib_with_egfr.csv")
TPM_FILE = Path("data/raw/rnaseq_merged_rsem_tpm_20260323.csv")
OUTPUT_FILE = Path("data/processed/osimertinib_with_egfr_expr.csv")
GENE = "EGFR"


def load_gene_row(path: Path, gene: str) -> pd.DataFrame:
    with path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        model_ids = next(reader)[3:]
        next(reader)  # model_name
        sources = next(reader)[3:]
        next(reader)  # header
        for row in reader:
            if row[0] == gene:
                values = row[3:]
                break
        else:
            raise ValueError(f"{gene} not found in {path}")
    df = pd.DataFrame({"model_id": model_ids, "source": sources,
                       "tpm": pd.to_numeric(values, errors="coerce")})
    return df.dropna(subset=["tpm"])


def collapse(df: pd.DataFrame) -> pd.DataFrame:
    def pick(g):
        s = g[g["source"] == "Sanger"]
        g = s if len(s) else g
        return pd.Series({"EGFR_LOG2TPM": g["tpm"].mean(), "EGFR_RNA_SOURCE": "/".join(sorted(set(g["source"])))})
    return df.groupby("model_id").apply(pick, include_groups=False).reset_index()


def main():
    data = pd.read_csv(INPUT_FILE, encoding="utf-8-sig")
    raw = load_gene_row(TPM_FILE, GENE)
    print(f"EGFR profiles: {len(raw)} (unique models: {raw['model_id'].nunique()})")
    expr = collapse(raw)

    out = data.merge(expr, left_on="SANGER_MODEL_ID", right_on="model_id", how="left").drop(columns="model_id")
    assert len(out) == len(data)

    missing = out["EGFR_LOG2TPM"].isna()
    median = out["EGFR_LOG2TPM"].median()
    out["EGFR_EXPR_MISSING"] = missing.astype(int)
    out["EGFR_LOG2TPM"] = out["EGFR_LOG2TPM"].fillna(median)

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT_FILE, index=False)

    print(f"Rows: {len(out)} | with EGFR expr: {(~missing).sum()} | imputed: {missing.sum()} (median={median:.4f})")
    print(out["EGFR_RNA_SOURCE"].value_counts(dropna=False).to_string())
    print("\nEGFR_LOG2TPM by EGFR_STATUS (measured only):")
    print(out[~missing].groupby("EGFR_STATUS")["EGFR_LOG2TPM"].describe()[["count", "mean", "50%"]].round(2).to_string())
    print("\nImputed cell lines:", ", ".join(out.loc[missing, "CELL_LINE_NAME"]))
    print(f"Saved -> {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
