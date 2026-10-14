"""
Step 1-3. Label EGFR mutation status per cell line using Cell Model Passports Mutations Summary.

EGFR_STATUS (4 classes; priority: Exon19del > L858R > Other > Wildtype):
- Exon19del : in-frame deletion / delins within exon 19 (codons 729-761, e.g. E746_A750del)
- L858R     : p.L858R
- Other     : any other EGFR coding mutation (atypical: G719S, L861Q, S768I, ...).
              Kept as an independent class and used in training (not merged into Wildtype).
- Wildtype  : Wildtype (Driver-negative) = no EGFR variant in Mutations Summary.
              The Mutations Summary lists only driver-annotated variants, so absence means
              "no EGFR driver mutation", not full-sequence wildtype.
Co-occurring T790M is recorded separately in EGFR_T790M.
One-hot columns EGFR_WT / EGFR_EX19DEL / EGFR_L858R / EGFR_OTHER are added for modeling.

NSCLC cohort (IS_NSCLC, used as Test Set in Step 1-5):
  TISSUE == 'Lung' AND CANCER_TYPE in {Non-Small Cell Lung Carcinoma, Squamous Cell Lung Carcinoma}
  i.e. lung tissue minus SCLC, additionally excluding pleural Mesothelioma and
  'Other Solid Cancers' (carcinoid / mucoepidermoid etc.), which are not NSCLC histologies.
"""

import re
import pandas as pd
from pathlib import Path

INPUT_FILE = Path("data/processed/osimertinib_with_metadata.csv")
MUT_FILE = Path("data/raw/mutations_summary_20260724.csv")
OUTPUT_FILE = Path("data/processed/osimertinib_with_egfr.csv")

EXON19_RANGE = range(729, 762)
NSCLC_TYPES = {"Non-Small Cell Lung Carcinoma", "Squamous Cell Lung Carcinoma"}
ONEHOT = {"Wildtype": "EGFR_WT", "Exon19del": "EGFR_EX19DEL", "L858R": "EGFR_L858R", "Other": "EGFR_OTHER"}


def is_exon19_del(p: str) -> bool:
    m = re.match(r"p\.[A-Z](\d+)_[A-Z](\d+)(del|delins)", p)
    return bool(m) and int(m.group(1)) in EXON19_RANGE and int(m.group(2)) in EXON19_RANGE


def classify(mutations: list[str]) -> str:
    if any(is_exon19_del(p) for p in mutations):
        return "Exon19del"
    if "p.L858R" in mutations:
        return "L858R"
    return "Other" if mutations else "Wildtype"


def main():
    df = pd.read_csv(INPUT_FILE, encoding="utf-8-sig")
    mut = pd.read_csv(MUT_FILE, low_memory=False)

    egfr = mut[(mut["gene_symbol"] == "EGFR") & (mut["coding"] == "t")]
    per_model = (
        egfr.groupby("model_id")["protein_mutation"]
        .apply(lambda s: sorted(set(s.dropna())))
        .rename("muts")
    )

    df = df.merge(per_model, left_on="SANGER_MODEL_ID", right_index=True, how="left")
    df["muts"] = df["muts"].apply(lambda x: x if isinstance(x, list) else [])
    df["EGFR_STATUS"] = df["muts"].apply(classify)
    df["EGFR_T790M"] = df["muts"].apply(lambda m: "p.T790M" in m)
    df["EGFR_MUTATIONS"] = df["muts"].apply(lambda m: ";".join(m))
    df = df.drop(columns="muts")
    for status, col in ONEHOT.items():
        df[col] = (df["EGFR_STATUS"] == status).astype(int)
    df["IS_NSCLC"] = (df["TISSUE"] == "Lung") & df["CANCER_TYPE"].isin(NSCLC_TYPES)

    df.to_csv(OUTPUT_FILE, index=False, encoding="utf-8-sig")

    nsclc = df[df["IS_NSCLC"]]
    print("[NSCLC cohort composition]")
    print(nsclc["CANCER_TYPE"].value_counts().to_string())
    print(f"Train (non-NSCLC): {len(df) - len(nsclc)} / Test (NSCLC): {len(nsclc)}\n")
    print(f"Saved -> {OUTPUT_FILE} ({len(df)} rows)\n")
    print("[Pan-Cancer EGFR_STATUS]")
    print(df["EGFR_STATUS"].value_counts().to_string())
    print(f"\n[NSCLC (n={len(nsclc)}) EGFR_STATUS]")
    vc = nsclc["EGFR_STATUS"].value_counts()
    print(pd.DataFrame({"n": vc, "%": (vc / len(nsclc) * 100).round(1)}).to_string())
    print("\n[EGFR-mutant cell lines]")
    cols = ["CELL_LINE_NAME", "CANCER_TYPE", "IS_NSCLC", "EGFR_STATUS", "EGFR_T790M", "EGFR_MUTATIONS", "LN_IC50"]
    print(df[df["EGFR_STATUS"] != "Wildtype"][cols]
          .sort_values(["EGFR_STATUS", "LN_IC50"]).to_string(index=False))


if __name__ == "__main__":
    main()
