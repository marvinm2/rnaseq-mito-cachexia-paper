#!/usr/bin/env python3
"""
Script 15: Supplementary muscle DEG tables.

Combines existing GO term outputs into two flat supplementary tables:

1. Muscle cell differentiation DEGs (35 genes) - from go_myogenesis_degs.csv
   (children of GO:0042692: striated + smooth muscle cell differentiation).

2. Skeletal muscle development DEGs (32 genes) - union of three children of
   GO:0007519 already saved as per-term files by script 09:
     - GO:0035914 skeletal muscle cell differentiation
     - GO:0048741 skeletal muscle fiber development
     - GO:0048630 skeletal muscle tissue growth

Pure dedup/union/join over existing artifacts; no GO walking, no annotation
fetching.

Outputs:
- results/tables/Supplementary_Table_Muscle_Cell_Differentiation_DEGs.csv
- results/tables/Supplementary_Table_Skeletal_Muscle_Development_DEGs.csv
"""

from pathlib import Path
import pandas as pd

BASE_DIR = Path(__file__).parent.parent
TABLES_DIR = BASE_DIR / "results" / "tables"

DEG_STATS_COLS = [
    "Ensembl_ID", "GeneSymbol", "baseMean", "log2FoldChange", "lfcSE",
    "stat", "pvalue", "padj", "regulation",
]


def load_deg_stats():
    df = pd.read_csv(TABLES_DIR / "mouse_degs_filtered.csv")
    df = df.rename(columns={
        "Unnamed: 0": "Ensembl_ID",
        "Unnamed: 1": "Ensembl_ID_dup",
        "Unnamed: 2": "GeneSymbol",
    })
    df["GeneSymbol_upper"] = df["GeneSymbol"].str.upper()
    return df


def build_muscle_cell_differentiation(deg):
    src = pd.read_csv(TABLES_DIR / "go_myogenesis_degs.csv")
    src["GeneSymbol_upper"] = src["GeneSymbol"].str.upper()

    subcats = (
        src.groupby("GeneSymbol_upper")["Subcategory"]
        .apply(lambda s: "; ".join(sorted(set(s))))
        .reset_index()
        .rename(columns={"Subcategory": "GO_subcategories"})
    )

    merged = deg.merge(subcats, on="GeneSymbol_upper", how="inner")
    out = merged[DEG_STATS_COLS + ["GO_subcategories"]].copy()
    return out.sort_values("log2FoldChange", ascending=False).reset_index(drop=True)


def build_skeletal_muscle_development(deg):
    children = {
        "GO:0035914": ("skeletal muscle cell differentiation", "go_GO_0035914_genes.csv"),
        "GO:0048741": ("skeletal muscle fiber development",    "go_GO_0048741_genes.csv"),
        "GO:0048630": ("skeletal muscle tissue growth",        "go_GO_0048630_genes.csv"),
    }

    rows = []
    for go_id, (name, fname) in children.items():
        g = pd.read_csv(TABLES_DIR / fname)
        for sym in g["gene_symbol"]:
            rows.append({"GeneSymbol_upper": sym.upper(), "Subcategory": name})
    long_df = pd.DataFrame(rows)

    subcats = (
        long_df.groupby("GeneSymbol_upper")["Subcategory"]
        .apply(lambda s: "; ".join(sorted(set(s))))
        .reset_index()
        .rename(columns={"Subcategory": "GO_subcategories"})
    )

    merged = deg.merge(subcats, on="GeneSymbol_upper", how="inner")
    out = merged[DEG_STATS_COLS + ["GO_subcategories"]].copy()
    return out.sort_values("log2FoldChange", ascending=False).reset_index(drop=True)


def sanity_check(df, name, expected_n):
    assert df["GeneSymbol"].is_unique, f"{name}: duplicate GeneSymbol"
    assert (df["padj"] < 0.01).all(), f"{name}: padj cutoff violated"
    assert (df["log2FoldChange"].abs() >= 0.5).all(), f"{name}: |log2FC| cutoff violated"
    assert df[["baseMean", "log2FoldChange", "padj"]].notna().all().all(), f"{name}: NaN in core stats"
    assert df["GO_subcategories"].str.len().gt(0).all(), f"{name}: empty GO_subcategories"
    assert len(df) == expected_n, f"{name}: expected {expected_n} rows, got {len(df)}"


def main():
    deg = load_deg_stats()
    print(f"Loaded {len(deg)} filtered DEGs.")

    mcd = build_muscle_cell_differentiation(deg)
    smd = build_skeletal_muscle_development(deg)

    sanity_check(mcd, "muscle cell differentiation", 35)
    sanity_check(smd, "skeletal muscle development", 32)

    p_mcd = TABLES_DIR / "Supplementary_Table_Muscle_Cell_Differentiation_DEGs.csv"
    p_smd = TABLES_DIR / "Supplementary_Table_Skeletal_Muscle_Development_DEGs.csv"
    mcd.to_csv(p_mcd, index=False)
    smd.to_csv(p_smd, index=False)

    print()
    print(f"Muscle cell differentiation DEGs: {len(mcd)} genes "
          f"(up={(mcd['regulation']=='up').sum()}, down={(mcd['regulation']=='down').sum()})")
    print(f"  -> {p_mcd}")
    print()
    print(f"Skeletal muscle development DEGs: {len(smd)} genes "
          f"(up={(smd['regulation']=='up').sum()}, down={(smd['regulation']=='down').sum()})")
    print(f"  -> {p_smd}")


if __name__ == "__main__":
    main()
