#!/usr/bin/env python3
"""
Script 01: Prepare and filter mouse DEG data for mitochondrial analysis
Author: Analysis pipeline
Date: 2025-10-06
Updated: 2025-11-03 (Mouse-only focus)
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Setup paths
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results" / "tables"

# Ensure output directory exists
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

def load_degs(filepath, model_name):
    """Load differential expression data"""
    print(f"\n{'='*60}")
    print(f"Loading {model_name} DEG data...")
    print(f"{'='*60}")

    # Try loading with default settings
    df = pd.read_csv(filepath, sep="\t")
    print(f"Total genes in file: {len(df):,}")
    print(f"Columns: {list(df.columns)}")

    # Check if numeric columns use comma as decimal separator
    # This is common in European Excel exports
    if model_name == "Mouse (OLCC)":
        numeric_cols = ['baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']
        for col in numeric_cols:
            if col in df.columns:
                # Replace comma with period and convert to numeric
                df[col] = df[col].astype(str).str.replace(',', '.', regex=False)
                df[col] = pd.to_numeric(df[col], errors='coerce')
        print("Converted European decimal format (comma) to standard format (period)")

    return df

def filter_degs(df, padj_threshold=0.01, log2fc_threshold=0.5):
    """Apply filtering criteria to DEG data"""
    print(f"\nApplying filters:")
    print(f"  - padj < {padj_threshold}")
    print(f"  - |log2FoldChange| >= {log2fc_threshold}")

    # Convert padj to numeric, coercing errors to NaN
    df = df.copy()
    df['padj'] = pd.to_numeric(df['padj'], errors='coerce')
    df['log2FoldChange'] = pd.to_numeric(df['log2FoldChange'], errors='coerce')

    # Drop rows with missing values
    initial_count = len(df)
    df = df.dropna(subset=['padj', 'log2FoldChange'])
    if len(df) < initial_count:
        print(f"Dropped {initial_count - len(df)} rows with missing values")

    # Filter by padj
    filtered = df[df['padj'] < padj_threshold].copy()
    print(f"After padj filter: {len(filtered):,} genes")

    # Filter by log2FC
    filtered = filtered[np.abs(filtered['log2FoldChange']) >= log2fc_threshold].copy()
    print(f"After log2FC filter: {len(filtered):,} genes")

    # Add regulation direction
    filtered['regulation'] = filtered['log2FoldChange'].apply(
        lambda x: 'up' if x > 0 else 'down'
    )

    # Count up/down
    up = (filtered['regulation'] == 'up').sum()
    down = (filtered['regulation'] == 'down').sum()
    print(f"\nUpregulated: {up:,} genes ({up/len(filtered)*100:.1f}%)")
    print(f"Downregulated: {down:,} genes ({down/len(filtered)*100:.1f}%)")

    return filtered

def summarize_degs(df):
    """Generate summary statistics"""
    summary = {
        'total_degs': len(df),
        'upregulated': (df['regulation'] == 'up').sum(),
        'downregulated': (df['regulation'] == 'down').sum(),
        'mean_log2fc': df['log2FoldChange'].mean(),
        'median_log2fc': df['log2FoldChange'].median(),
        'min_log2fc': df['log2FoldChange'].min(),
        'max_log2fc': df['log2FoldChange'].max(),
        'mean_padj': df['padj'].mean(),
        'median_padj': df['padj'].median()
    }
    return summary

def main():
    """Main analysis pipeline"""
    print("\n" + "="*60)
    print("MITOCHONDRIAL ANALYSIS: Data Preparation")
    print("="*60)

    # Load mouse DEGs
    mouse_file = DATA_DIR / "mouse_degs.txt"
    mouse_df = load_degs(mouse_file, "Mouse (OLCC)")

    # Filter mouse DEGs
    mouse_filtered = filter_degs(mouse_df, padj_threshold=0.01, log2fc_threshold=0.5)

    # Save filtered DEGs
    output_file = RESULTS_DIR / "mouse_degs_filtered.csv"
    mouse_filtered.to_csv(output_file, index=False)
    print(f"\nSaved filtered DEGs to: {output_file}")

    # Generate and save summary
    summary = summarize_degs(mouse_filtered)
    summary_df = pd.DataFrame([summary])
    summary_file = RESULTS_DIR / "mouse_degs_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    print(f"Saved summary to: {summary_file}")

    # Print final summary
    print("\n" + "="*60)
    print("DATA PREPARATION COMPLETE")
    print("="*60)
    print(f"Mouse DEGs (filtered): {len(mouse_filtered):,}")
    print(f"\nOutput files saved to: {RESULTS_DIR}")

if __name__ == "__main__":
    main()
