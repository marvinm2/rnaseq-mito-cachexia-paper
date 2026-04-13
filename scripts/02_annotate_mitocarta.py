#!/usr/bin/env python3
"""
Script 02: Annotate DEGs with MitoCarta 3.0
Author: Analysis pipeline
Date: 2025-10-06
"""

import pandas as pd
import numpy as np
from pathlib import Path
import re

# Setup paths
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results" / "tables"

def load_mitocarta(filepath):
    """Load MitoCarta 3.0 database"""
    print(f"\n{'='*60}")
    print("Loading MitoCarta 3.0 database...")
    print(f"{'='*60}")

    mc = pd.read_csv(filepath)
    print(f"Total mitochondrial genes in MitoCarta: {len(mc):,}")
    print(f"Columns: {list(mc.columns)[:10]}... ({len(mc.columns)} total)")

    # Show distribution of pathways
    if 'MitoCarta3.0_MitoPathways' in mc.columns:
        pathways = mc['MitoCarta3.0_MitoPathways'].dropna()
        print(f"\nGenes with pathway annotation: {len(pathways):,}")

    return mc

def extract_ensembl_id(gene_id_str):
    """Extract clean Ensembl ID from string"""
    if pd.isna(gene_id_str):
        return None
    # Handle various formats: ENSMUSG00000000001 or similar
    match = re.search(r'ENSMUSG\d+', str(gene_id_str))
    if match:
        return match.group(0)
    return str(gene_id_str).strip()

def annotate_degs_with_mitocarta(degs_df, mitocarta_df, model_name="Mouse"):
    """Annotate DEGs with MitoCarta information"""
    print(f"\n{'='*60}")
    print(f"Annotating {model_name} DEGs with MitoCarta...")
    print(f"{'='*60}")

    # Prepare DEG identifiers
    # The DEG file has Ensembl IDs in different columns
    if 'Unnamed: 0' in degs_df.columns:
        degs_df['EnsemblID'] = degs_df['Unnamed: 0'].apply(extract_ensembl_id)
    else:
        print("Warning: Could not find Ensembl ID column in DEG data")
        degs_df['EnsemblID'] = degs_df.index

    # Prepare MitoCarta identifiers
    if 'MouseGeneID' in mitocarta_df.columns:
        # Clean up the mouse gene IDs
        mitocarta_df['MouseGeneID_clean'] = mitocarta_df['MouseGeneID'].apply(
            lambda x: str(int(x)) if pd.notna(x) and x != '' else None
        )
        print(f"MitoCarta mouse genes: {mitocarta_df['MouseGeneID_clean'].notna().sum():,}")

        # Note: Now using Mouse MitoCarta 3.0 directly (not human with orthologs)

    print(f"\nTotal DEGs to annotate: {len(degs_df):,}")

    # Try matching by Ensembl ID pattern
    # First, let's check what identifiers are available
    print("\nSample DEG identifiers:")
    print(degs_df[['EnsemblID']].head())

    print("\nSample MitoCarta identifiers:")
    if 'Symbol' in mitocarta_df.columns:
        print(mitocarta_df[['Symbol', 'MouseGeneID']].head())

    # Create a mapping from Ensembl ID to MitoCarta info
    # We'll use the gene symbol approach since we have Ensembl IDs

    # Extract gene symbols from DEG file
    if 'Unnamed: 2' in degs_df.columns:
        degs_df['GeneSymbol'] = degs_df['Unnamed: 2']
    elif 'Unnamed: 1' in degs_df.columns:
        # Try the second unnamed column
        degs_df['GeneSymbol'] = degs_df['Unnamed: 1']

    print(f"\nSample gene symbols from DEGs:")
    if 'GeneSymbol' in degs_df.columns:
        print(degs_df['GeneSymbol'].head())

    # Match DEGs to MitoCarta by gene symbol (most reliable for mouse)
    if 'Symbol' in mitocarta_df.columns and 'GeneSymbol' in degs_df.columns:
        # Merge on gene symbol
        mitocarta_symbols = mitocarta_df['Symbol'].str.upper().unique()
        print(f"\nUnique gene symbols in MitoCarta: {len(mitocarta_symbols):,}")

        # Case-insensitive matching
        degs_df['GeneSymbol_upper'] = degs_df['GeneSymbol'].str.upper()

        # Merge with MitoCarta
        annotated = degs_df.merge(
            mitocarta_df,
            left_on='GeneSymbol_upper',
            right_on=mitocarta_df['Symbol'].str.upper(),
            how='left',
            suffixes=('', '_mitocarta')
        )

        # Identify mitochondrial genes
        annotated['is_mitochondrial'] = annotated['MitoCarta3.0_List'].notna()
        n_mito = annotated['is_mitochondrial'].sum()

        print(f"\n{'='*60}")
        print(f"ANNOTATION RESULTS")
        print(f"{'='*60}")
        print(f"Total DEGs: {len(annotated):,}")
        print(f"Mitochondrial DEGs: {n_mito:,} ({n_mito/len(annotated)*100:.1f}%)")

        # Get mitochondrial DEGs only
        mito_degs = annotated[annotated['is_mitochondrial']].copy()

        # Count up/down regulation
        if 'regulation' in mito_degs.columns:
            n_up = (mito_degs['regulation'] == 'up').sum()
            n_down = (mito_degs['regulation'] == 'down').sum()
            print(f"  - Upregulated: {n_up:,} ({n_up/len(mito_degs)*100:.1f}%)")
            print(f"  - Downregulated: {n_down:,} ({n_down/len(mito_degs)*100:.1f}%)")

        return annotated, mito_degs
    else:
        print("ERROR: Could not find matching columns for annotation")
        return degs_df, pd.DataFrame()

def parse_mitocarta_pathways(pathway_str):
    """Parse MitoCarta pathway string into categories"""
    if pd.isna(pathway_str):
        return []

    # Split by pipe
    pathways = str(pathway_str).split('|')
    pathways = [p.strip() for p in pathways]

    categories = set()
    for pathway in pathways:
        # Extract main category
        if 'Complex I' in pathway:
            categories.add('Complex I')
        if 'Complex II' in pathway:
            categories.add('Complex II')
        if 'Complex III' in pathway:
            categories.add('Complex III')
        if 'Complex IV' in pathway:
            categories.add('Complex IV')
        if 'Complex V' in pathway or 'ATP synth' in pathway:
            categories.add('Complex V')
        if 'TCA cycle' in pathway or 'Citric' in pathway:
            categories.add('TCA cycle')
        if 'ribosom' in pathway.lower():
            categories.add('Mitochondrial ribosome')
        if 'Translation' in pathway:
            categories.add('Translation factors')
        if 'SLC25' in pathway or 'carrier' in pathway.lower():
            categories.add('Transporters')

    return list(categories)

def main():
    """Main analysis pipeline"""
    print("\n" + "="*60)
    print("MITOCHONDRIAL ANALYSIS: MitoCarta Annotation")
    print("="*60)

    # Load MitoCarta
    mitocarta_file = DATA_DIR / "mitocarta.csv"
    mitocarta_df = load_mitocarta(mitocarta_file)

    # Load filtered DEGs
    mouse_degs_file = RESULTS_DIR / "mouse_degs_filtered.csv"
    mouse_degs = pd.read_csv(mouse_degs_file)
    print(f"\nLoaded {len(mouse_degs):,} filtered mouse DEGs")

    # Annotate mouse DEGs
    mouse_annotated, mouse_mito = annotate_degs_with_mitocarta(
        mouse_degs, mitocarta_df, "Mouse"
    )

    # Parse pathway categories
    if len(mouse_mito) > 0:
        mouse_mito['pathway_categories'] = mouse_mito['MitoCarta3.0_MitoPathways'].apply(
            parse_mitocarta_pathways
        )

        # Show pathway distribution
        print(f"\n{'='*60}")
        print("Pathway Distribution")
        print(f"{'='*60}")

        all_categories = []
        for cats in mouse_mito['pathway_categories']:
            all_categories.extend(cats)

        if all_categories:
            pathway_counts = pd.Series(all_categories).value_counts()
            print(pathway_counts.to_string())

    # Save results
    # Full annotated file
    output_file = RESULTS_DIR / "mouse_degs_annotated_full.csv"
    mouse_annotated.to_csv(output_file, index=False)
    print(f"\nSaved full annotated DEGs to: {output_file}")

    # Mitochondrial DEGs only
    mito_output = RESULTS_DIR / "mitochondrial_degs.csv"
    mouse_mito.to_csv(mito_output, index=False)
    print(f"Saved mitochondrial DEGs to: {mito_output}")

    # Create a clean summary table for the manuscript
    if len(mouse_mito) > 0:
        summary_cols = [
            'GeneSymbol', 'log2FoldChange', 'padj', 'regulation',
            'MitoCarta3.0_MitoPathways', 'MitoCarta3.0_SubMitoLocalization',
            'pathway_categories'
        ]
        available_cols = [col for col in summary_cols if col in mouse_mito.columns]

        summary = mouse_mito[available_cols].copy()
        summary = summary.sort_values('padj')

        summary_output = RESULTS_DIR / "mitochondrial_degs_summary.xlsx"
        summary.to_excel(summary_output, index=False)
        print(f"Saved summary table to: {summary_output}")

    print("\n" + "="*60)
    print("MITOCARTA ANNOTATION COMPLETE")
    print("="*60)

if __name__ == "__main__":
    main()
