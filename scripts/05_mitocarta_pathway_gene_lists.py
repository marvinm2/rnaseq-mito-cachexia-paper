#!/usr/bin/env python3
"""
MitoCarta Pathway Gene Extraction
==================================

This script extracts gene lists for each MitoCarta pathway across all 3 hierarchical levels
to enable Fisher's exact statistical testing.

The MitoCarta3.0_MitoPathways column uses pipe "|" to separate multiple pathways
and ">" to separate hierarchical levels within a pathway.

Example pathway annotation:
"OXPHOS > Complex I > CI subunits | OXPHOS > OXPHOS subunits"

This gene belongs to:
- Level 0: OXPHOS
- Level 1: Complex I, OXPHOS subunits
- Level 2: CI subunits

Output:
-------
Creates gene list CSV files for each of the 104 pathways:
- 7 Level 0 pathways
- 20 Level 1 pathways
- 77 Level 2 pathways

Saved to: results/tables/pathway_genes/

Author: Bioinformatics Team
Date: December 2025
"""

import pandas as pd
from pathlib import Path
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')


def pathway_matches(pathway_string, target_pathway, level):
    """
    Check if target pathway appears at specified level in pathway string.

    Parameters
    ----------
    pathway_string : str
        MitoCarta pathway annotation (e.g., "OXPHOS > Complex I > CI subunits")
    target_pathway : str
        Pathway name to search for
    level : int
        Hierarchy level (0, 1, or 2)

    Returns
    -------
    bool
        True if target pathway found at specified level

    Examples
    --------
    >>> pathway_matches("OXPHOS > Complex I > CI subunits", "OXPHOS", 0)
    True
    >>> pathway_matches("OXPHOS > Complex I > CI subunits", "Complex I", 1)
    True
    >>> pathway_matches("OXPHOS > Complex I > CI subunits", "CI subunits", 2)
    True
    >>> pathway_matches("OXPHOS > Complex I > CI subunits", "OXPHOS", 1)
    False
    """
    if pd.isna(pathway_string):
        return False

    # Split by '|' to handle multiple pathways
    pathways = pathway_string.split('|')

    for pathway in pathways:
        # Split by '>' to get hierarchy levels
        parts = [p.strip() for p in pathway.split('>')]

        # Check if target pathway appears at correct level
        if level < len(parts):
            if parts[level] == target_pathway:
                return True

    return False


def extract_pathway_genes(mitocarta_df, pathway_csv_file, level):
    """
    Extract genes for each pathway in the pathway CSV file.

    Parameters
    ----------
    mitocarta_df : pd.DataFrame
        MitoCarta data with 'Symbol' and 'MitoCarta3.0_MitoPathways' columns
    pathway_csv_file : Path
        Path to pathway CSV (e.g., mitocarta_pathways_level_0.csv)
    level : int
        Hierarchy level (0, 1, or 2)

    Returns
    -------
    dict
        Dictionary: {pathway_name: [sorted list of gene symbols]}
    """
    # Load pathway CSV to get pathway names
    pathways_df = pd.read_csv(pathway_csv_file)

    # Detect pathway column name (different for Level 2)
    if 'pathway' in pathways_df.columns:
        pathway_col = 'pathway'
    elif 'level2' in pathways_df.columns:
        pathway_col = 'level2'
    elif 'level1' in pathways_df.columns:
        pathway_col = 'level1'
    elif 'level0' in pathways_df.columns:
        pathway_col = 'level0'
    else:
        # Try first column
        pathway_col = pathways_df.columns[0]

    print(f"\nProcessing {len(pathways_df)} pathways from Level {level} (column: '{pathway_col}')...")

    pathway_genes = {}

    for idx, row in pathways_df.iterrows():
        pathway_name = row[pathway_col]
        matching_genes = []

        # Search all MitoCarta genes for this pathway
        for _, mito_row in mitocarta_df.iterrows():
            pathway_str = mito_row['MitoCarta3.0_MitoPathways']

            if pathway_matches(pathway_str, pathway_name, level):
                matching_genes.append(mito_row['Symbol'])

        # Sort genes alphabetically
        pathway_genes[pathway_name] = sorted(matching_genes)

        print(f"  {pathway_name}: {len(matching_genes)} genes")

    return pathway_genes


def save_pathway_gene_lists(pathway_genes, output_dir, level):
    """
    Save gene lists to CSV files.

    Parameters
    ----------
    pathway_genes : dict
        Dictionary: {pathway_name: [gene_list]}
    output_dir : Path
        Directory to save gene list files
    level : int
        Hierarchy level (0, 1, or 2)

    Returns
    -------
    list
        List of saved file paths
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    saved_files = []

    for pathway_name, genes in pathway_genes.items():
        # Create safe filename
        safe_name = pathway_name.replace(' ', '_').replace(',', '').replace('/', '_')
        filename = f"level_{level}_{safe_name}_genes.csv"
        filepath = output_dir / filename

        # Create DataFrame with gene list
        df = pd.DataFrame({
            'gene_symbol': genes,
            'pathway': pathway_name,
            'level': f'level_{level}',
            'n_genes': len(genes)
        })

        # Save to CSV
        df.to_csv(filepath, index=False)
        saved_files.append(filepath)

    print(f"\n  Saved {len(saved_files)} gene list files for Level {level}")

    return saved_files


def main():
    """Main execution function."""

    # Define paths
    base_dir = Path(__file__).parent.parent
    mitocarta_file = base_dir / "data" / "mitocarta.csv"
    tables_dir = base_dir / "results" / "tables"
    output_dir = base_dir / "results" / "tables" / "pathway_genes"

    print("="*70)
    print("MitoCarta Pathway Gene Extraction")
    print("="*70)
    print(f"\nInput: {mitocarta_file}")
    print(f"Output: {output_dir}")

    # Load MitoCarta data
    print("\nLoading MitoCarta 3.0 data...")
    mitocarta_df = pd.read_csv(mitocarta_file)
    print(f"Loaded {len(mitocarta_df)} mitochondrial genes")

    # Check required columns
    if 'Symbol' not in mitocarta_df.columns or 'MitoCarta3.0_MitoPathways' not in mitocarta_df.columns:
        print("ERROR: Required columns not found in MitoCarta CSV")
        print(f"Available columns: {mitocarta_df.columns.tolist()}")
        return

    # Process each level
    pathway_files = {
        0: tables_dir / "mitocarta_pathways_level_0.csv",
        1: tables_dir / "mitocarta_pathways_level_1.csv",
        2: tables_dir / "mitocarta_pathways_level_2_detailed.csv"
    }

    total_files = 0
    all_pathway_genes = {}

    for level, pathway_file in pathway_files.items():
        print(f"\n{'='*70}")
        print(f"LEVEL {level}: {pathway_file.name}")
        print(f"{'='*70}")

        # Extract genes for this level
        pathway_genes = extract_pathway_genes(mitocarta_df, pathway_file, level)

        # Save gene lists
        saved_files = save_pathway_gene_lists(pathway_genes, output_dir, level)
        total_files += len(saved_files)

        # Store for summary
        all_pathway_genes[f'level_{level}'] = pathway_genes

    # Print summary
    print(f"\n{'='*70}")
    print("Extraction Complete!")
    print(f"{'='*70}")
    print(f"\nTotal gene list files created: {total_files}")
    print(f"\nBreakdown by level:")
    for level in range(3):
        n_pathways = len(all_pathway_genes[f'level_{level}'])
        print(f"  Level {level}: {n_pathways} pathways")

    # Validation: check total gene counts
    print(f"\n{'='*70}")
    print("Validation Checks")
    print(f"{'='*70}")

    for level in range(3):
        level_pathways = all_pathway_genes[f'level_{level}']

        # Count unique genes across all pathways at this level
        all_genes_at_level = set()
        for genes in level_pathways.values():
            all_genes_at_level.update(genes)

        print(f"\nLevel {level}:")
        print(f"  Total pathways: {len(level_pathways)}")
        print(f"  Unique genes across all pathways: {len(all_genes_at_level)}")

        # Show pathway with most/fewest genes
        pathway_sizes = {k: len(v) for k, v in level_pathways.items()}
        if pathway_sizes:
            max_pathway = max(pathway_sizes, key=pathway_sizes.get)
            min_pathway = min(pathway_sizes, key=pathway_sizes.get)
            print(f"  Largest pathway: {max_pathway} ({pathway_sizes[max_pathway]} genes)")
            print(f"  Smallest pathway: {min_pathway} ({pathway_sizes[min_pathway]} genes)")

    print(f"\n{'='*70}")
    print(f"All gene lists saved to: {output_dir}/")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
