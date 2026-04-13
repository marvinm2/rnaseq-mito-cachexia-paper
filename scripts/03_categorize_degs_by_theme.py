#!/usr/bin/env python3
"""
Script 07: Comprehensive Three-Theme Pathway Analysis
Identifies genes and pathways for:
1. Proteolysis (E3 ligases, autophagy, proteasome)
2. Myogenesis (muscle regeneration, myogenic factors)
3. Mitochondrial function (already done)

Gene Sources:
- Proteolysis: Literature-based gene lists
  - Fbxo32/Trim63: Bodine et al. (2001) Science 294:1704-1708
  - Autophagy: Mizushima & Komatsu (2011) Cell 147:728-741
  - Mitophagy: Youle & Narendra (2011) Nat Rev Mol Cell Biol 12:9-14
- Myogenesis: Literature-based gene lists
  - MRFs: Buckingham & Rigby (2014) Dev Cell 28:225-238
  - Pax7: Relaix & Zammit (2012) Development 139:2845-2856
  - MEF2: Potthoff & Olson (2007) Development 134:4131-4140

See docs/GENE_SOURCES.md for complete citations and traceability.

Author: Analysis pipeline
Date: 2025-10-06
Updated: 2025-11-26 (added gene source citations)
"""

import pandas as pd
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Setup paths
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results" / "tables"

def get_proteolysis_genes():
    """
    Define proteolysis-related gene sets
    Covers: E3 ubiquitin ligases, autophagy, proteasome
    """
    proteolysis_genes = {
        # Muscle-specific E3 ubiquitin ligases (atrogenes)
        'E3_ligases_muscle': [
            'Fbxo32',  # Atrogin-1/MAFbx - THE atrogene
            'Trim63',  # MuRF1 - THE other atrogene
            'Fbxo30',  # MUSA1
            'Fbxo31',  # MUSA2
            'Nedd4',   # NEDD4
            'Itch',    # ITCH
            'Mdm2',    # MDM2
            'Murf2',   # MuRF2
            'Murf3',   # MuRF3
            'Trim32',  # TRIM32
        ],

        # General E3 ligases
        'E3_ligases_general': [
            'Ube2d1', 'Ube2d2', 'Ube2d3',  # E2 ubiquitin-conjugating enzymes
            'Ube3a', 'Ube3b', 'Ube3c',  # E3 ligases
            'Rnf1', 'Rnf2', 'Rnf4',  # RING finger proteins
            'Cbl', 'Cblb', 'Cblc',  # Cbl family
        ],

        # Autophagy genes
        'autophagy': [
            'Atg5', 'Atg7', 'Atg12', 'Atg16l1',  # ATG conjugation system
            'Becn1',  # Beclin 1 - autophagy initiator
            'Map1lc3a', 'Map1lc3b',  # LC3 - autophagosome marker
            'Ulk1', 'Ulk2',  # ULK kinases
            'Gabarap', 'Gabarapl1', 'Gabarapl2',  # GABARAP family
            'Sqstm1',  # p62 - autophagy cargo receptor
            'Bnip3', 'Bnip3l',  # Mitophagy receptors
            'Fundc1',  # Mitophagy receptor
            'Pink1', 'Prkn',  # PINK1-Parkin pathway
        ],

        # Proteasome subunits
        'proteasome': [
            'Psma1', 'Psma2', 'Psma3', 'Psma4', 'Psma5', 'Psma6', 'Psma7',  # 20S alpha
            'Psmb1', 'Psmb2', 'Psmb3', 'Psmb4', 'Psmb5', 'Psmb6', 'Psmb7',  # 20S beta
            'Psmb8', 'Psmb9', 'Psmb10',  # Immunoproteasome
            'Psmc1', 'Psmc2', 'Psmc3', 'Psmc4', 'Psmc5', 'Psmc6',  # 19S ATPase
            'Psmd1', 'Psmd2', 'Psmd3', 'Psmd4', 'Psmd11', 'Psmd12', 'Psmd13',  # 19S non-ATPase
        ],

        # Calpains (calcium-dependent proteases)
        'calpains': [
            'Capn1', 'Capn2', 'Capn3',  # Calpains
            'Capns1', 'Capns2',  # Calpain small subunits
            'Cast',  # Calpastatin (inhibitor)
        ],

        # Caspases (apoptotic proteases)
        'caspases': [
            'Casp3', 'Casp6', 'Casp7',  # Executioner caspases
            'Casp8', 'Casp9', 'Casp10',  # Initiator caspases
        ],
    }

    # Flatten into single list
    all_proteolysis = []
    for category, genes in proteolysis_genes.items():
        all_proteolysis.extend(genes)

    return proteolysis_genes, list(set(all_proteolysis))

def get_myogenesis_genes():
    """
    Define myogenesis/muscle regeneration gene sets
    Covers: Myogenic factors, muscle structure, satellite cells
    """
    myogenesis_genes = {
        # Myogenic regulatory factors (MRFs)
        'myogenic_factors': [
            'Myod1',  # MyoD - master regulator
            'Myog',   # Myogenin - differentiation
            'Myf5',   # Myf5 - specification
            'Myf6',   # MRF4 - differentiation
            'Pax3', 'Pax7',  # Paired box TFs - satellite cells
        ],

        # Myosin heavy chains
        'myosins_heavy': [
            'Myh1',  # Type IIX
            'Myh2',  # Type IIA
            'Myh3',  # Embryonic
            'Myh4',  # Type IIB
            'Myh7',  # Type I (slow)
            'Myh8',  # Perinatal
        ],

        # Myosin light chains
        'myosins_light': [
            'Myl1', 'Myl2', 'Myl3', 'Myl4', 'Myl6',
        ],

        # Actin and actin-binding proteins
        'actin_cytoskeleton': [
            'Acta1',  # Skeletal muscle actin
            'Actc1',  # Cardiac actin
            'Actg1',  # Cytoplasmic gamma actin
            'Tpm1', 'Tpm2', 'Tpm3',  # Tropomyosins
            'Tnni1', 'Tnni2', 'Tnni3',  # Troponin I
            'Tnnt1', 'Tnnt2', 'Tnnt3',  # Troponin T
            'Tnnc1', 'Tnnc2',  # Troponin C
        ],

        # Sarcomere and structural proteins
        'sarcomere': [
            'Ttn',  # Titin
            'Neb',  # Nebulin
            'Des',  # Desmin
            'Vim',  # Vimentin
            'Myom1', 'Myom2', 'Myom3',  # Myomesin
        ],

        # Satellite cell markers
        'satellite_cells': [
            'Pax7',  # Satellite cell marker
            'Myf5',  # Activated satellite cells
            'Cd34',  # Satellite cell surface marker
            'Vcam1',  # Satellite cell marker
            'Itga7',  # Alpha7 integrin
            'Cxcr4',  # Chemokine receptor
        ],

        # Growth factors and signaling
        'growth_signaling': [
            'Igf1', 'Igf2',  # IGF signaling
            'Mstn',  # Myostatin (negative regulator)
            'Gdf11',  # GDF11 (myostatin-like)
            'Fst',  # Follistatin (myostatin inhibitor)
            'Hgf',  # Hepatocyte growth factor
        ],

        # Cell cycle regulators (for proliferation)
        'cell_cycle': [
            'Cdkn1a', 'Cdkn1b', 'Cdkn1c',  # p21, p27, p57
            'Ccnd1', 'Ccnd2', 'Ccnd3',  # Cyclin D
            'Ccne1', 'Ccne2',  # Cyclin E
        ],
    }

    # Flatten into single list
    all_myogenesis = []
    for category, genes in myogenesis_genes.items():
        all_myogenesis.extend(genes)

    return myogenesis_genes, list(set(all_myogenesis))

def identify_deg_themes(all_degs):
    """
    Categorize DEGs into the three themes
    """
    print("\n" + "="*70)
    print("CATEGORIZING DEGs BY BIOLOGICAL THEME")
    print("="*70)

    # Get gene lists
    proteolysis_categories, proteolysis_list = get_proteolysis_genes()
    myogenesis_categories, myogenesis_list = get_myogenesis_genes()

    # Load mitochondrial genes
    mito_degs = pd.read_csv(RESULTS_DIR / "mitochondrial_degs.csv")
    mito_list = mito_degs['GeneSymbol'].str.upper().tolist()

    # Normalize gene symbols
    all_degs['GeneSymbol_upper'] = all_degs['GeneSymbol'].str.upper()
    proteolysis_set = set([g.upper() for g in proteolysis_list])
    myogenesis_set = set([g.upper() for g in myogenesis_list])
    mito_set = set(mito_list)

    # Categorize
    all_degs['theme_proteolysis'] = all_degs['GeneSymbol_upper'].isin(proteolysis_set)
    all_degs['theme_myogenesis'] = all_degs['GeneSymbol_upper'].isin(myogenesis_set)
    all_degs['theme_mitochondrial'] = all_degs['GeneSymbol_upper'].isin(mito_set)

    # Count overlaps and exclusives
    proteolysis_degs = all_degs[all_degs['theme_proteolysis']]
    myogenesis_degs = all_degs[all_degs['theme_myogenesis']]
    mitochondrial_degs = all_degs[all_degs['theme_mitochondrial']]

    print(f"\nTheme Distribution:")
    print(f"  Proteolysis genes: {len(proteolysis_degs)} / {len(proteolysis_set)} ({len(proteolysis_degs)/len(all_degs)*100:.1f}% of DEGs)")
    print(f"  Myogenesis genes: {len(myogenesis_degs)} / {len(myogenesis_set)} ({len(myogenesis_degs)/len(all_degs)*100:.1f}% of DEGs)")
    print(f"  Mitochondrial genes: {len(mitochondrial_degs)} / {len(mito_set)} ({len(mitochondrial_degs)/len(all_degs)*100:.1f}% of DEGs)")

    # Check overlaps
    overlaps = all_degs[
        (all_degs['theme_proteolysis'] | all_degs['theme_myogenesis'] | all_degs['theme_mitochondrial'])
    ]
    print(f"\nTotal DEGs in at least one theme: {len(overlaps)} ({len(overlaps)/len(all_degs)*100:.1f}%)")

    # Regulation direction by theme
    for theme_name, theme_col in [('Proteolysis', 'theme_proteolysis'),
                                   ('Myogenesis', 'theme_myogenesis'),
                                   ('Mitochondrial', 'theme_mitochondrial')]:
        theme_degs = all_degs[all_degs[theme_col]]
        if len(theme_degs) > 0:
            n_up = (theme_degs['regulation'] == 'up').sum()
            n_down = (theme_degs['regulation'] == 'down').sum()
            print(f"\n{theme_name}:")
            print(f"  Upregulated: {n_up} ({n_up/len(theme_degs)*100:.1f}%)")
            print(f"  Downregulated: {n_down} ({n_down/len(theme_degs)*100:.1f}%)")

    return proteolysis_degs, myogenesis_degs, proteolysis_categories, myogenesis_categories

def save_theme_gene_lists(proteolysis_degs, myogenesis_degs):
    """Save gene lists for each theme"""
    print("\n" + "="*70)
    print("SAVING THEME-SPECIFIC GENE LISTS")
    print("="*70)

    # Save proteolysis genes
    prot_file = RESULTS_DIR / "proteolysis_degs.csv"
    proteolysis_degs.to_csv(prot_file, index=False)
    print(f"\nProteolysis DEGs: {prot_file}")
    print(f"  Total: {len(proteolysis_degs)}")

    # Save myogenesis genes
    myo_file = RESULTS_DIR / "myogenesis_degs.csv"
    myogenesis_degs.to_csv(myo_file, index=False)
    print(f"\nMyogenesis DEGs: {myo_file}")
    print(f"  Total: {len(myogenesis_degs)}")

    # Save gene symbols for g:Profiler
    all_deg_symbols_file = RESULTS_DIR / "all_degs_for_gprofiler.txt"
    all_degs = pd.read_csv(RESULTS_DIR / "mouse_degs_filtered.csv")

    # Handle column name issues - gene symbols are in Unnamed: 2
    if 'GeneSymbol' not in all_degs.columns:
        if 'Unnamed: 2' in all_degs.columns:
            all_degs['GeneSymbol'] = all_degs['Unnamed: 2']
        elif 'Unnamed: 1' in all_degs.columns:
            all_degs['GeneSymbol'] = all_degs['Unnamed: 1']

    # Get unique gene symbols
    gene_symbols = all_degs['GeneSymbol'].dropna().unique()
    with open(all_deg_symbols_file, 'w') as f:
        f.write('\n'.join(gene_symbols))

    print(f"\nAll DEGs for g:Profiler: {all_deg_symbols_file}")
    print(f"  Total unique genes: {len(gene_symbols)}")
    print("\nInstructions:")
    print("1. Go to https://biit.cs.ut.ee/gprofiler/gost")
    print("2. Copy gene list from file above")
    print("3. Select organism: Mus musculus")
    print("4. Data sources: GO:BP, GO:MF, GO:CC, KEGG, REAC, WikiPathways")
    print("5. Run analysis and download results")

def create_summary_statistics(proteolysis_degs, myogenesis_degs, proteolysis_categories, myogenesis_categories):
    """Create summary statistics by subcategory"""
    print("\n" + "="*70)
    print("THEME SUBCATEGORY BREAKDOWN")
    print("="*70)

    # Proteolysis subcategories
    print("\nProteolysis subcategories:")
    for category, genes in proteolysis_categories.items():
        genes_upper = [g.upper() for g in genes]
        found = proteolysis_degs[proteolysis_degs['GeneSymbol'].str.upper().isin(genes_upper)]
        if len(found) > 0:
            n_up = (found['regulation'] == 'up').sum()
            n_down = (found['regulation'] == 'down').sum()
            print(f"  {category}: {len(found)}/{len(genes)} DEGs ({n_down}↓, {n_up}↑)")
            if len(found) <= 10:  # Show genes if small list
                print(f"    Genes: {', '.join(found['GeneSymbol'].tolist())}")

    # Myogenesis subcategories
    print("\nMyogenesis subcategories:")
    for category, genes in myogenesis_categories.items():
        genes_upper = [g.upper() for g in genes]
        found = myogenesis_degs[myogenesis_degs['GeneSymbol'].str.upper().isin(genes_upper)]
        if len(found) > 0:
            n_up = (found['regulation'] == 'up').sum()
            n_down = (found['regulation'] == 'down').sum()
            print(f"  {category}: {len(found)}/{len(genes)} DEGs ({n_down}↓, {n_up}↑)")
            if len(found) <= 10:  # Show genes if small list
                print(f"    Genes: {', '.join(found['GeneSymbol'].tolist())}")

def main():
    """Main pipeline"""
    print("\n" + "="*70)
    print("COMPREHENSIVE THREE-THEME CACHEXIA ANALYSIS")
    print("="*70)
    print("\nThemes to analyze:")
    print("1. Proteolysis (E3 ligases, autophagy, proteasome)")
    print("2. Myogenesis (muscle regeneration, myogenic factors)")
    print("3. Mitochondrial function (already analyzed)")

    # Load all DEGs
    print("\nLoading DEG data...")
    all_degs = pd.read_csv(RESULTS_DIR / "mouse_degs_filtered.csv")

    # Handle column name issues - gene symbols are in Unnamed: 2
    if 'GeneSymbol' not in all_degs.columns:
        if 'Unnamed: 2' in all_degs.columns:
            all_degs['GeneSymbol'] = all_degs['Unnamed: 2']
        elif 'Unnamed: 1' in all_degs.columns:
            all_degs['GeneSymbol'] = all_degs['Unnamed: 1']

    print(f"Total DEGs: {len(all_degs)}")

    # Identify genes in each theme
    proteolysis_degs, myogenesis_degs, prot_cats, myo_cats = identify_deg_themes(all_degs)

    # Save results
    save_theme_gene_lists(proteolysis_degs, myogenesis_degs)

    # Summary statistics
    create_summary_statistics(proteolysis_degs, myogenesis_degs, prot_cats, myo_cats)

    print("\n" + "="*70)
    print("COMPREHENSIVE ANALYSIS COMPLETE")
    print("="*70)
    print("\nNext steps:")
    print("1. Submit all_degs_for_gprofiler.txt to g:Profiler web interface")
    print("2. Download g:Profiler results as CSV")
    print("3. Run script 08 to generate comprehensive figure")

if __name__ == "__main__":
    main()
