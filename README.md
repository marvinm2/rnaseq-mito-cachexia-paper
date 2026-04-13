# Mitochondrial impairment in cachexia — RNA-seq figure pipeline

Scripts to reproduce the RNA-seq figures for the cachexia / mitochondrial dysfunction manuscript. Mouse model: 344P tumor-bearing (cachexia) vs sham, bulk RNA-seq, DEGs at padj < 0.01 and |log2FC| ≥ 0.5.

## Output

Running the pipeline end-to-end writes these files to `results/figures/` (PDF + PNG + SVG for each):

- `figure_go_myogenesis_proteolysis` — combined myogenesis + proteolysis GO panels (main figure)
- `figure_volcano_mitochondrial` — genome-wide volcano with mitochondrial genes highlighted
- `figure_oxphos_hierarchy` — OXPHOS three-level hierarchy
- `figure_myogenesis_network` — TF-target network of 56 literature-curated myogenesis genes
- `figure_appendix_mitocarta_hierarchy` — 3-level MitoCarta hierarchy (appendix)

The composite manuscript figure combining the volcano, OXPHOS panel and a PathVisio pathway is assembled outside this repo (Illustrator/Inkscape). The PathVisio pathway panel is not generated here.

## Install

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Fetch inputs

The repo ships no data. Populate `data/` with:

```bash
# 1. Mouse DEG file (raw DESeq2 output, tab-separated, 11,599 genes)
#    Place at: data/mouse_degs.txt
#    Source: internal — not publicly hosted; contact the authors.

# 2. MitoCarta 3.0 (Mouse)
curl -L -o data/Mouse.MitoCarta3.0.xls \
  https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Mouse.MitoCarta3.0.xls
# Convert the 'A Mouse MitoCarta3.0' sheet to CSV and save as data/mitocarta.csv
# (the scripts read mitocarta.csv, not the .xls directly).

# 3. GO ontology
mkdir -p data/go_ontology
curl -L -o data/go_ontology/go-basic.obo \
  http://purl.obolibrary.org/obo/go/go-basic.obo

# 4. TFLink mouse interactions (needed for the myogenesis network figure)
mkdir -p data/tflink_mice.tsv
curl -L -o /tmp/tflink_mice.tsv.gz \
  https://cdn.netbiol.org/tflink/download_files/TFLink_Mus_musculus_interactions_All_simpleFormat_v1.0.tsv.gz
gunzip -c /tmp/tflink_mice.tsv.gz \
  > data/tflink_mice.tsv/TFLink_Mus_musculus_interactions_All_simpleFormat_v1.0.tsv
```

## Run

Execute in order — each script writes into `results/tables/` or `results/figures/`:

```bash
# Stage 1: data prep
python scripts/01_filter_degs.py
python scripts/02_annotate_mitocarta.py
python scripts/03_categorize_degs_by_theme.py

# Stage 2: statistics
python scripts/04_mitocarta_pathway_hierarchy.py
python scripts/05_mitocarta_pathway_gene_lists.py
python scripts/06_mitocarta_pathway_stats.py
python scripts/07_oxphos_gene_lists.py
python scripts/08_oxphos_stats.py
python scripts/09_go_fetch_myogenesis_proteolysis.py
python scripts/10_go_stats.py

# Stage 3: figures
python scripts/11_figure_main_go_panels.py          # figure_go_myogenesis_proteolysis + figure_volcano_mitochondrial
python scripts/12_figure_oxphos.py                  # figure_oxphos_hierarchy
python scripts/13_figure_myogenesis_network.py      # figure_myogenesis_network
python scripts/14_figure_mitocarta_three_levels.py  # figure_appendix_mitocarta_hierarchy
```

Script 09 calls the QuickGO REST API; script 04's g:Profiler step also hits a public endpoint. Both need network access.

## Repo layout

```
scripts/     Analysis + figure scripts (run in the order above)
gene_sets/   literature_56_myogenesis_genes.csv — hand-curated gene list
data/        .gitignored — populated by the fetch step above
results/     .gitignored — populated by running the scripts
```
