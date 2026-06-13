# Blood Transcriptome Bridges Orthogonal Immunotypes and Gut Microbiome Enterotypes in Humans

This repository contains the analysis code, source data, and generated results for the study
**"Blood transcriptome bridges orthogonal immunotypes and gut microbiome enterotypes in humans,"**
accepted for publication in *Cell Reports*.

The workflow integrates clinical parameters, gut microbial abundances, blood immune-cell
composition, whole-blood RNA-seq, and bulk CD4+/CD8+ T cell receptor sequencing. Patient
similarities are combined with similarity network fusion (SNF), then analyzed with spectral
clustering, regression models, one-vs-rest statistical testing, enrichment analyses, and downstream
visualizations.

## Pipeline Overview

The analysis is organized around the main notebook and reusable Python/R helper scripts:

1. **Metadata processing:** Loads clinical metadata and harmonizes donor identifiers, cohort labels,
   CMV IgG, BMI, vitamin D, smoking, anxiety, depression, and long-COVID fields. Metadata are used
   for regression and cluster characterization, not for building patient similarities.
2. **TCR-seq processing:** Loads MiXCR clone files and TCR summary statistics from bulk CD4+ and
   CD8+ repertoires, standardizes TCR V/J gene names, aggregates clone-level data, and incorporates
   module or specificity annotations.
3. **CyTOF processing:** Loads CyTOF/FlowSOM metacluster summaries and maps metacluster IDs to
   curated immune-cell labels and metafeatures.
4. **Microbiome processing:** Loads 16S abundance, sample, and taxonomy tables, filters sparse taxa,
   builds relative/absolute abundance matrices, and summarizes eigentaxa by PCA where used.
5. **RNA-seq processing:** Merges readcount files, normalizes counts through DESeq2 via `rpy2`,
   maps Ensembl IDs to gene names, and aggregates genes into BloodGen3 Blood Transcriptome Modules
   using the first principal component per module.
6. **Unimodal analysis:** Runs clustering, eigengap plots, cluster heatmaps, logistic regressions,
   one-vs-rest tests, microbiome alpha-diversity plots, and GSEA where applicable.
7. **SNF fusion analysis:** Combines modality-specific similarity matrices for multimodal analyses
   such as TCR/RNA, CyTOF/RNA, CyTOF/TCR, microbiome/RNA, and trimodal fusion. The manuscript
   interprets RNA-seq as the bridge between largely orthogonal immunotype and enterotype structure.
8. **Network and summary visualizations:** Builds graph visualizations, chord diagrams, bubble
   plots, Venn plots, PCoA plots, spider plots, and validation plots from processed tables.

## Repository Structure

| Path | Purpose |
|------|---------|
| `activ_covid-and-omics.ipynb` | Main Python analysis notebook. |
| `microbiome_validation.ipynb` | Microbiome validation workflow for extra samples. |
| `src/modules/` | Python helpers for data loading, DESeq2, clustering, testing, plotting, GSEA, and network generation. |
| `R_scripts/` | R scripts for reviewer/figure visualizations, microbiome plots, CyTOF metafeature plots, chord plots, and sensitivity analyses. |
| `data/` | Raw and processed input data for clinical metadata, RNA-seq, TCR-seq, CyTOF, microbiome, and validation workflows. |
| `data/processed_data/` | Intermediate tables generated or reused by downstream Python and R analyses. |
| `results/` | Generated plots, graph HTML files, GSEA outputs, cluster heatmaps, and modality-specific result folders. |
| `icons/` | Icons used in interactive network plots. |

## Requirements

### Python

The analysis was developed with **Python 3.11.9**. Install the pinned Python dependencies with:

```bash
python -m pip install -r requirements.txt
```

Key Python packages include:

- `numpy==1.24.2`
- `pandas==1.5.3`
- `scipy==1.14.1`
- `statsmodels==0.14.3`
- `scikit-learn==1.5.2`
- `scikit-bio==0.6.2`
- `gseapy==1.1.3`
- `snfpy==0.2.2`
- `rpy2==3.5.16`
- `matplotlib==3.7.1`
- `seaborn==0.11.2`
- `plotly==5.24.1`
- `pyvis==0.3.2`

### R

The R scripts were developed with **R 4.3.3**. Packages used by the scripts include:

- `tidyverse==2.0.0`
- `tidyr==1.3.1`
- `dplyr==1.1.4`
- `reshape2==1.4.4`
- `Matrix==1.6-5`
- `ROCR==1.0-11`
- `vegan==2.6-8`
- `kernlab==0.9-33`
- `multidiffabundance==0.0.1`
- `circlize==0.4.16`
- `ggplot2==3.5.1`
- `ggvenn==0.1.10`
- `RColorBrewer==1.1-3`
- `gridExtra==2.3`

## How to Run

The main workflow is notebook-first:

```bash
jupyter notebook activ_covid-and-omics.ipynb
```

Run the notebook sections in order. The notebook sets the core data and output directories, loads
the Python helper modules in `src/modules/`, writes processed tables under `data/processed_data/`,
and saves modality-specific outputs under `results/`.

The microbiome validation workflow is separate:

```bash
jupyter notebook microbiome_validation.ipynb
```

R scripts in `R_scripts/` generate additional plots and validation outputs from the processed data.
Several scripts resolve the repository root automatically when run with `Rscript --file`, from
RStudio with the script open, or from an existing repository-root working directory.

```r
source("R_scripts/adonis.R")
source("R_scripts/sensitivity_aitchison.R")
source("R_scripts/chord.R")
source("R_scripts/bubble_plotting.R")
```

## Input Data

The repository includes inputs for the analysis:

| Input area | Location |
|------------|----------|
| Clinical metadata | `data/ERC_meta.xlsx` and validation metadata CSV/XLSX files |
| RNA-seq readcounts | `data/RNA_seq/readcounts_*.txt` |
| TCR clone and summary files | `data/TCR_seq/run*/` and `data/TCR_seq/TCR_stats.txt` |
| CyTOF FlowSOM and cell-type mapping | `data/CyTOF/FlowSOM/` and `data/CyTOF/CellType_mapping/` |
| Microbiome abundance, sample, and taxonomy tables | `data/microbiome_analysis/` |
| Extra microbiome validation data | `data/Extra_microbiome_samples/` |
| Processed handoff tables for R plots | `data/processed_data/` |

The accepted manuscript reports 394 participants with at least one data layer. Per-modality
discovery analyses include whole-blood RNA-seq, bulk TCR-seq, CyTOF mass cytometry, and 16S gut
microbiome data; several validation scripts use independent microbiome or flow-cytometry cohorts.

The code assumes the existing file names and folder layout. Renaming data files or moving generated
tables will require corresponding notebook/script updates.

## Outputs

Primary generated outputs are written to:

| Output area | Location |
|-------------|----------|
| Processed multimodal tables | `data/processed_data/` |
| Modality-specific clustering, regression, GSEA, network, and heatmap outputs | `results/*/` |
| R-script figures and validation plots | `results/R_scripts_plots/` |
| Interactive graph visualizations | `results/*/graph.html` |
| Microbiome significance tables | `results/*/sig_microbiome.csv` |

Many generated outputs are tracked in the repository. Before making source-code changes, check
`git status` and keep regenerated notebooks, processed data, and figures separate from code or
documentation edits.
