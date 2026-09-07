# TCGA PanCancer Atlas mechanostress and lipid-droplet analysis
[![DOI](https://zenodo.org/badge/1359819872.svg)](https://doi.org/10.5281/zenodo.22640778)

This repository contains the R code used to generate BRCA
scatter plots and four-cohort correlation heatmaps based on
TCGA PanCancer Atlas 2018 gene-expression data.

## Analysis code

The complete analysis is implemented in:

`analysis/TCGA_BRCA_scatter_4cohort_heatmap_COMPACT.R`

## Cohorts

- CRC: coadread_tcga_pan_can_atlas_2018
- LUAD: luad_tcga_pan_can_atlas_2018
- BRCA: brca_tcga_pan_can_atlas_2018
- SKCM: skcm_tcga_pan_can_atlas_2018

Heatmap order:

CRC, LUAD, BRCA, SKCM

## Data source

Gene-expression data are downloaded through the cBioPortal API
from TCGA PanCancer Atlas 2018 studies.

All complete samples in each RNA-seq sample list are used.
No tumor-sample filter is applied.

## Scoring

Hierarchical mean Z-score calculation is used for:

- Mechanostress
- LD buffering
- Turnover competence

Turnover-adjusted LD excess is calculated within each cohort as:

Adjusted LD excess =
Z(residuals(lm(LD_buffering ~ Turnover)))

## Required R packages

- httr
- jsonlite
- ggplot2
- writexl
- scales
- renv

## Running the analysis

Open `TCGA_Mechanostress_LD_analysis.Rproj` in RStudio.

Restore the package environment:

```r
install.packages("renv")
renv::restore(prompt = FALSE)
```

Run the complete analysis:

```r
source(
  "analysis/TCGA_BRCA_scatter_4cohort_heatmap_COMPACT.R",
  echo = TRUE
)
```

Results are written automatically to:

`TCGA_compact_scatter_heatmap_results/`

## Statistical analysis

Associations are evaluated using two-sided Spearman correlation.

Three heatmap versions are generated:

1. Color only
2. Spearman rho and raw P-value
3. Spearman rho and raw-P significance stars

## Citation

Lee, Jeong Uk. (2026). *TCGA Mechanostress and Lipid-Droplet Analysis* (Version v1.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.22640779
