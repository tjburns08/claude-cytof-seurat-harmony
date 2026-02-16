# Batch Effect Detection and Correction in CyTOF Data Using Seurat and Harmony

Can we detect and correct batch effects in mass cytometry (CyTOF) data using Seurat, which was built for single-cell RNA-seq? This experiment used Claude Code to test the approach on a multi-center CyTOF dataset using Seurat for data management and Harmony for batch correction, with k-nearest neighbor (kNN) entropy and per-cluster entropy as the quantitative metrics.

## What's here

- **`final_report.Rmd`** — Narrative report with figures and discussion (no code). Start here.
- **`analysis_report.Rmd`** — Full reproducible analysis with code. Requires the data (see below).
- **`analysis.R`** — Standalone R script that runs the complete pipeline.
- **`reasoning_traces.md`** — Decision log explaining every methodological choice.
- **`output/`** — Pre-computed figures (PNGs and PDFs).

## Data

The analysis uses 36 FCS files from a multi-center CyTOF study (6 centers, 3 conditions, 2 replicates per condition per center). The dataset comes from Leipold et al, J. Immunol. Methods 2018, the full text of which can be found [here](https://pmc.ncbi.nlm.nih.gov/articles/PMC5805584/) To reproduce:

1. Place FCS files (from the paper above) in a `data/` directory at the project root
2. File names are expected in the format `Center_N_Condition_Replicate_....fcs`
3. Run `analysis.R` or knit `analysis_report.Rmd`

## Key findings

- Batch effects are detectable using per-cell kNN entropy and a simple UMAP visualization.
- Harmony partially corrected the batch effects (mean kNN entropy 58.4% → 70.9% of theoretical max)
- Biological signal preservation was not assessed, an important gap for production use that must be addressed in subsequent experiments
- See the full discussion in `final_report.Rmd` and the code in `analysis_report.Rmd`

## Dependencies

```r
install.packages(c("Seurat", "harmony", "ggplot2", "patchwork"))
BiocManager::install("flowCore")
```

## Author

Tyler Burns and Claude Code + Opus 4.6
