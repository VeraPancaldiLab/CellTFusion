# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Package Overview

CellTFusion is an R package for integrating immune-cell type
deconvolution with transcription factor (TF)–gene regulatory networks to
characterize immune cell states in the tumor microenvironment using bulk
RNAseq data.

## Development Commands

This is a standard R package using the devtools workflow:

``` r

devtools::load_all()      # Load package in development
devtools::document()      # Regenerate NAMESPACE and Rd files from roxygen2
devtools::check()         # Full package check (CRAN-style)
devtools::test()          # Run unit tests (testthat)
devtools::install()       # Install package locally
devtools::build()         # Build source tarball
```

Run a single test file:

``` r

testthat::test_file("tests/testthat/test-<name>.R")
```

Launch the Shiny app:

``` r

CellTFusion::launch_app()
```

Documentation uses **roxygen2** with Markdown enabled
(`Roxygen: list(markdown = TRUE)`). Always run
[`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
after modifying roxygen comments.

## Architecture

### Pipeline Flow

The package implements a multi-step pipeline:

1.  **Normalization** — log-TPM normalization of raw counts
2.  **Cell-type deconvolution** — multiple algorithms via `multideconv`
    (Quantiseq, Epidish, DeconRNASeq, DWLS, CIBERSORTx)
3.  **TF activity inference** — `viper` + CollecTRI/Dorothea regulons
    via `decoupleR`
4.  **TF module construction** — WGCNA-based weighted co-expression
    networks
5.  **Pathway scoring** — PROGENy-based pathway activities
6.  **Cell group identification** — supervised and unsupervised
    clustering
7.  **Statistical analysis** — clinical trait associations, enrichment
    analyses

### Key Source Files

- `R/CellTFusion.R` — All core logic (~4,500 lines); single monolithic
  file containing every exported function
- `inst/shiny/server.R` / `ui.R` — Shiny app backend and frontend
- `vignettes/CellTFusion.Rmd` — Primary usage tutorial

### Main Exported Functions

| Function | Purpose |
|----|----|
| [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md) | Full pipeline wrapper |
| [`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md) | TF activity scoring via viper |
| [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md) | WGCNA module construction |
| [`compute.pathway.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md) | PROGENy pathway scoring |
| [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md) | Unsupervised/supervised cell group clustering |
| [`compute.test.set()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.test.set.md) | Apply trained model to new dataset |
| [`identify_hub_TFs()`](https://verapancaldilab.github.io/CellTFusion/reference/identify_hub_TFs.md) | Identify driver TFs from modules |
| [`compute_factor_gsea()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_factor_gsea.md) | GSEA on latent factors |
| [`compute.metadata.association()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.metadata.association.md) | Clinical trait association + visualization |
| [`compute.modules.relationship()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md) | Correlate TF modules with pathways |
| [`compute.modules.enrichment()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.enrichment.md) | Pathway enrichment of TF modules |

### Tutorial Data (in `data/`)

Pre-built `.rda` objects for examples and testing: `raw.counts.tuto`,
`counts.norm.tuto`, `traitdata.tuto`, `tfs.tuto`, `network.tuto`,
`deconv.tuto`, `deconv_subgroups.tuto`.

## Dependencies

Heavy dependencies are in `Imports` (always loaded): `multideconv`,
`decoupleR`, `viper`, `WGCNA`, `clusterProfiler`, `limma`, `ggplot2`,
and many tidyverse packages.

ML classifiers (`caret`, `C50`, `glmnet`, `randomForest`, `xgboost`) are
in `Suggests` (optional, used in supervised analysis).

## CI/CD

GitHub Actions (`.github/workflows/pkgdown.yaml`) builds and deploys the
pkgdown documentation site to `gh-pages` on push to `main`. The pkgdown
site config is in `_pkgdown.yml`.
