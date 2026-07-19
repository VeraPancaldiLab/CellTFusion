# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

CellTFusion is an R package for integrating immune-cell type deconvolution with transcription factor (TF)–gene regulatory networks to characterize immune cell states in the tumor microenvironment using bulk RNAseq data.

## Development Commands

This is a standard R package using the devtools workflow:

```r
devtools::load_all()      # Load package in development
devtools::document()      # Regenerate NAMESPACE and Rd files from roxygen2
devtools::check()         # Full package check (CRAN-style)
devtools::test()          # Run unit tests (testthat)
devtools::install()       # Install package locally
devtools::build()         # Build source tarball
```

There is currently no `tests/testthat/` directory in the repo (no unit tests exist yet), even though `testthat` (edition 3) is configured in `DESCRIPTION` and listed in `Suggests`. Once tests are added, run a single file with:
```r
testthat::test_file("tests/testthat/test-<name>.R")
```

Launch the Shiny app (there is no exported `launch_app()` wrapper — run it directly):
```r
shiny::runApp(system.file("shiny", package = "CellTFusion"))
# or, from a source checkout:
shiny::runApp("inst/shiny")
```

Documentation uses **roxygen2** with Markdown enabled (`Roxygen: list(markdown = TRUE)`). Always run `devtools::document()` after modifying roxygen comments.

## Architecture

### Pipeline Flow

The package implements a multi-step pipeline (see `vignettes/CellTFusion.Rmd` for the authoritative step/function/article mapping):

1. **Normalization** — log-TPM normalization of raw counts (`ADImpute::NormalizeTPM`)
2. **Cell-type deconvolution** — multiple algorithms via `multideconv::compute.deconvolution()` (Quantiseq, Epidish, DeconRNASeq, DWLS, CIBERSORTx)
3. **TF activity inference** — `compute.TFs.activity()`: CollecTRI/Dorothea/ARACNe regulons scored via `decoupleR` consensus methods
4. **TF module construction** — `compute.WTCNA()`: WGCNA-based weighted co-expression networks
5. **Pathway scoring** — `compute.pathway.activity()`: PROGENy-based pathway activities via `decoupleR::run_mlm()`
6. **Cell group construction** — `construct_cell_groups()`: supervised and unsupervised clustering
7. **Latent factor extraction** — `compute.latent_factors()`: NMF-based latent factors; `compute_cells_niches()` for cell niche derivation
8. **TME state characterisation** — Hallmark GSEA per factor (`compute_factor_gsea()`), meta-program derivation/mapping (`derive_meta_programs()`, `map_factors_to_metaprograms()`), and TME subtype annotation against Bagaev et al. (2021) MFP subtypes (`map_factors_to_TME()`, `annotate_metaprograms_TME()`)
9. **Statistical analysis** — clinical trait association (`scores.stat.analysis()` and the `scores.*` family), survival analysis (`compute.survival.analysis()`)
10. **Test-set / batch projection** — apply a trained model to new data (`compute.test.set()`, `project_test_factors()`) or run multi-cohort analysis via `batch = TRUE` in `CellTFusion()`

### Key Source Files

- `R/CellTFusion.R` — All core logic (~5,050 lines); single monolithic file containing every exported function
- `inst/shiny/server.R` / `ui.R` — Shiny app backend and frontend
- `vignettes/CellTFusion.Rmd` — Getting-started overview with the full pipeline step table
- `vignettes/articles/` — In-depth tutorials: `01-feature-computation`, `02-cell-groups`, `03-tme-states`, `04-analysis`, `06-machine-learning`, `07-batch-analysis`

### Main Exported Functions

| Function | Purpose |
|---|---|
| `CellTFusion()` | Full pipeline wrapper (unsupervised, supervised, and multi-cohort/batch modes) |
| `compute.TFs.activity()` | TF activity scoring via `decoupleR` |
| `compute.WTCNA()` | WGCNA module construction |
| `compute.pathway.activity()` | PROGENy pathway scoring |
| `construct_cell_groups()` | Unsupervised/supervised cell group clustering |
| `compute.latent_factors()` | NMF-based latent factor extraction from cell group scores |
| `compute.test.set()` / `project_test_factors()` | Apply a trained model / project a test set onto trained NMF factors |
| `identify_hub_TFs()` | Identify driver TFs from modules |
| `compute_factor_gsea()` | Hallmark GSEA on latent factors |
| `derive_meta_programs()` / `map_factors_to_metaprograms()` | Derive and map latent factors to TCGA meta-programs |
| `map_factors_to_TME()` / `annotate_metaprograms_TME()` | Annotate factors/meta-programs with Bagaev et al. (2021) TME (MFP) subtypes |
| `compute.metadata.association()` | Clinical trait association + visualization |
| `compute.survival.analysis()` | Survival analysis (S3 method; requires `survival`, `survminer`, `gridExtra`) |
| `compute.modules.relationship()` | Correlate TF modules with pathways |
| `compute.modules.enrichment()` | Pathway enrichment of TF modules |

This is a curated subset — see `NAMESPACE` or the [pkgdown reference](https://verapancaldilab.github.io/CellTFusion/reference/index.html) for the full list, which also includes the `scores.*` statistical-test family and lower-level helpers (`cell.groups.computation()`, `classify.deconvolution()`, `create_tfs_modules()`, etc.).

### Tutorial Data (in `data/`)

Pre-built `.rda` objects for examples and testing: `raw.counts.tuto`, `counts.norm.tuto`, `traitdata.tuto`, `tfs.tuto`, `network.tuto`, `deconv.tuto`, `deconv_subgroups.tuto`.

## Dependencies

Heavy dependencies are in `Imports` (always loaded): `multideconv`, `decoupleR`, `WGCNA`, `limma`, `GSVA`, `ggplot2`, and several tidyverse packages (`dplyr`, `tidyr`, `tibble`, `stringr`, `purrr`).

`viper`, `clusterProfiler`, `dorothea`, `OmnipathR`, and ML classifiers (`caret`, `C50`, `glmnet`, `randomForest`, `xgboost`) are all in `Suggests` (optional, loaded only when the relevant code path is used).

## CI/CD

GitHub Actions (`.github/workflows/pkgdown.yaml`) builds and deploys the pkgdown documentation site to `gh-pages` on push to `main`. The pkgdown site config is in `_pkgdown.yml`.
