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

Unit tests (`testthat`, edition 3) are in `tests/testthat/`: they use
the tutorial data, need no download and take about 5 minutes
(`helper-data.R` builds the shared objects). Run a single file with:

``` r

testthat::test_file("tests/testthat/test-<name>.R")
```

Launch the Shiny app (there is no exported `launch_app()` wrapper — run
it directly):

``` r

shiny::runApp(system.file("shiny", package = "CellTFusion"))
# or, from a source checkout:
shiny::runApp("inst/shiny")
```

Documentation uses **roxygen2** with Markdown enabled
(`Roxygen: list(markdown = TRUE)`). Always run
[`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
after modifying roxygen comments.

## Architecture

### Pipeline Flow

The package implements a multi-step pipeline (see
`vignettes/CellTFusion.Rmd` for the authoritative step/function/tutorial
mapping):

1.  **Normalization** — log-TPM normalization of raw counts
    ([`ADImpute::NormalizeTPM`](https://rdrr.io/pkg/ADImpute/man/NormalizeTPM.html))
2.  **Cell-type deconvolution** — multiple algorithms via
    [`multideconv::compute.deconvolution()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.html)
    (Quantiseq, Epidish, DeconRNASeq, DWLS, CIBERSORTx)
3.  **TF activity inference** —
    [`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md):
    CollecTRI/Dorothea/ARACNe regulons scored via `decoupleR`
    (`method = "ulm"` by default; `"mlm"`, `"wsum"`, `"consensus"`)
4.  **TF module construction** —
    [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md):
    WGCNA-based weighted co-expression networks
5.  **Pathway scoring** —
    [`compute.pathway.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md):
    PROGENy-based pathway activities via `decoupleR` (`method = "ulm"`
    by default)
6.  **Cell group construction** —
    [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md):
    TF module–deconvolution correlation, dendrogram cutting and CCA
    composite scores
7.  **Latent factor extraction** —
    [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md):
    NMF-based latent factors;
    [`compute_cells_niches()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_cells_niches.md)
    for cell niche derivation
8.  **Statistical analysis** — clinical trait association
    ([`scores.stat.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.stat.analysis.md)
    and the `scores.*` family), survival analysis
    ([`compute.survival.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.survival.analysis.md))
9.  **Test-set / batch projection** — apply a trained model to new data
    ([`compute.test.set()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.test.set.md)
    for cell group scores;
    [`project_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/project_factors.md)
    to project those scores onto trained NMF factors), run leakage-aware
    cross-validation with `pipeML`
    ([`prepare_celltfusion_folds()`](https://verapancaldilab.github.io/CellTFusion/reference/prepare_celltfusion_folds.md)),
    or run multi-cohort analysis via `batch` (cohort of each sample) in
    [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)

### Key Source Files

- `R/CellTFusion.R` — All core logic (~4,700 lines); single monolithic
  file containing every exported function
- `inst/shiny/server.R` / `ui.R` — Shiny app backend and frontend
- `vignettes/CellTFusion.Rmd` — Getting-started overview (installation,
  quick start, pipeline step table, tutorial list)
- `vignettes/a1_feature_computation.Rmd` … `a6_batch_analysis.Rmd` —
  In-depth tutorials (feature computation, cell groups, statistical
  analysis, machine learning, multi-cohort; there is no `a3` since the
  TME state functions were removed on 2026-10-16). They are package
  vignettes (same layout as `multideconv`); heavy code chunks use
  `eval = FALSE` and figures live in `vignettes/figures/`. The pkgdown
  article menu is defined in `_pkgdown.yml`.

### Main Exported Functions

| Function | Purpose |
|----|----|
| [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md) | Full pipeline wrapper (single-cohort and multi-cohort/batch modes) |
| [`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md) | TF activity scoring via `decoupleR` |
| [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md) | WGCNA module construction |
| [`compute.pathway.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md) | PROGENy pathway scoring |
| [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md) | Cell group construction and composite scores |
| [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md) | NMF-based latent factor extraction from cell group scores |
| [`compute.test.set()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.test.set.md) / [`project_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/project_factors.md) | Score new samples on trained cell groups / project cell group scores onto trained NMF factors |
| [`prepare_celltfusion_folds()`](https://verapancaldilab.github.io/CellTFusion/reference/prepare_celltfusion_folds.md) | `pipeML` fold construction function (`fold_construction_fun`): runs [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md) within each CV fold (cell groups only, `compute_factors = FALSE`) and scores the test samples with [`compute.test.set()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.test.set.md); the ML features are the cell group scores (classification and survival) |
| [`identify_hub_TFs()`](https://verapancaldilab.github.io/CellTFusion/reference/identify_hub_TFs.md) | Identify driver TFs from modules |
| [`compute.metadata.association()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.metadata.association.md) | Clinical trait association + visualization |
| [`compute.survival.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.survival.analysis.md) | Kaplan-Meier / log-rank survival analysis (requires `survival`, `survminer`, `gridExtra`) |
| [`compute.modules.relationship()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md) | Correlate TF modules with other continuous features (pathways, deconvolution features); returns correlations, p-values and significant pairs |
| [`compute.modules.enrichment()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.enrichment.md) | Pathway enrichment of TF modules |

This is a curated subset — see `NAMESPACE` or the [pkgdown
reference](https://verapancaldilab.github.io/CellTFusion/reference/index.html)
for the full list, which also includes the `scores.*` statistical-test
family and lower-level helpers
([`cell.groups.computation()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.computation.md),
[`compute_composite_score()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_composite_score.md),
etc.).

### Tutorial Data (in `data/`)

Pre-built `.rda` objects for examples and testing: `raw.counts.tuto`,
`counts.norm.tuto`, `traitdata.tuto`, `tfs.tuto`, `network.tuto`,
`deconv.tuto`, `deconv_subgroups.tuto`. `data-raw/tutorial_data.R`
rebuilds the ones that depend on the package defaults (`traitdata.tuto`
with overall survival, `tfs.tuto`, `network.tuto`).

## Dependencies

Heavy dependencies are in `Imports` (always loaded): `multideconv`,
`decoupleR`, `WGCNA`, `GSVA`, `ggplot2`, and several tidyverse packages
(`dplyr`, `tidyr`, `tibble`). `multideconv` functions are imported via
`@importFrom` in `R/CellTFusion.R`.

`OmnipathR` (needed by `decoupleR` to fetch CollecTRI and PROGENy),
`dorothea` (Dorothea TF collection), `RcppML` (NMF),
`ReactomePA`/`enrichplot`
([`compute.modules.enrichment()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.enrichment.md)),
`survival`/`survminer`/`gridExtra` (survival analysis), `caret`
(machine-learning vignette) and `shiny` (the app in `inst/shiny`) are in
`Suggests` (optional, loaded only when the relevant code path is used).

## CI/CD

GitHub Actions (`.github/workflows/pkgdown.yaml`) builds and deploys the
pkgdown documentation site to `gh-pages` on push to `main`. The pkgdown
site config is in `_pkgdown.yml`.
