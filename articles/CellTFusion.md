# Getting Started with CellTFusion

``` r

library(CellTFusion)
#> 
#> 
```

## Overview

`CellTFusion` integrates immune cell-type deconvolution with
transcription factor (TF)–gene regulatory networks to characterize
immune cell states in the tumor microenvironment from bulk RNA-seq data.

Starting from a raw count matrix, the pipeline produces **latent
factors** — compact representations of the TME landscape that can be
tested for clinical associations, mapped to known cancer meta-programs,
and used as features for machine learning.

## Pipeline steps

| Step | Function | Article |
|----|----|----|
| Cell-type deconvolution | [`multideconv::compute.deconvolution()`](https://rdrr.io/pkg/multideconv/man/compute.deconvolution.html) | [Feature Computation](https://verapancaldilab.github.io/CellTFusion/articles/articles/01-feature-computation.md) |
| TF activity inference | [`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md) | [Feature Computation](https://verapancaldilab.github.io/CellTFusion/articles/articles/01-feature-computation.md) |
| TF module construction | [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md) | [Feature Computation](https://verapancaldilab.github.io/CellTFusion/articles/articles/01-feature-computation.md) |
| Pathway activity scoring | [`compute.pathway.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md) | [Feature Computation](https://verapancaldilab.github.io/CellTFusion/articles/articles/01-feature-computation.md) |
| Cell group construction | [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md) | [Cell Groups & Latent Factors](https://verapancaldilab.github.io/CellTFusion/articles/articles/02-cell-groups.md) |
| Latent factor extraction | [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md) | [Cell Groups & Latent Factors](https://verapancaldilab.github.io/CellTFusion/articles/articles/02-cell-groups.md) |
| Cell niche derivation | [`compute_cells_niches()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_cells_niches.md) | [Cell Groups & Latent Factors](https://verapancaldilab.github.io/CellTFusion/articles/articles/02-cell-groups.md) |
| Hallmark GSEA per factor | [`compute_factor_gsea()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_factor_gsea.md) | [TME State Characterisation](https://verapancaldilab.github.io/CellTFusion/articles/articles/03-tme-states.md) |
| Meta-program mapping | [`map_factors_to_metaprograms()`](https://verapancaldilab.github.io/CellTFusion/reference/map_factors_to_metaprograms.md) | [TME State Characterisation](https://verapancaldilab.github.io/CellTFusion/articles/articles/03-tme-states.md) |
| TME subtype annotation | [`map_factors_to_TME()`](https://verapancaldilab.github.io/CellTFusion/reference/map_factors_to_TME.md) | [TME State Characterisation](https://verapancaldilab.github.io/CellTFusion/articles/articles/03-tme-states.md) |
| Clinical association testing | [`scores.stat.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.stat.analysis.md) | [Statistical Analysis](https://verapancaldilab.github.io/CellTFusion/articles/articles/04-analysis.md) |
| Survival analysis | [`compute.survival.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.survival.analysis.md) | [Statistical Analysis](https://verapancaldilab.github.io/CellTFusion/articles/articles/04-analysis.md) |
| Test-set projection | [`project_test_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/project_test_factors.md) | [Machine Learning](https://verapancaldilab.github.io/CellTFusion/articles/articles/06-machine-learning.md) |
| ML model training | `pipeML::compute_features.training.ML()` | [Machine Learning](https://verapancaldilab.github.io/CellTFusion/articles/articles/06-machine-learning.md) |
| Multi-cohort (batch) analysis | `batch = TRUE` in [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md) | [Batch/Multi-cohort Analysis](https://verapancaldilab.github.io/CellTFusion/articles/articles/07-batch-analysis.md) |

The full pipeline — every step above run in the right order in a single
call — is available through the
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
wrapper function; see the package
[README](https://VeraPancaldiLab.github.io/CellTFusion/) for usage
examples, including unsupervised, supervised, and multi-cohort modes.

## Quick start

Load the pre-packaged example data and run the full pipeline in one
call:

``` r

raw.counts <- CellTFusion::raw.counts.tuto
traitdata  <- CellTFusion::traitdata.tuto
```

``` r

res <- CellTFusion(
  raw.counts    = raw.counts,
  normalized    = TRUE,
  coldata       = traitdata,
  task          = "unsupervised",
  deconv_methods = c("Quantiseq", "Epidish"),
  cancer_type   = "skcm",
  corr          = 0.7,
  pval          = 0.05,
  file_name     = "Tutorial",
  return        = TRUE
)

# Latent factor scores — use for stat tests and ML
head(res$Latent_spaces$Z)

# TME state annotations
head(res$TME_states)
```

Follow the articles linked in the table above for a step-by-step
explanation of each stage.

## Installation

``` r

remotes::install_github("VeraPancaldiLab/CellTFusion")
```

## References
