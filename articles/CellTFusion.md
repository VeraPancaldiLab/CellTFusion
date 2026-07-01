# Getting Started with CellTFusion

``` r

library(CellTFusion)
#> 
#> 
```

### Overview

`CellTFusion` integrates immune cell-type deconvolution with
transcription factor (TF)–gene regulatory networks to characterize
immune cell states in the tumor microenvironment from bulk RNA-seq data.

The pipeline takes raw counts and optional clinical metadata as input,
and produces **cell group scores** — biologically coherent clusters of
samples defined by shared TF activity and cell-type composition
patterns.

### Pipeline steps

| Step | Function | Article |
|----|----|----|
| Cell-type deconvolution | [`multideconv::compute.deconvolution()`](https://rdrr.io/pkg/multideconv/man/compute.deconvolution.html) | [Feature Computation](https://verapancaldilab.github.io/CellTFusion/articles/articles/01-feature-computation.md) |
| TF activity inference | [`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md) | [Feature Computation](https://verapancaldilab.github.io/CellTFusion/articles/articles/01-feature-computation.md) |
| TF module construction | [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md) | [Feature Computation](https://verapancaldilab.github.io/CellTFusion/articles/articles/01-feature-computation.md) |
| Pathway activity scoring | [`compute.pathway.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md) | [Feature Computation](https://verapancaldilab.github.io/CellTFusion/articles/articles/01-feature-computation.md) |
| Cell group construction | [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md) | [Cell Group Construction](https://verapancaldilab.github.io/CellTFusion/articles/articles/02-cell-groups.md) |
| Statistical association | [`scores.stat.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.stat.analysis.md) | [Statistical Analysis](https://verapancaldilab.github.io/CellTFusion/articles/articles/03-analysis.md) |
| Full pipeline (one call) | [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md) | [One-step Pipeline](https://verapancaldilab.github.io/CellTFusion/articles/articles/04-one-step-pipeline.md) |
| ML model training | `pipeML::compute_features.training.ML()` | [Machine Learning](https://verapancaldilab.github.io/CellTFusion/articles/articles/05-machine-learning.md) |

### Quick start

Load the pre-packaged example data and run the full pipeline in one
call:

``` r

raw.counts <- CellTFusion::raw.counts.tuto
traitdata  <- CellTFusion::traitdata.tuto
```

``` r

res <- CellTFusion(
  raw.counts     = raw.counts,
  normalized     = TRUE,
  coldata        = traitdata,
  trait          = "Best.Confirmed.Overall.Response",
  trait.positive = "CR",
  deconv_methods = c("Quantiseq", "Epidish"),
  file_name      = "TestRun",
  corr           = 0.7,
  pval           = 0.05,
  return         = FALSE
)
```

For a step-by-step walkthrough of each stage, follow the articles listed
in the table above.

### Installation

``` r

# Install from GitHub
remotes::install_github("VeraPancaldiLab/CellTFusion")
```

## References
