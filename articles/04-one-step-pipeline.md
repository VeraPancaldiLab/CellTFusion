# One-step Pipeline

``` r

library(CellTFusion)
#> 
#> 
```

For convenience,
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
is an all-in-one wrapper that automates the entire pipeline — computing
features, running intermediate analyses, and returning cell group scores
— in a single call.

The individual steps (feature computation, cell group construction,
statistical analysis) are described in their respective articles:

- [Feature
  Computation](https://verapancaldilab.github.io/CellTFusion/articles/01-feature-computation.md)
- [Cell Group
  Construction](https://verapancaldilab.github.io/CellTFusion/articles/02-cell-groups.md)
- [Statistical
  Analysis](https://verapancaldilab.github.io/CellTFusion/articles/03-analysis.md)

## Running the full pipeline

``` r

library(CellTFusion)

raw.counts <- CellTFusion::raw.counts.tuto
traitdata  <- CellTFusion::traitdata.tuto

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
  high_corr_groups = 0.85,
  return         = FALSE
)
```

## Key parameters

| Parameter | Description |
|----|----|
| `raw.counts` | Raw or normalized count matrix (genes × samples) |
| `normalized` | Set `TRUE` if counts are already log-normalized |
| `coldata` | Data frame of clinical metadata (optional) |
| `trait` | Column in `coldata` to use for supervised analysis |
| `trait.positive` | Positive class label for the trait |
| `deconv_methods` | Deconvolution algorithms to use |
| `corr` | Correlation threshold for grouping cell types |
| `pval` | P-value threshold for significance filters |
| `high_corr_groups` | Threshold for merging highly correlated cell groups |
| `return` | If `TRUE`, returns the result object instead of writing to disk |

## Output

[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
returns a named list containing:

| Element                    | Description                          |
|----------------------------|--------------------------------------|
| `$Cell_groups`             | Cell group scores and composition    |
| `$Processed_deconvolution` | Reduced deconvolution feature matrix |
| `$TF_modules`              | WTCNA network output                 |
| `$Pathways`                | PROGENy pathway activity scores      |
