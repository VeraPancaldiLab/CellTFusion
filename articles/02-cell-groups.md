# Cell Group Construction and Latent Factors

``` r

library(CellTFusion)
#> 
#> 
```

This article covers the core analytical steps of CellTFusion:
identifying cell groups from TF modules and deconvolution data,
extracting latent factors, and deriving cell niches. These steps require
outputs from the [Feature
Computation](https://verapancaldilab.github.io/CellTFusion/articles/01-feature-computation.md)
article:

- `network` — WTCNA TF modules from
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)
- `dt` — reduced deconvolution groups from
  [`multideconv::compute.deconvolution.analysis()`](https://rdrr.io/pkg/multideconv/man/compute.deconvolution.analysis.html)
- `latent_spaces` — NMF latent factors from
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md)

## Construct cell groups

Cell groups are identified by correlating TF module eigengenes with
deconvolution subgroups and cutting the resulting dendrogram. The
function takes only two primary inputs — `network` and `dt` — both
produced in the feature computation step.

``` r

cell_groups <- construct_cell_groups(
  network           = network,
  dt                = dt,
  batch             = NULL,          # pass batch vector if multi-cohort
  pval              = 0.05,
  clustering.method = "ward.D2",
  n_perm            = 999,
  return_dendrogram = FALSE
)
```

For multi-cohort data, pass the batch vector used in
[`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md):

``` r

batch_vec   <- traitdata[, "Cohort"]
cell_groups <- construct_cell_groups(network, dt, batch = batch_vec, pval = 0.05)
```

### Output structure

[`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
returns a three-element list:

| Element | Description                                            |
|---------|--------------------------------------------------------|
| `[[1]]` | Data frame of composite scores (samples × cell groups) |
| `[[2]]` | Named list of cell-type compositions per group         |
| `[[3]]` | Named list of TF loading vectors per group             |

Results are also written to `Results/Cell.groups.composition.csv` and
`Results/Cell.groups.scores.csv`.

## Extract latent factors

Cell group scores are decomposed into NMF latent factors to obtain a
compact, non-negative representation of the TME landscape. The `rank`
argument controls the number of factors; if `NULL`, it is selected
automatically.

``` r

latent_spaces <- compute.latent_factors(
  X         = cell_groups[[1]],
  rank      = NULL,
  seed      = 123,
  file_name = "Tutorial",
  return    = TRUE
)
```

The result is a list whose `$Z` element is a samples × factors matrix —
the primary input for downstream statistical tests and machine learning.

``` r

# Access the factor scores matrix (samples x factors)
head(latent_spaces$Z)
```

## Extract cell niches

Cell niches are derived from latent factors by enriching each factor for
specific cell-type signals from the deconvolution output. This step
characterises each latent factor with a biological cell-type context.

``` r

cells_niches <- compute_cells_niches(
  latent_spaces   = latent_spaces,
  dt              = dt,
  cell.groups     = cell_groups,
  enrich_thresh   = 1.5,
  quantile_cutoff = 0.7,
  return          = TRUE,
  file_name       = "Tutorial"
)
```

## Next steps

With `latent_spaces` in hand, you can:

- **Characterise TME states** via GSEA and meta-program mapping → [TME
  States](https://verapancaldilab.github.io/CellTFusion/articles/03-tme-states.md)
- **Test clinical associations** of latent factors → [Statistical
  Analysis](https://verapancaldilab.github.io/CellTFusion/articles/04-analysis.md)
- **Train ML models** on `latent_spaces$Z` → [Machine
  Learning](https://verapancaldilab.github.io/CellTFusion/articles/06-machine-learning.md)
