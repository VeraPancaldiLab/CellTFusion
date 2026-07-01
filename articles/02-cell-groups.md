# Cell Group Construction

``` r

library(CellTFusion)
#> 
#> 
```

Cell group construction is the core step of the CellTFusion pipeline.
Using previously computed features (TF modules, deconvolution groups,
and pathway activities), we identify clusters of samples that share
similar biological activity patterns as reflected by TF module scores.

This step requires the outputs from the [Feature
Computation](https://verapancaldilab.github.io/CellTFusion/articles/01-feature-computation.md)
article:

- `counts.norm` — log-normalized expression matrix
- `tfs` — TF activity matrix
- `deconv` — cell-type deconvolution proportions
- `network` — WTCNA TF modules
- `dt` — reduced deconvolution groups
- `traitdata` — clinical metadata

## Supervised cell group construction

In the supervised mode, clinical traits guide the identification of cell
groups associated with specific phenotypes (e.g., treatment response). A
trait column and its positive class must be provided:

``` r

cell_groups <- construct_cell_groups(
  counts.norm, tfs, deconv, network, dt,
  traitdata,
  pval     = 0.05,
  trait    = "Best.Confirmed.Overall.Response",
  positive = "CR"
)
```

## Unsupervised cell group construction

Alternatively, cell groups can be identified purely from molecular
features, without using clinical annotations:

``` r

cell_groups <- construct_cell_groups(
  counts.norm, tfs, deconv, network, dt,
  pval = 0.05
)
```

## Output structure

[`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
returns a list with three elements:

| Element | Description                                        |
|---------|----------------------------------------------------|
| `[[1]]` | Data frame of cell group scores (samples × groups) |
| `[[2]]` | Named list of cell-type compositions per group     |
| `[[3]]` | Named list of loading vectors per group            |

Results are also written to `Results/Cell.groups.composition.csv` and
`Results/Cell.groups.scores.csv`.
