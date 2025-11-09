# Compute projected cell group scores on an independent cohort

This function computes the projection of previously defined cell groups
onto an independent test cohort by applying the same composite scoring
strategy used in the training data.

## Usage

``` r
compute_cell_groups_signatures(
  deconv_res,
  cell_groups,
  features,
  deconvolution_test,
  tfs.module.network
)
```

## Arguments

- deconv_res:

  A list resulting from
  [`multideconv::compute.deconvolution.analysis()`](https://rdrr.io/pkg/multideconv/man/compute.deconvolution.analysis.html)
  on the training cohort. This list should include a component named
  "Deconvolution subgroups composition".

- cell_groups:

  A list output from
  [`cell.groups.computation()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.computation.md),
  containing cell group scores, feature compositions, and loadings.

- features:

  A character vector of selected cell group names to be projected on the
  test data.

- deconvolution_test:

  A matrix or data frame with deconvolution features from an independent
  test cohort. Raw (unprocessed) deconvolution output is expected.

- tfs.module.network:

  A list from
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)
  containing the TF module network information, used to calculate
  composite scores for each group.

## Value

A data frame with samples as rows and projected cell group scores as
columns. Each column corresponds to one of the selected cell groups in
`features`.

## Examples

``` r
if (FALSE) { # \dontrun{
cell.groups.projected <- compute_cell_groups_signatures(
    deconv_res = deconv_training,
    cell_groups = cell.groups,
    features = selected_features,
    deconvolution_test = deconv_test,
    tfs.module.network = network
)
} # }
```
