# Project cell group scores onto trained NMF latent factors

Keeps the positive part of the cell group scores of new samples (as in
training) and projects them onto the latent space learned during
training using the Moore-Penrose pseudoinverse of the scaled basis
matrix W. Negative values of the projection are set to 0 (on the
tutorial data the result is the same as an exact non-negative
projection). Cell groups missing in `scores` are tolerated with a
warning naming them: the projection uses the other cell groups (on the
tutorial data, with 20 of 144 cell groups missing, the factors still
correlate 0.999 with the full projection).

## Usage

``` r
project_factors(latent_factors, scores)
```

## Arguments

- latent_factors:

  A list returned by
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md),
  containing at minimum `W` (features x rank) and the RcppML model
  object with the `nmf_model$d` scaling vector.

- scores:

  A samples x cell groups matrix of signed CCA composite scores for the
  test cohort (for new samples, computed with
  [`compute.test.set()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.test.set.md)).
  Column names must match the cell groups used in training.

## Value

A non-negative matrix of size samples x rank containing the projected
NMF factor scores for the test cohort.

## Details

The RcppML NMF decomposition is A = W %*% diag(d) %*% H, not A = W %\*%
H. The diagonal scaling vector d must be absorbed into W before
computing the pseudoinverse, otherwise the projected scores are on the
wrong scale.

## Examples

``` r
data("network.tuto")
data("deconv_subgroups.tuto")
cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Dendrogram plots saved to: Results/Dendrogram_color_clusters.pdf
#> 44 of 897 features were not assigned to any cluster by the tree cut and are not part of any cell group (pink: 2, yellow: 8, red: 2, blue: 6, brown: 4, turquoise: 12, black: 4, green: 6)
latent <- compute.latent_factors(cell_groups$Cell_groups, return = FALSE)
#> Selecting the rank with bootstrap NMF (100 runs per rank, ranks 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)...
#> Selected rank: 5 (largest rank with all factors robust and none redundant)
#> Running RcppML NMF with rank = 5...
projected <- project_factors(latent, cell_groups$Cell_groups[1:10, ])
head(projected)
#>                      Factor1     Factor2     Factor3      Factor4      Factor5
#> SAM7f0d9cc7f001 0.0000000000 0.000000000 0.016288028 0.0071367700 1.297582e-02
#> SAM4305ab968b90 0.0001494108 0.003315123 0.001295361 0.0000000000 0.000000e+00
#> SAMcf018fee2acd 0.0000000000 0.004341070 0.015464817 0.0058149063 6.239588e-03
#> SAMcc4675f394a1 0.0019439433 0.000000000 0.016317355 0.0000000000 5.190626e-03
#> SAM49f9b2e57aa5 0.0008856227 0.003758303 0.002254951 0.0000000000 5.060087e-03
#> SAM2e7aa8fa0ab3 0.0006108836 0.005278085 0.002614708 0.0003165547 1.984736e-05
```
