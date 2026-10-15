# Project cell group scores onto trained NMF latent factors

Applies the positive-negative split to new samples and projects them
onto the latent space learned during training using the Moore-Penrose
pseudoinverse of the scaled basis matrix W. Negative values of the
projection are set to 0 (on the tutorial data the result is the same as
an exact non-negative projection). Cell groups missing in `scores` are
tolerated with a warning naming them: the projection uses the other cell
groups (on the tutorial data, with 20 of 144 cell groups missing, the
factors still correlate 0.999 with the full projection).

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
  test cohort. Column names must match those used in training (before
  the \_pos/\_neg suffix was added).

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
#> Estimating optimal rank via elbow on reconstruction MSE (ranks 2:8)...
#> MSE by rank: k2=0.00129  k3=0.00114  k4=0.00105  k5=0.000978  k6=0.000933  k7=0.000893  k8=0.000857
#> Selected rank (elbow): 3
#> Running RcppML NMF with rank = 3...
projected <- project_factors(latent, cell_groups$Cell_groups[1:10, ])
head(projected)
#>                      Factor1      Factor2     Factor3
#> SAM7f0d9cc7f001 0.0001801689 0.0083329667 0.005554023
#> SAM4305ab968b90 0.0197346289 0.0000000000 0.000000000
#> SAMcf018fee2acd 0.0018939304 0.0063142198 0.005344026
#> SAMcc4675f394a1 0.0048871177 0.0087565089 0.000000000
#> SAM49f9b2e57aa5 0.0086966300 0.0006169756 0.001544994
#> SAM2e7aa8fa0ab3 0.0046110226 0.0017782760 0.001931307
```
