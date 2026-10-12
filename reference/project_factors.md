# Project cell group scores onto trained NMF latent factors

Applies the positive-negative split to new samples and projects them
onto the latent space learned during training using the Moore-Penrose
pseudoinverse of the scaled basis matrix W. Negative values of the
projection are set to 0 (on the tutorial data the result is the same as
an exact non-negative projection). Cell groups missing in `scores_test`
are tolerated with a warning naming them: the projection uses the other
cell groups (on the tutorial data, with 20 of 144 cell groups missing,
the factors still correlate 0.999 with the full projection).

## Usage

``` r
project_factors(latent_spaces, scores_test)
```

## Arguments

- latent_spaces:

  A list returned by
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md),
  containing at minimum `W` (features x rank) and the RcppML model
  object with the `nmf_model$d` scaling vector.

- scores_test:

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
#> 20 of 764 features were not assigned to any cluster by the tree cut and are not part of any cell group (blue: 6, pink: 2, red: 2, yellow: 8, black: 2)
latent <- compute.latent_factors(cell_groups$Cell_groups, return = FALSE)
#> Estimating optimal rank via elbow on reconstruction MSE (ranks 2:8)...
#> MSE by rank: k2=0.00127  k3=0.000964  k4=0.000838  k5=0.000756  k6=0.00071  k7=0.000675  k8=0.000632
#> Selected rank (elbow): 3
#> Running RcppML NMF with rank = 3...
projected <- project_factors(latent, cell_groups$Cell_groups[1:10, ])
head(projected)
#>                     Factor1     Factor2      Factor3
#> SAM7f0d9cc7f001 0.000000000 0.007094249 9.496681e-03
#> SAM4305ab968b90 0.019944131 0.000000000 1.019721e-05
#> SAMcf018fee2acd 0.003743711 0.002179128 5.920653e-03
#> SAMcc4675f394a1 0.007208812 0.007875054 0.000000e+00
#> SAM49f9b2e57aa5 0.008377438 0.000000000 2.906892e-03
#> SAM2e7aa8fa0ab3 0.003268215 0.003791681 0.000000e+00
```
