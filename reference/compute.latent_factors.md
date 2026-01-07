# Full NMF pipeline for latent immune states (single cohort)

Full NMF pipeline for latent immune states (single cohort)

## Usage

``` r
compute.latent_factors(X, batch = NULL, seed = 123)
```

## Arguments

- X:

  Numeric matrix of size (samples × cell_groups)

- seed:

  Random seed (default = 123)

- K_range:

  Vector of K to try (default = 2:6)

- remove_low_var:

  Logical, whether to remove low-variance groups (default = TRUE)

- min_var:

  Minimum variance threshold for filtering (default = 1e-5)

- nrun:

  Number of NMF runs per K (default = 10)

## Value

List with: \$best_K = suggested number of latent factors \$W = sample ×
latent states matrix for training \$H = latent states × cell group
contributions \$reconstruction_errors = reconstruction errors for each K
\$consensus = list of consensus matrices for each K
