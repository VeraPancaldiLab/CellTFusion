# Compute latent factors from cell-group features

Compute latent factors from cell-group features

## Usage

``` r
# S3 method for class 'latent_factors'
compute(X, batch = NULL, seed = 123)
```

## Arguments

- X:

  Numeric matrix of size samples x cell groups.

- batch:

  Optional vector indicating batch assignment for samples.

- seed:

  Random seed used to initialize model fitting.

## Value

A list with latent representation and loadings:

- `Z`: Sample-level latent factors.

- `W`: Feature weights per latent factor.

- `expanded_Z`: Expanded latent features used downstream.
