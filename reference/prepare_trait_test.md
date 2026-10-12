# Prepare the inputs of a score-trait test

Internal helper of the `scores.*` tests. Takes the score matrix (first
element of `scores`), matches the rows of `coldata` to its samples by
name, leaves out samples with a missing trait (with a message) and
checks that the trait is a grouping variable.

## Usage

``` r
prepare_trait_test(scores, coldata, trait, binary = FALSE, test = "test")
```

## Arguments

- scores:

  A list whose first element is a samples x features score matrix, the
  matrix itself, or the output of
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md)
  (its element `Z` is used).

- coldata:

  A data frame (or tibble) of sample annotations.

- trait:

  Character. Name of the column of `coldata` to test.

- binary:

  Logical. If TRUE, the trait must have exactly two groups.

- test:

  Character. Name of the test, used in the messages.

## Value

A list with `scores` (the input as a list), `X` (score matrix of the
samples with a trait value) and `trait` (factor of their trait values).
