# Fisher's exact test for score-trait association

Fisher's exact test for score-trait association

## Usage

``` r
scores.fisher.test(scores, coldata, trait, pval = 0.05)
```

## Arguments

- scores:

  A list, NMF output from
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md),
  or a score matrix. When a list, the first element must be a samples x
  features score matrix. Continuous scores are binarised at the median
  into High/Low groups.

- coldata:

  A data frame containing the clinical or experimental traits.

- trait:

  Character. Name of the column in `coldata` to test with Fisher's exact
  test.

- pval:

  Numeric. P-value threshold for significance (default 0.05).

## Value

A list containing the significant features after Fisher test.
Additionally, it saves corresponding barplot visualizations in the
"Results/" folder.
