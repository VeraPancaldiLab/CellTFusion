# Build the result of a score-trait test

Internal helper of the `scores.*` tests: adjusts the p-values, selects
the significant features and returns the significant part of `scores`
with a table of the tests of all features.

## Usage

``` r
finish_trait_test(input, statistic, p, test, statistic_name, pval, p_adjust)
```

## Arguments

- input:

  The output of
  [`prepare_trait_test()`](https://verapancaldilab.github.io/CellTFusion/reference/prepare_trait_test.md).

- statistic:

  Numeric vector with the test statistic of each feature.

- p:

  Numeric vector with the p-value of each feature.

- test:

  Character. Name of the test.

- statistic_name:

  Character. Name of the statistic (e.g. "t").

- pval:

  Numeric. Significance threshold for the adjusted p-values.

- p_adjust:

  Character. Method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html).

## Value

A list with `result` (see the `scores.*` tests), `sig` (indices of the
significant features) and `adjusted_text` (adjusted p-value of each
feature for the plot subtitles, empty without adjustment).
