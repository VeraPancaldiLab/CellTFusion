# Fisher test using cell groups scores

Fisher test using cell groups scores

## Usage

``` r
cell.groups.fisher.test(cell.groups, coldata, trait, pval = 0.05)
```

## Arguments

- cell.groups:

  A list where the first element is a data frame of cell group scores,
  and the second element contains metadata or labels for these groups.

- coldata:

  A data frame containing the clinical or experimental traits.

- trait:

  Character. Name of the column in `coldata` to test with Fisher's exact
  test.

- pval:

  Numeric. P-value threshold for significance (default 0.05).

## Value

A list containing the significant cell groups after Fisher test.
Additionally, it saves corresponding barplot visualizations in the
"Results/" folder.
