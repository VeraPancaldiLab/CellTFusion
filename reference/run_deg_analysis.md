# Run differential expression analysis with edgeR/limma-voom

Filters low-expression genes, applies TMM normalization, runs voom
transformation, fits a linear model, and returns the top differentially
expressed genes via
[`limma::topTable`](https://rdrr.io/pkg/limma/man/toptable.html).

## Usage

``` r
run_deg_analysis(counts, coldata, group_col, ref_level = NULL)
```

## Arguments

- counts:

  A raw count matrix (genes x samples).

- coldata:

  A data frame of sample metadata whose row names match the column names
  of `counts`.

- group_col:

  Character. Name of the column in `coldata` used as the grouping factor
  for differential expression.

- ref_level:

  Character or `NULL`. Reference level for the group factor. If `NULL`,
  the default factor ordering is used.

## Value

A data frame of differentially expressed genes (p.adj \< 0.05) as
returned by
[`limma::topTable`](https://rdrr.io/pkg/limma/man/toptable.html), with
columns `logFC`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`, and `B`.
