# Build a Hallmarks x factors NES matrix from GSEA results

Combines the per-factor GSEA outputs from
[`compute_factor_gsea()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_factor_gsea.md)
into a single matrix. Hallmarks not significant in a given factor are
filled with 0.

## Usage

``` r
build_nes_matrix(gsea_results)
```

## Arguments

- gsea_results:

  Output list from
  [`compute_factor_gsea()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_factor_gsea.md),
  containing a `GSEA_results` element (named list of `enrichResult`
  objects, one per factor).

## Value

A numeric matrix of NES values with Hallmarks as rows and NMF factors as
columns. Missing Hallmark-factor combinations are set to 0.
