# Module enrichment

Module enrichment

## Usage

``` r
module_enrich(tpm.counts, module_color, hub_genes, tfs_universe)
```

## Arguments

- tpm.counts:

  A matrix with normalized counts (genes as rows and samples as columns)

- module_color:

  A character vector with TF module colors.

- hub_genes:

  List of hub TFs per module.

- tfs_universe:

  A matrix with TF-gene interactions

## Value

Reactome results
