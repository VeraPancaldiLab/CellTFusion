# Run Reactome pathway enrichment for a single TF module

Given the hub TFs of a module, extracts their target genes from the
TF-gene universe, keeps the most variable targets, and runs Reactome
over-representation analysis using the full gene expression matrix as
the background universe.

## Usage

``` r
module_enrich(tpm.counts, module_color, hub_genes, tfs_universe)
```

## Arguments

- tpm.counts:

  A numeric matrix of normalized expression values (genes x samples).

- module_color:

  Character. Name of the TF module (color label) to analyze.

- hub_genes:

  A list as returned by
  [`identify_hub_TFs()`](https://verapancaldilab.github.io/CellTFusion/reference/identify_hub_TFs.md),
  where the first element maps module names to vectors of hub TF gene
  symbols.

- tfs_universe:

  A data frame of TF-target interactions with at minimum columns
  `source` (TF) and `target` (gene).

## Value

A `ReactomePA` enrichResult object restricted to pathways with adjusted
p-value \< 0.05, or `NULL` if no enrichment is found.
