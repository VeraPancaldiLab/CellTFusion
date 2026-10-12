# Run Reactome pathway enrichment for a single TF module

Given the hub TFs of a module, extracts their target genes from the
TF-gene universe, keeps the 20% most variable targets, and runs Reactome
over-representation analysis
([`ReactomePA::enrichPathway()`](https://rdrr.io/pkg/ReactomePA/man/enrichPathway.html),
human) using all the genes of the expression matrix as the background
universe. Gene symbols are converted to Entrez IDs with `org.Hs.eg.db`;
symbols without an Entrez ID are left out (a message gives how many).
Pathways with a Benjamini-Hochberg adjusted p-value \< 0.05 are kept. TF
names are matched after
[`make.names()`](https://rdrr.io/r/base/make.names.html) on both sides,
so the syntactic names returned by
[`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md)
(e.g. `NKX2.1`) match the names of the network (`NKX2-1`).

## Usage

``` r
module_enrich(tpm.counts, module_color, hub_genes, tfs_universe, min_genes = 5)
```

## Arguments

- tpm.counts:

  A numeric matrix of normalized expression values (genes x samples),
  with gene symbols as row names.

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

- min_genes:

  Integer. Minimum number of genes to test (20% most variable targets);
  with fewer, no enrichment is computed and a message is given. Default
  is 5.

## Value

A `ReactomePA` enrichResult object restricted to pathways with adjusted
p-value \< 0.05, or `NULL` if no enrichment is found, if none of the
targets is in the row names of `tpm.counts`, or if fewer than
`min_genes` genes would be tested.

## Details

The result describes the genes regulated by the hub TFs of the module:
any set of TFs has hundreds of targets that carry the pathways those TFs
regulate, so random sets of TFs also give enriched pathways. It should
be read as what the hub TFs of the module regulate, not as evidence that
the module itself is special.
