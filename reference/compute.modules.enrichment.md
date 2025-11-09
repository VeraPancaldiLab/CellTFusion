# Compute TF module enrichment using directed target genes

This function performs enrichment analysis for transcription factor (TF)
modules using known TF-target interactions. For each module, it
identifies hub TFs and retrieves their known targets from the CollecTRI
database. It then conducts an over-representation analysis (ORA) against
the Reactome pathway database. To reduce redundancy, only unique
pathways per module are retained by filtering out overlaps between
modules.

## Usage

``` r
compute.modules.enrichment(RNA.tpm, hub_tfs)
```

## Arguments

- RNA.tpm:

  A numeric matrix of normalized gene expression values. Rows are gene
  symbols, columns are samples. This matrix is used to restrict the
  universe of genes used in enrichment and validate TF-target gene
  presence.

- hub_tfs:

  A list of hub TFs per module, typically obtained using
  [`identify_hub_TFs()`](https://verapancaldilab.github.io/CellTFusion/reference/identify_hub_TFs.md).
  This should be a list of two elements: one with named hub TFs per
  module, and another with their module membership or additional
  metadata.

## Value

No object is returned. For each module with significant enrichment
(p-value \< 0.05), a dot plot is saved in the `Results/` directory as a
PDF file named `Module <color>.pdf`. If no enrichment is found for a
module, a message is printed and no file is saved for that module.

## Details

The function uses the
[`decoupleR::get_collectri()`](https://saezlab.github.io/decoupleR/reference/get_collectri.html)
function to obtain TF-target relationships, and
`clusterProfiler::enrichPathway()` for ORA using the Reactome database.
Pathways shared between multiple modules are filtered using a Venn
diagram-based comparison to retain only module-specific results.

## Examples

``` r
if (FALSE) { # \dontrun{
# Identify hub TFs from a TF activity matrix and WGCNA network
hub_tfs <- identify_hub_TFs(t(tfs_activity), wgcna_network, MM_thresh = 0.8, degree_thresh = 0.9)

# Perform module-level Reactome enrichment
compute.modules.enrichment(counts.norm, hub_tfs)
} # }
```
