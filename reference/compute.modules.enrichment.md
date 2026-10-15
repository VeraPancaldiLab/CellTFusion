# Compute TF module enrichment using directed target genes

This function performs enrichment analysis for transcription factor (TF)
modules using known TF-target interactions. For each module, it takes
the hub TFs and retrieves their known targets from the TF-target network
(CollecTRI by default). It then conducts an over-representation analysis
(ORA) against the Reactome pathway database. To reduce redundancy, only
unique pathways per module are retained by filtering out overlaps
between modules.

## Usage

``` r
compute.modules.enrichment(
  counts_norm,
  hub_tfs,
  TF.collection = "CollecTRI",
  universe = NULL,
  cancer_type = NULL,
  top_n = 10,
  file_name = NULL
)
```

## Arguments

- counts_norm:

  A numeric matrix of normalized gene expression values. Rows are gene
  symbols, columns are samples. This matrix is used to restrict the
  universe of genes used in enrichment and validate TF-target gene
  presence.

- hub_tfs:

  A list of hub TFs per module, typically obtained using
  [`identify_hub_TFs()`](https://verapancaldilab.github.io/CellTFusion/reference/identify_hub_TFs.md).
  This should be a list of two elements: one with named hub TFs per
  module, and another with their module membership or additional
  metadata. Only the modules in this list (those with hub TFs) are
  analysed.

- TF.collection:

  Character. The TF-target network from which the targets of the hub TFs
  are taken: `"CollecTRI"` (default), `"Dorothea"` (confidence A and B)
  or `"ARACNE"`, read as in
  [`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md).
  Use the collection the TF activity was computed with.

- universe:

  Optional. A user-specified data frame of TF-target interactions
  (columns `source` and `target`), used instead of `TF.collection`.

- cancer_type:

  Optional character. Cancer type used to locate the ARACNe network
  (`input/ARACNE/<cancer_type>/network/network.txt`), only used when
  `TF.collection = "ARACNE"`. If `NULL`, the network is auto-detected
  when only one exists under `input/ARACNE/`.

- top_n:

  Integer. Number of pathways shown in the dot plot of each module.
  Default is 10. All the pathways are in the returned tables.

- file_name:

  Optional character suffix for the names of the saved plots.

## Value

Invisibly, a list with two elements, each a named list with one table
per module (the result table of
[`ReactomePA::enrichPathway()`](https://rdrr.io/pkg/ReactomePA/man/enrichPathway.html),
or `NULL` for a module without enrichment):

- `enrichment`: the module-specific pathways (adjusted p-value \< 0.05,
  and not significant in any other module).

- `enrichment_all`: all the significant pathways of each module, before
  removing those shared between modules.

For each module with module-specific pathways, a dot plot is also saved
in the `Results/` directory as a PDF file named
`Enrichment_Reactome_Module_<color>_<file_name>.pdf`
(`Enrichment_Reactome_Module_<color>.pdf` without `file_name`). If no
enrichment is found for a module, a message is given and no file is
saved for that module. An error is raised if none of the hub TFs is in
the TF-target network, or if none of their targets is found in the row
names of `counts_norm` (they must be gene symbols).

## Details

The targets of the hub TFs are taken from the TF-target network (for
CollecTRI and Dorothea, cached in
`Results/TF_target_collection_<TF.collection>.csv`), and
[`ReactomePA::enrichPathway()`](https://rdrr.io/pkg/ReactomePA/man/enrichPathway.html)
is used for ORA with the Reactome database. The ORA input is the 20%
most variable targets of the module's hub TFs; a module with fewer than
5 such genes is skipped with a message. The background is all the genes
of `counts_norm`. Any set of TFs has targets enriched in the pathways
those TFs regulate, so the results describe what the hub TFs of each
module regulate, not evidence that the module itself is special.
Pathways enriched in more than one module are removed, so only
module-specific results are kept. TF names are matched in their
syntactic form too (e.g. `NKX2.1`, as returned by
[`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md),
matches `NKX2-1`).

## Examples

``` r
if (FALSE) { # \dontrun{
# Identify the hub TFs of a TF network from compute.WTCNA()
hub_tfs <- identify_hub_TFs(wgcna_network, mm_thresh = 0.8, degree_thresh = 0.9)

# Perform module-level Reactome enrichment
enr <- compute.modules.enrichment(counts.norm, hub_tfs)
head(enr$enrichment[[1]])
} # }
```
