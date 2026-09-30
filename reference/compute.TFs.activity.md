# Compute Transcription Factor (TF) activity

Infers transcription factor (TF) activity from a gene expression matrix
with
[`decoupleR::decouple()`](https://saezlab.github.io/decoupleR/reference/decouple.html),
keeping the `consensus` score (ensemble of the decoupleR statistics;
Badia-i-Mompel et al., 2022). The TF-target network can be provided by
the user, obtained from OmnipathR resources (CollecTRI or Dorothea), or
read from an ARACNe-inferred network.

## Usage

``` r
compute.TFs.activity(
  RNA.counts,
  TF.collection = "CollecTRI",
  min_targets_size = 5,
  universe = NULL,
  cancer.type = NULL,
  scale = TRUE,
  return = TRUE,
  file.name = NULL
)
```

## Arguments

- RNA.counts:

  A gene expression matrix with genes as rows and samples as columns.
  The matrix should be normalized (e.g., TPM, log2CPM, etc.).

- TF.collection:

  Character. The source of the TF-target network. Options are
  `"CollecTRI"` (default), `"Dorothea"`, or `"ARACNE"`.

  - `"CollecTRI"` and `"Dorothea"` (confidence A and B) use prebuilt
    collections from OmnipathR. Each collection is cached in its own
    file, `Results/TF_target_collection_<TF.collection>.csv`, and reused
    on later calls.

  - `"ARACNE"` reads a tab-separated network file with `Regulator` and
    `Target` columns from
    `input/ARACNE/<cancer.type>/network/network.txt` (relative to the
    working directory). The mode of regulation of each edge is the sign
    of the Spearman correlation between TF and target expression.

- min_targets_size:

  Integer. Minimum number of target genes per regulon (passed to
  [`decoupleR::decouple()`](https://saezlab.github.io/decoupleR/reference/decouple.html)
  as `minsize`). Default is 5.

- universe:

  Optional. A user-specified data frame of TF-target interactions
  (columns `source`, `target`, `mor`). If not provided, the network is
  fetched based on `TF.collection`. Ignored when
  `TF.collection = "ARACNE"`.

- cancer.type:

  Optional character. TCGA cancer type abbreviation used to locate the
  ARACNe network (only used when `TF.collection = "ARACNE"`). If `NULL`,
  the network is auto-detected when only one `network.txt` exists under
  `input/ARACNE/`.

- scale:

  Logical. If TRUE (default), z-score scales the TF activity matrix
  across samples.

- return:

  Logical; if TRUE, saves matrix in Results/ folder. Default is TRUE.

- file.name:

  Optional character suffix used when writing the TF activity matrix to
  disk.

## Value

A data frame of inferred (and, if `scale = TRUE`, scaled) TF activity
scores, with samples as rows and TFs as columns. Column names are made
syntactically valid with
[`make.names()`](https://rdrr.io/r/base/make.names.html).

## References

Badia-i-Mompel, P. et al. (2022). decoupleR: ensemble of computational
methods to infer biological activities from omics data. *Bioinformatics
Advances*, 2(1), vbac016. https://doi.org/10.1093/bioadv/vbac016

Tuerei, D., Korcsmaros, T., & Saez-Rodriguez, J. (2016). OmniPath:
guidelines and gateway for literature-curated signaling pathway
resources. *Nature Methods*, 13(12), 966-967.
https://doi.org/10.1038/nmeth.4077

Garcia-Alonso, L. et al. (2019). Benchmark and integration of resources
for the estimation of human transcription factor activities. *Genome
Research*. https://doi.org/10.1101/gr.240663.118

Lachmann, A. et al. (2016). ARACNe-AP: gene network reverse engineering
through adaptive partitioning inference of mutual information.
*Bioinformatics*, 32(14), 2233-2235.
https://doi.org/10.1093/bioinformatics/btw216

Margolin, A.A. et al. (2006). ARACNE: an algorithm for the
reconstruction of gene regulatory networks in a mammalian cellular
context. *BMC Bioinformatics*, 7(Suppl 1), S7.
https://doi.org/10.1186/1471-2105-7-S1-S7

## Examples

``` r
data("counts.norm.tuto")
tfs_activity <- compute.TFs.activity(counts.norm.tuto)
#> Warning: One or more parsing issues, call `problems()` on your data frame for details,
#> e.g.:
#>   dat <- vroom(...)
#>   problems(dat)
#> Error in if (.keep) . else select(., -!!evs_col): argument is of length zero
```
