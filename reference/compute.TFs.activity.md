# Compute Transcription Factor (TF) activity

Infers transcription factor (TF) activity from a gene expression matrix
with
[`decoupleR::decouple()`](https://saezlab.github.io/decoupleR/reference/decouple.html)
(Badia-i-Mompel et al., 2022). By default the univariate linear model
(`"ulm"`) is used. The TF-target network can be provided by the user,
obtained from prebuilt collections (CollecTRI from OmnipathR or Dorothea
from the `dorothea` package), or read from an ARACNe-inferred network.

## Usage

``` r
compute.TFs.activity(
  RNA.counts,
  TF.collection = "CollecTRI",
  min_targets_size = 5,
  universe = NULL,
  cancer.type = NULL,
  method = "ulm",
  scale = TRUE,
  return = TRUE,
  file.name = NULL
)
```

## Arguments

- RNA.counts:

  A gene expression matrix with genes (symbols) as rows and samples as
  columns, with unique sample names. The matrix should be log-normalized
  (e.g., log2(TPM + 1), log2CPM, etc.). A warning is given if it has
  values above 20, which log2(TPM + 1) cannot reach.

- TF.collection:

  Character. The source of the TF-target network. Options are
  `"CollecTRI"` (default), `"Dorothea"`, or `"ARACNE"`.

  - `"CollecTRI"` uses the prebuilt collection from OmnipathR and
    `"Dorothea"` (confidence A and B) the one from the `dorothea`
    package. Each collection is cached in its own file,
    `Results/TF_target_collection_<TF.collection>.csv`, and reused on
    later calls: a message gives the number of TFs and interactions of
    the cached file and the date it was saved. A cached file that is
    empty or lacks the columns `source`, `target` and `mor` is not used
    (the collection is obtained again, with a warning). Delete the file
    to obtain the collection again.

  - `"ARACNE"` reads a tab-separated network file with `Regulator` and
    `Target` columns from
    `input/ARACNE/<cancer.type>/network/network.txt` (relative to the
    working directory). The mode of regulation of each edge is read from
    a `mor` column of the file (1 activation, -1 repression), which
    should be the sign of the Spearman correlation between TF and target
    computed on the cohort used to build the network. If the file has no
    `mor` column, a warning is given and the signs are computed as the
    sign of the Spearman correlation between TF and target in the
    samples of `RNA.counts`; they then depend on the samples given (e.g.
    they can differ between cohorts, or between a training and a test
    set).

- min_targets_size:

  Integer. Minimum number of target genes per regulon (passed to
  [`decoupleR::decouple()`](https://saezlab.github.io/decoupleR/reference/decouple.html)
  as `minsize`). Default is 5.

- universe:

  Optional. A user-specified data frame of TF-target interactions
  (columns `source`, `target`, `mor`). The three columns are required:
  `mor` is the mode of regulation (1 activation, -1 repression; use 1
  for every edge if the network has no sign). If not provided, the
  network is fetched based on `TF.collection`. Ignored when
  `TF.collection = "ARACNE"`.

- cancer.type:

  Optional character. TCGA cancer type abbreviation used to locate the
  ARACNe network (only used when `TF.collection = "ARACNE"`). If `NULL`,
  the network is auto-detected when only one `network.txt` exists under
  `input/ARACNE/`.

- method:

  Character. The `decoupleR` statistic used to score the TFs: `"ulm"`
  (univariate linear model, default), `"mlm"` (multivariate linear
  model), `"wsum"` (weighted sum) or `"consensus"` (consensus of the
  three).

- scale:

  Logical. If TRUE (default), z-score scales the TF activity matrix
  across samples. The scaled scores of a sample then depend on the other
  samples given in the same call; use `scale = FALSE` to get scores that
  can be compared between separate calls (e.g. a training and a test
  set). With fewer than three samples the scores are not scaled (a
  message is given).

- return:

  Logical; if TRUE, saves matrix in Results/ folder. Default is TRUE.

- file.name:

  Optional character suffix used when writing the TF activity matrix to
  disk (`Results/TF_matrix_<file.name>.csv`; without it,
  `Results/TF_matrix.csv`).

## Value

A data frame of inferred (and, if `scale = TRUE`, scaled) TF activity
scores, with samples as rows and TFs as columns. Column names are made
syntactically valid with
[`make.names()`](https://rdrr.io/r/base/make.names.html) (e.g.
`NKX2.1`), in the returned table and in the saved file.

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
if (FALSE) { # \dontrun{
# Downloads the CollecTRI network from OmniPath
data("counts.norm.tuto")
tfs_activity <- compute.TFs.activity(counts.norm.tuto)
} # }
```
