# Compute Transcription Factor (TF) activity

Infers transcription factor (TF) activity from a gene expression matrix
using the VIPER algorithm (Alvarez et al., 2016). The function requires
a TF-target gene regulatory network, which can be provided by the user
or obtained from OmnipathR resources such as CollecTRI or Dorothea.
ARACNE-inferred networks are also supported.

## Usage

``` r
# S3 method for class 'TFs.activity'
compute(
  RNA.counts,
  TF.collection = "CollecTRI",
  min_targets_size = 5,
  universe = NULL,
  cancer.type = NULL,
  cores = 3,
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

  - `"CollecTRI"` and `"Dorothea"` use prebuilt collections from
    OmnipathR.

  - `"ARACNE"` allows user input of a custom network file in a 3-column
    format: `regulator`, `target`, and `mutual information`.

- min_targets_size:

  Integer. Minimum number of target genes per regulon required for TF
  activity inference. Default is 5.

- universe:

  Optional. A user-specified data frame of TF-target interactions. If
  not provided, the function will fetch the relevant network based on
  the `TF.collection` argument.

- cancer.type:

  Optional character. Cancer type label used when caching the TF
  collection.

- cores:

  Integer. Number of cores used by VIPER inference. Default is 4.

- scale:

  Logical. If TRUE (default), z-score scales the TF activity matrix
  across samples.

- return:

  Logical; if TRUE, saves matrix in Results/ folder. Default is TRUE.

- file.name:

  Optional character suffix used when writing the TF activity matrix to
  disk.

## Value

A data frame of inferred and scaled TF activity scores, with samples as
rows and TFs as columns.

## References

Alvarez, M. et al. (2016). Functional characterization of somatic
mutations in cancer using network-based inference of protein activity.
*Nature Genetics*, 48(8), 838-847. https://doi.org/10.1038/ng.3593

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
tfs_activity <- compute.TFs.activity(counts.norm.tuto, cores = 1)
#> Warning: One or more parsing issues, call `problems()` on your data frame for details,
#> e.g.:
#>   dat <- vroom(...)
#>   problems(dat)
#> Error in if (.keep) . else select(., -!!evs_col): argument is of length zero
```
