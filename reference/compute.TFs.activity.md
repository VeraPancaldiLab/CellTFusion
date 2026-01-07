# Compute Transcription Factor (TF) activity

Infers transcription factor (TF) activity from a gene expression matrix
using the VIPER algorithm (Alvarez et al., 2016). The function requires
a TF–target gene regulatory network, which can be provided by the user
or obtained from OmnipathR resources such as CollecTRI or Dorothea.
ARACNE-inferred networks are also supported.

## Usage

``` r
compute.TFs.activity(
  RNA.counts,
  TF.collection = "CollecTRI",
  min_targets_size = 5,
  tfs.pruned = FALSE,
  universe = NULL,
  return = TRUE
)
```

## Arguments

- RNA.counts:

  A gene expression matrix with genes as rows and samples as columns.
  The matrix should be normalized (e.g., TPM, log2CPM, etc.).

- TF.collection:

  Character. The source of the TF–target network. Options are
  `"CollecTRI"` (default), `"Dorothea"`, or `"ARACNE"`.

  - `"CollecTRI"` and `"Dorothea"` use prebuilt collections from
    OmnipathR.

  - `"ARACNE"` allows user input of a custom network file in a 3-column
    format: `regulator`, `target`, and `mutual information`.

- min_targets_size:

  Integer. Minimum number of target genes per regulon required for TF
  activity inference. Default is 5.

- tfs.pruned:

  Logical. Whether to prune TF regulons to limit the number of target
  genes, which helps reduce bias introduced by TFs with large regulons.
  If `TRUE`, the user will be prompted to input a maximum size for
  regulons. Default is `FALSE`.

- universe:

  Optional. A user-specified data frame of TF–target interactions. If
  not provided, the function will fetch the relevant network based on
  the `TF.collection` argument.

- return:

  Logical; if TRUE, saves matrix in Results/ folder. Default is TRUE.

## Value

A data frame of inferred and scaled TF activity scores, with samples as
rows and TFs as columns.

## References

Alvarez, M. et al. (2016). Functional characterization of somatic
mutations in cancer using network-based inference of protein activity.
*Nature Genetics*, 48(8), 838–847. https://doi.org/10.1038/ng.3593

Türei, D., Korcsmáros, T., & Saez-Rodriguez, J. (2016). OmniPath:
guidelines and gateway for literature-curated signaling pathway
resources. *Nature Methods*, 13(12), 966–967.
https://doi.org/10.1038/nmeth.4077

Garcia-Alonso, L. et al. (2019). Benchmark and integration of resources
for the estimation of human transcription factor activities. *Genome
Research*. https://doi.org/10.1101/gr.240663.118

Lachmann, A. et al. (2016). ARACNe-AP: gene network reverse engineering
through adaptive partitioning inference of mutual information.
*Bioinformatics*, 32(14), 2233–2235.
https://doi.org/10.1093/bioinformatics/btw216

Margolin, A.A. et al. (2006). ARACNE: an algorithm for the
reconstruction of gene regulatory networks in a mammalian cellular
context. *BMC Bioinformatics*, 7(Suppl 1), S7.
https://doi.org/10.1186/1471-2105-7-S1-S7

## Examples

``` r
data("counts.norm.tuto")
tfs_activity <- compute.TFs.activity(counts.norm.tuto)
#> Error in if (.keep) . else select(., -!!evs_col): argument is of length zero
```
