# Compute pathway activity scores

This function computes pathway activity scores from log-normalized gene
expression data with `decoupleR` (Badia-i-Mompel et al., 2022), based on
the PROGENy resource (Schubert et al., 2018). By default the univariate
linear model (`"ulm"`) is used. Optionally, it also performs Gene Set
Variation Analysis (GSVA) using hallmark signatures or any user-provided
gene sets.

## Usage

``` r
compute.pathway.activity(
  RNA.tpm,
  gene_sets = NULL,
  paths = NULL,
  method = "ulm",
  scale = TRUE,
  batch = NULL,
  return = TRUE,
  file.name = NULL
)
```

## Arguments

- RNA.tpm:

  A numeric matrix of log-normalized gene expression values (e.g.
  log2(TPM + 1)) with genes (symbols) as rows and samples as columns. A
  warning is given if it has values above 20, which log2(TPM + 1) cannot
  reach.

- gene_sets:

  A named list of gene sets (character vectors of gene symbols; e.g.,
  hallmark signatures or user-defined sets). If provided, GSVA scores
  will be computed for these sets. Default is `NULL`.

- paths:

  A data frame describing the pathway-gene interactions (columns
  `source`, `target`, `weight`) for use with PROGENy. If `NULL`, the
  human PROGENy resource (top 500 genes) is used and cached in
  `Results/Pathways_collection_PROGENy.csv` (also when `return = FALSE`,
  so that it is downloaded only once). A message gives the number of
  pathways and interactions of the cached file and the date it was
  saved; a cached file that is empty or lacks the needed columns is not
  used (the collection is downloaded again, with a warning). Delete the
  file to download the collection again.

- method:

  Character. The `decoupleR` statistic used to score the PROGENy
  pathways: `"ulm"` (univariate linear model, default), `"mlm"`
  (multivariate linear model), `"wsum"` (weighted sum) or `"consensus"`
  (consensus of the three). The gene sets of `gene_sets` are always
  scored with GSVA.

- scale:

  Logical. If TRUE (default), each pathway is z-scored across the
  samples. The scaled scores of a sample then depend on the other
  samples given in the same call; use `scale = FALSE` to get scores that
  can be compared between separate calls (e.g. a training and a test
  set). With fewer than three samples the scores are not scaled (a
  message is given).

- batch:

  Optional vector with the cohort of each sample, in the column order of
  `RNA.tpm` (if it is named, it is matched to the samples by name
  instead). It is always treated as categorical. If given, the samples
  are split by cohort, the pathway activities (PROGENy and, if
  requested, GSVA) are computed and scaled (see `scale`) separately in
  each cohort, and the results are joined in the sample order of
  `RNA.tpm`; the output has the same form as without `batch`. Pathways
  or gene sets that cannot be scored in every cohort are left out (with
  a message). Default is `NULL` (single cohort); a vector with only one
  cohort is treated the same way.

- return:

  Logical; if TRUE, saves the score matrices in Results/ folder. Default
  is TRUE.

- file.name:

  Optional character suffix used when writing output CSV files
  (`Results/Pathway_matrix_PROGENy_<file.name>.csv` and
  `Results/Pathway_matrix_GSVA_<file.name>.csv`; without it,
  `Results/Pathway_matrix_PROGENy.csv` and
  `Results/Pathway_matrix_GSVA.csv`).

## Value

If `gene_sets` is `NULL`, a data frame of PROGENy pathway activity
scores (samples as rows, pathways as columns), scaled if `scale = TRUE`.
If `gene_sets` is provided, a list with two elements:

- `PROGENy`: A data frame of PROGENy pathway activity scores.

- `GSVA`: A data frame of GSVA scores based on the provided gene sets.

Column names are made syntactically valid with
[`make.names()`](https://rdrr.io/r/base/make.names.html) (e.g.
`JAK.STAT`), in the returned tables and in the saved files. Pathways
with fewer than 5 of their genes in the data, and gene sets that GSVA
cannot score, are left out with a message.

## References

Schubert M, Klinger B, Kluenemann M, Sieber A, Uhlitz F, Sauer S,
Garnett MJ, Bluethgen N, Saez-Rodriguez J. Perturbation-response genes
reveal signaling footprints in cancer gene expression. Nature
Communications. 2018.
[doi:10.1038/s41467-017-02391-6](https://doi.org/10.1038/s41467-017-02391-6)

Badia-i-Mompel, P. et al. (2022). decoupleR: ensemble of computational
methods to infer biological activities from omics data. *Bioinformatics
Advances*, 2(1), vbac016. https://doi.org/10.1093/bioadv/vbac016

## Examples

``` r
if (FALSE) { # \dontrun{
# Compute only PROGENy activities (downloads the PROGENy model from OmniPath)
data("counts.norm.tuto")
pathways <- compute.pathway.activity(counts.norm.tuto)
} # }
```
