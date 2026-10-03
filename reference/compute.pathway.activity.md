# Compute pathway activity scores

This function computes pathway activity scores from normalized gene
expression data using a multivariate linear model (MLM) based on the
PROGENy resource (Schubert et al., 2018). Optionally, it also performs
Gene Set Variation Analysis (GSVA) using hallmark signatures or any
user-provided gene sets.

## Usage

``` r
compute.pathway.activity(
  RNA.tpm,
  gene_sets = NULL,
  paths = NULL,
  return = TRUE,
  file.name = NULL
)
```

## Arguments

- RNA.tpm:

  A numeric matrix of normalized gene expression values with genes as
  rows and samples as columns.

- gene_sets:

  A named list of gene sets (e.g., hallmark signatures or user-defined
  sets). If provided, GSVA scores will be computed for these sets.
  Default is `NULL`.

- paths:

  A data frame describing the pathway-gene interactions (columns
  `source`, `target`, `weight`) for use with PROGENy. If `NULL`, the
  human PROGENy resource (top 500 genes) is used and cached in
  `Results/Pathways_collection_PROGENy.csv`.

- return:

  Logical; if TRUE, saves matrices in Results/ folder. Default is TRUE.

- file.name:

  Optional character suffix used when writing output CSV files.

## Value

If `gene_sets` is `NULL`, a scaled data frame of PROGENy pathway
activity scores (samples as rows, pathways as columns). If `gene_sets`
is provided, a list with two elements:

- `PROGENy`: A scaled data frame of PROGENy pathway activity scores.

- `GSVA`: A scaled data frame of GSVA scores based on the provided gene
  sets.

Column names are made syntactically valid with
[`make.names()`](https://rdrr.io/r/base/make.names.html).

## References

Schubert M, Klinger B, Kluenemann M, Sieber A, Uhlitz F, Sauer S,
Garnett MJ, Bluethgen N, Saez-Rodriguez J. Perturbation-response genes
reveal signaling footprints in cancer gene expression. Nature
Communications. 2018.
[doi:10.1038/s41467-017-02391-6](https://doi.org/10.1038/s41467-017-02391-6)

## Examples

``` r
if (FALSE) { # \dontrun{
# Compute only PROGENy activities (downloads the PROGENy model from OmniPath)
data("counts.norm.tuto")
pathways <- compute.pathway.activity(counts.norm.tuto)
} # }
```
