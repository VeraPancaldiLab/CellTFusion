# Computes TF-modules pathway activities scores

This function computes pathway activity scores from normalized gene
expression data using a multivariate linear model (MLM) based on the
PROGENy resource (Schubert et al., 2018). Optionally, it also performs
Gene Set Variation Analysis (GSVA) using hallmark signatures or any
user-provided gene sets.

## Usage

``` r
# S3 method for class 'pathway.activity'
compute(
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

  A list of gene sets (e.g., hallmark signatures or user-defined sets).
  If provided, GSVA scores will be computed for these sets. Default is
  `NULL`.

- paths:

  A data frame describing the pathway-gene interactions for use with
  PROGENy. If `NULL`, the human PROGENy resource (top 500 genes) will be
  used by default.

- return:

  Logical; if TRUE, saves matrices in Results/ folder. Default is TRUE.

- file.name:

  Optional character suffix used when writing output CSV files.

## Value

If `gene_sets` is `NULL`, a scaled matrix of PROGENy pathway activity
scores (samples as rows, pathways as columns). If `gene_sets` is
provided, a list with two elements:

- `sample_acts_progeny`: A scaled matrix of PROGENy pathway activity
  scores.

- `sample_acts_gsva`: A scaled matrix of GSVA scores based on the
  provided gene sets.

## References

Schubert M, Klinger B, Kluenemann M, Sieber A, Uhlitz F, Sauer S,
Garnett MJ, Bluethgen N, Saez-Rodriguez J. Perturbation-response genes
reveal signaling footprints in cancer gene expression. Nature
Communications. 2018.
[doi:10.1038/s41467-017-02391-6](https://doi.org/10.1038/s41467-017-02391-6)

## Examples

``` r
# Compute only PROGENy activities
data("counts.norm.tuto")
pathways <- compute.pathway.activity(counts.norm.tuto)
#> Warning: 'OmnipathR::get_annotation_resources' is deprecated.
#> Use 'annotation_resources' instead.
#> See help("Deprecated")
#> Warning: 'OmnipathR::import_omnipath_annotations' is deprecated.
#> Use 'annotations' instead.
#> See help("Deprecated")
#> Warning: cannot open file 'Results/Pathways_collection_PROGENy.csv': No such file or directory
#> Error in file(file, ifelse(append, "a", "w")): cannot open the connection
```
