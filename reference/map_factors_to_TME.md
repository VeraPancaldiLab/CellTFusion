# Annotate NMF factors with Bagaev et al. (2021) MFP subtypes

For a single cancer type, matches TCGA patients present in the NMF
factor score matrix `Z` to the Bagaev annotation, then computes Spearman
correlations between each factor and each one-hot-encoded MFP subtype
(IE, IE/F, F, D). The best-matching subtype per factor is returned;
factors with a maximum absolute correlation below 0.15 are labelled
`"uncharacterized"`.

## Usage

``` r
map_factors_to_TME(cancer_name, Z, plot = TRUE, file_name = NULL)
```

## Arguments

- cancer_name:

  Character. Cancer type abbreviation matching the `TCGA_project` column
  of the internal annotation (case-insensitive, e.g. `"skcm"`).

- Z:

  Numeric matrix. Samples x factors NMF score matrix (row names = TCGA
  barcodes).

- plot:

  Logical. If TRUE (default), saves a boxplot of factor scores by MFP
  group.

- file_name:

  Optional character. File path prefix for saving output plots.

## Value

A data frame with columns `factor`, `best_MFP`, `kw_pval`, `median_IE`,
`median_IEF`, `median_F`, `median_D`, `n_samples`, or `NULL` if fewer
than 10 patients are matched.
