# Annotate NMF factors with Bagaev et al. (2021) MFP subtypes

For a single cancer type, matches TCGA patients present in the NMF
factor score matrix `Z` to the Bagaev et al. (2021) MFP annotation
shipped with the package, then tests each factor across the four MFP
subtypes (IE, IE/F, F, D) with a Kruskal-Wallis test. Factors with p \<
0.05 are labelled with the subtype that has the highest median factor
score; the others are labelled `"uncharacterized"`.

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

  Logical. If TRUE (default), saves violin/boxplots of factor scores by
  MFP group to `Results/TME_factors_MFP_<cancer_name>_<file_name>.pdf`.

- file_name:

  Optional character suffix for saving output plots.

## Value

A data frame with columns `factor`, `best_MFP`, `kw_pval`, `median_IE`,
`median_IEF`, `median_F`, `median_D`, `n_samples`, or `NULL` (with a
warning) if the cancer type is not annotated or fewer than 10 patients
are matched.
