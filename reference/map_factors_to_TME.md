# Annotate NMF factors with Bagaev et al. (2021) MFP subtypes

For a single cancer type, matches TCGA patients present in the NMF
factor score matrix `scores` to the Bagaev et al. (2021) MFP annotation
shipped with the package, then tests each factor across the four MFP
subtypes (IE, IE/F, F, D) with a Kruskal-Wallis test. Factors with p \<
0.05 (fixed threshold, not corrected for the number of factors) are
labelled with the subtype that has the highest median factor score; the
others are labelled `"uncharacterized"`.

## Usage

``` r
map_factors_to_TME(cancer_type, scores, return = TRUE, file_name = NULL)
```

## Arguments

- cancer_type:

  Character. Cancer type abbreviation matching the `TCGA_project` column
  of the internal annotation (case-insensitive, e.g. `"skcm"`).

- scores:

  Numeric matrix. Samples x factors NMF score matrix (row names = TCGA
  barcodes).

- return:

  Logical. If TRUE (default), saves violin/boxplots of factor scores by
  MFP group to `Results/TME_factors_MFP_<cancer_type>_<file_name>.pdf`
  (`cancer_type` in lower case).

- file_name:

  Optional character suffix for saving output plots.

## Value

A data frame with columns `factor`, `best_MFP`, `kw_pval` (3 significant
digits), `median_IE`, `median_IEF`, `median_F`, `median_D`, `n_samples`,
or `NULL` (with a warning) if the cancer type is not annotated or fewer
than 10 patients are matched.

## Details

Samples are matched to patients by the first three parts of their
barcode (`TCGA-XX-XXXX`; barcodes written with "." instead of "-" are
accepted). When a patient has several samples, a tumour sample (sample
type code 01 to 09, e.g. `-01A`) is kept rather than a normal one (e.g.
`-11A`), and a message gives how many samples were left out.

## Examples

``` r
if (FALSE) { # \dontrun{
# scores: latent factors of a TCGA cohort, samples named by their TCGA barcodes
tme <- map_factors_to_TME(cancer_type = "BLCA", scores = res$Latent_spaces$Z)
} # }
```
