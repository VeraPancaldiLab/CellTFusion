# Map study factors to TCGA meta-programs

For each study NMF factor, scores it against each TCGA meta-program by
computing the mean NES of the meta-program's Hallmarks in that factor.
The meta-program with the highest mean NES is the best match. If Bagaev
MFP annotations are available (either embedded in the loaded
meta-program object or in
`inst/extdata/bagaev_factor_annotations.RData`), a `mfp_label` column is
appended to the output.

## Usage

``` r
map_factors_to_metaprograms(
  gsea_study,
  cancer_type,
  mp_file = NULL,
  plot = TRUE,
  file_name = NULL
)
```

## Arguments

- gsea_study:

  Output from compute_factor_gsea() on study cohort.

- cancer_type:

  Character. TCGA cancer type abbreviation (e.g., `"blca"`, `"brca"`,
  `"cesc"`, `"chol"`, `"coad"`, `"skcm"`) identifying which pre-built
  TCGA meta-program file to load from `inst/extdata/`.

- mp_file:

  Optional character. Path to a custom meta-program RData file. If NULL,
  the pre-built file for `cancer_type` is loaded from `inst/extdata/`.

- plot:

  Logical. If TRUE (default), saves a heatmap of factor-to-meta-program
  scores.

- file_name:

  Optional character. File path prefix for saving output plots.

## Value

Data frame with one row per study factor: factor, best_MP, best_score,
all_scores, active_hallmarks
