# Map study factors to TCGA meta-programs

For each study NMF factor, scores it against each TCGA meta-program by
computing the mean NES of the meta-program's Hallmarks in that factor.
The meta-program with the highest positive mean NES is the best match.
If the meta-program reference has a `TME_subtype` column (see
[`annotate_metaprograms_TME()`](https://verapancaldilab.github.io/CellTFusion/reference/annotate_metaprograms_TME.md)),
it is appended to the output.

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

  Character. TCGA cancer type abbreviation identifying which pre-built
  TCGA meta-program reference shipped with the package to use.
  Available: `"blca"`, `"luad"`, `"skcm"`. Ignored if `mp_file` is
  given.

- mp_file:

  Optional. Either a meta-program data frame (as returned by
  [`derive_meta_programs()`](https://verapancaldilab.github.io/CellTFusion/reference/derive_meta_programs.md)
  /
  [`annotate_metaprograms_TME()`](https://verapancaldilab.github.io/CellTFusion/reference/annotate_metaprograms_TME.md))
  or the path to an RData file containing such an object named
  `meta_programs`. If NULL, the pre-built reference for `cancer_type` is
  used.

- plot:

  Logical. If TRUE (default), saves a barplot of factor-to-meta-program
  scores to `Results/Factor_MP_mapping_<file_name>.pdf`.

- file_name:

  Optional character suffix for saving output plots.

## Value

A list with:

- `factor_mapping`: Data frame with one row per study factor and columns
  `best_MP`, `factor`, `best_score`, `all_scores` (all meta-program
  scores as `"MP:score"` pairs) and, if available, `TME_subtype`.

- `reference`: The meta-program reference data frame used.
