# Annotate meta-programs with Bagaev TME subtypes

Assigns each meta-program a TME subtype (IE, IE/F, F, D) by weighted
majority vote across the TCGA factors with a positive score for it: the
scores of the factors are summed per MFP subtype and the subtype with
the highest total is assigned. Factors without a subtype in
`factor_tme_df` do not vote (a message lists them).

## Usage

``` r
annotate_metaprograms_TME(meta_programs_df, factor_tme_df, factors_mp_df)
```

## Arguments

- meta_programs_df:

  Data frame output of
  [`derive_meta_programs()`](https://verapancaldilab.github.io/CellTFusion/reference/derive_meta_programs.md),
  with columns `meta_program` and `hallmarks`.

- factor_tme_df:

  Data frame output of
  [`map_factors_to_TME()`](https://verapancaldilab.github.io/CellTFusion/reference/map_factors_to_TME.md),
  with columns `factor` and `best_MFP`.

- factors_mp_df:

  The `factor_mapping` data frame returned by
  [`map_factors_to_metaprograms()`](https://verapancaldilab.github.io/CellTFusion/reference/map_factors_to_metaprograms.md)
  run on the TCGA factors themselves, with columns `factor` and
  `all_scores`.

## Value

`meta_programs_df`, rows in the same order, with an additional (or
replaced) `TME_subtype` column (`"IE"`, `"IE/F"`, `"F"`, `"D"`, or
`"uncharacterized"` when no factor with a subtype has a positive score
for the meta-program, the vote is tied, or the winning label is
`"uncharacterized"`).

## Examples

``` r
if (FALSE) { # \dontrun{
# On a TCGA cohort: factors, their GSEA, meta-programs, and the MFP subtype of each factor
gsea <- compute_factor_gsea(counts_norm = counts.norm, scores = res$Latent_spaces$Z, return = FALSE)
meta_programs <- derive_meta_programs(gsea)
factors_mp <- map_factors_to_metaprograms(gsea, mp_file = meta_programs,
                                          return = FALSE)$factor_mapping
tme <- map_factors_to_TME(cancer_type = "BLCA", Z = res$Latent_spaces$Z)
meta_programs <- annotate_metaprograms_TME(meta_programs, tme, factors_mp)
} # }
```
