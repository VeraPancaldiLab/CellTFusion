# Annotate meta-programs with Bagaev TME subtypes

Assigns each meta-program a TME subtype (IE, IE/F, F, D) by majority
vote across the TCGA factors that were mapped to it.

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

  Data frame output of
  [`map_factors_to_metaprograms()`](https://verapancaldilab.github.io/CellTFusion/reference/map_factors_to_metaprograms.md)
  run on the TCGA factors themselves, with columns `factor` and
  `best_MP`.

## Value

`meta_programs_df` with an additional `TME_subtype` column.
