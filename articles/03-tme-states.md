# TME State Characterisation

``` r

library(CellTFusion)
#> 
#> 
```

Once latent factors have been extracted (see [Cell Group
Construction](https://verapancaldilab.github.io/CellTFusion/articles/02-cell-groups.md)),
each factor can be functionally annotated by linking it to known
biological programs. This article covers:

1.  **Hallmark GSEA** — associate each latent factor with MSigDB
    Hallmark gene sets
2.  **Meta-program mapping** — compare factors to cancer-type-specific
    reference programs derived from TCGA
3.  **TME subtype annotation** — annotate factors with established TME
    immune subtypes

These steps require:

- `counts.norm` — log-normalized expression matrix (genes × samples)
- `latent_spaces` — output of
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md);
  the `$Z` element is the samples × factors matrix

------------------------------------------------------------------------

## Hallmark GSEA per latent factor

[`compute_factor_gsea()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_factor_gsea.md)
fits a multivariate limma model using the latent factor scores as
covariates, ranks genes by moderated t-statistic per factor, and runs
Hallmark GSEA via `fgsea`. A dot plot of the top enriched pathways is
saved per factor.

``` r

gsea_results <- compute_factor_gsea(
  RNA.tpm     = counts.norm,        # genes × samples expression matrix
  features_df = latent_spaces$Z,    # samples × factors (from compute.latent_factors)
  plot_dot    = TRUE,
  top_n       = 10,
  file_name   = "Tutorial",
  width       = 8,
  height      = 10
)
```

The result is a list with two elements:

| Element         | Description                                            |
|-----------------|--------------------------------------------------------|
| `$DE_results`   | Named list of limma DEG tables, one per factor         |
| `$GSEA_results` | Named list of fgsea result data frames, one per factor |

------------------------------------------------------------------------

## Map factors to cancer-type meta-programs

[`map_factors_to_metaprograms()`](https://verapancaldilab.github.io/CellTFusion/reference/map_factors_to_metaprograms.md)
compares each factor’s Hallmark NES profile against a reference set of
meta-programs derived from TCGA bulk RNA-seq data. The reference covers
multiple cancer types and allows cross-study interpretation of the
identified TME states.

Supported `cancer_type` values: `"skcm"` (melanoma), `"blca"` (bladder
cancer), `"luad"` (lung adenocarcinoma).

``` r

mp_mapping <- map_factors_to_metaprograms(
  gsea_study  = gsea_results,
  cancer_type = "skcm",        # match to your cancer type
  plot        = TRUE,
  file_name   = "Tutorial"
)
```

The result contains:

| Element | Description |
|----|----|
| `$factor_mapping` | Data frame mapping each study factor to the closest reference meta-program |
| `$reference` | The reference meta-program NES matrix used for comparison |

------------------------------------------------------------------------

## Annotate factors with TME immune subtypes

[`map_factors_to_TME()`](https://verapancaldilab.github.io/CellTFusion/reference/map_factors_to_TME.md)
uses TCGA sample-level TME subtype annotations to score each latent
factor for enrichment in specific immune environments (e.g.,
immune-desert, immune-excluded, inflamed).

``` r

tme_annotation <- map_factors_to_TME(
  cancer_name = "skcm",
  Z           = latent_spaces$Z,
  plot        = TRUE,
  file_name   = "Tutorial"
)
```

------------------------------------------------------------------------

## Derive meta-programs from GSEA (unsupervised)

If you have run CellTFusion across multiple cohorts and want to derive
consensus meta-programs from scratch rather than mapping to TCGA
references, use
[`derive_meta_programs()`](https://verapancaldilab.github.io/CellTFusion/reference/derive_meta_programs.md).
This clusters factors by their Hallmark NES profiles to identify
recurring biological programs.

``` r

meta_programs <- derive_meta_programs(
  gsea_results = gsea_results,
  k            = NULL,       # number of clusters; NULL for automatic selection
  file_name    = "Tutorial",
  plot         = TRUE
)
```

------------------------------------------------------------------------

## Putting it all together

In practice, this full annotation sequence is run automatically when
using the
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
wrapper. The outputs are accessible directly from the result object:

``` r

res <- CellTFusion(raw.counts = raw.counts, cancer_type = "skcm", ...)

# Access GSEA + meta-program mapping results
res$TME_states           # factor-to-meta-program mapping table
res$Metaprograms_reference  # reference NES matrix

# Access latent factors used as input to GSEA
head(res$Latent_spaces$Z)
```
