# Statistical Analysis

``` r

library(CellTFusion)
#> 
#> 
```

This article demonstrates how to test associations between CellTFusion
latent factors and clinical variables. The primary input is
`res$Latent_spaces` — the NMF latent factor object returned by
[`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md)
or the
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
wrapper.

Load example data:

``` r

raw.counts <- CellTFusion::raw.counts.tuto
traitdata  <- CellTFusion::traitdata.tuto
```

## Association between latent factors and clinical traits

[`scores.stat.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.stat.analysis.md)
tests the association between latent factor scores and a clinical
variable. It accepts the full `latent_spaces` object (not just `$Z`) as
its `scores` argument.

Supported `method` values: `"fisher"`, `"wilcox"`, `"anova"`,
`"kruskal"`, `"ttest"`.

``` r

# Run CellTFusion to obtain latent spaces
res <- CellTFusion(
  raw.counts  = raw.counts,
  normalized  = TRUE,
  cancer_type = "skcm",
  return      = TRUE
)

# Wilcoxon rank-sum test (binary response variable)
sig_factors <- scores.stat.analysis(
  scores  = res$Latent_spaces,
  coldata = traitdata,
  trait   = "Best.Confirmed.Overall.Response",
  method  = "wilcox",
  pval    = 0.05
)

# Fisher's exact test
sig_factors_fisher <- scores.stat.analysis(
  scores  = res$Latent_spaces,
  coldata = traitdata,
  trait   = "Best.Confirmed.Overall.Response",
  method  = "fisher",
  pval    = 0.05
)

# ANOVA (multi-level categorical variable)
sig_factors_anova <- scores.stat.analysis(
  scores  = res$Latent_spaces,
  coldata = traitdata,
  trait   = "Best.Confirmed.Overall.Response",
  method  = "anova",
  pval    = 0.05
)
```

The function returns a list of significant factor results along with
visualisations (box plots, violin plots) saved to `Results/`.

## Direct access to factor scores

If you want to work with the factor scores directly (e.g., for custom
plots or downstream modelling), access the `$Z` matrix:

``` r

# samples x factors matrix
factor_scores <- data.frame(res$Latent_spaces$Z)
head(factor_scores)
```

## TF module — clinical trait association

To visualise associations between TF module eigengenes (rather than
latent factors) and all available clinical traits simultaneously, use
[`compute.metadata.association()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.metadata.association.md).
It runs Pearson correlation for continuous traits and ANOVA for
categorical traits, and saves a labelled heatmap and violin plots to
`Results/`:

``` r

compute.metadata.association(
  tfs.modules = res$TF_network[[1]],
  coldata     = traitdata,
  pval        = 0.05,
  file.name   = "Tutorial",
  width       = 10
)
```

## TF modules — pathway correlation

To explore the relationship between TF module scores and pathway
activities:

``` r

compute.modules.relationship(
  matA      = res$TF_network[[1]],
  matB      = res$Pathways_scores,
  file_name = "Pathways_vs_TF_modules",
  width     = 15
)
```
