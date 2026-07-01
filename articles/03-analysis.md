# Statistical Analysis

``` r

library(CellTFusion)
#> 
#> 
```

Once cell groups have been constructed (see [Cell Group
Construction](https://verapancaldilab.github.io/CellTFusion/articles/02-cell-groups.md)),
you can evaluate their association with clinical traits using a range of
statistical tests.

## Statistical association tests

The unified function
[`scores.stat.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.stat.analysis.md)
accepts:

- Cell group scores from
  [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
- NMF latent factors from
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md)
- Any samples × features numeric matrix

Supported methods are `"fisher"`, `"anova"`, `"wilcox"`, `"kruskal"`,
and `"ttest"`.

``` r

# Fisher's exact test
res_fisher <- scores.stat.analysis(
  scores  = cell_groups,
  coldata = traitdata,
  trait   = "Best.Confirmed.Overall.Response",
  method  = "fisher",
  pval    = 0.05
)

# ANOVA
res_anova <- scores.stat.analysis(
  scores  = cell_groups,
  coldata = traitdata,
  trait   = "Best.Confirmed.Overall.Response",
  method  = "anova",
  pval    = 0.05
)

# Wilcoxon rank-sum test
res_wilcox <- scores.stat.analysis(
  scores  = cell_groups,
  coldata = traitdata,
  trait   = "Best.Confirmed.Overall.Response",
  method  = "wilcox",
  pval    = 0.05
)

# Kruskal–Wallis test
res_kruskal <- scores.stat.analysis(
  scores  = cell_groups,
  coldata = traitdata,
  trait   = "Best.Confirmed.Overall.Response",
  method  = "kruskal",
  pval    = 0.05
)
```

## Latent factor analysis

[`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md)
decomposes cell group scores into NMF latent factors, which can then be
passed directly to
[`scores.stat.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.stat.analysis.md):

``` r

nmf <- compute.latent_factors(cell_groups)

res_nmf <- scores.stat.analysis(
  scores  = nmf,
  coldata = traitdata,
  trait   = "Best.Confirmed.Overall.Response",
  method  = "anova",
  pval    = 0.05
)
```

## Clinical metadata association

To explore associations between TF module scores and all available
clinical traits simultaneously, use
[`compute.metadata.association()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.metadata.association.md).
It applies Pearson correlation for continuous traits and ANOVA for
categorical traits, and saves a labeled heatmap and violin plots to
`Results/`:

``` r

compute.metadata.association(
  tfs.modules = network[[1]],
  coldata     = traitdata,
  pval        = 0.05,
  file.name   = "Tutorial",
  width       = 10
)
```
