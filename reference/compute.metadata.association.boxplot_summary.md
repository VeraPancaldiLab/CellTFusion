# Summarize TF module-trait associations as ANOVA boxplot grids

For each categorical trait in `coldata`, runs a one-way ANOVA for each
TF module, performs Tukey HSD post-hoc tests for significant modules,
and saves a multi-panel SVG boxplot grid to `Results/`.

## Usage

``` r
compute.metadata.association.boxplot_summary(
  tfs.modules,
  coldata,
  pval = 0.05,
  file.name,
  ncol = 5,
  p_adjust = "none"
)
```

## Arguments

- tfs.modules:

  A numeric matrix or data frame of TF module scores (samples x
  modules), typically the first element (`TFs module matrix`) of the
  output from
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md).

- coldata:

  A data frame of sample metadata. Only character, factor and logical
  columns are used as traits. Traits with a single level, or with as
  many levels as samples, are skipped with a message; samples with a
  missing value are left out for that trait.

- pval:

  Numeric. ANOVA (adjusted) p-value threshold for significance
  (`<= pval`). Default 0.05.

- file.name:

  Character. Base name appended to output SVG file names
  (`Results/ANOVA_boxplot_summary_<file.name>_<trait>.svg`, or
  `Results/ANOVA_boxplot_summary_<trait>.svg` without `file.name`). In
  `<trait>`, characters other than letters, digits, ".", "*", "-" and
  spaces (e.g. "/") are replaced by "*".

- ncol:

  Integer. Number of columns in the boxplot facet grid. Default 5. The
  size of the figure is computed from the number of panels, and the
  y-axis of each panel from its data.

- p_adjust:

  Character. Method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) to
  correct the ANOVA p-values across the modules, for each trait. Default
  `"none"`.

## Value

Invisibly, a data frame with the ANOVA of every module for every trait
tested (`trait`, `module`, `F`, `p`, `p_adjusted`), or `NULL` when no
trait could be tested. Called mainly for its side effect (saves one SVG
file per categorical trait with at least one significant module).
