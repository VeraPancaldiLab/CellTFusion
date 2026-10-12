# Compute associations between TF module scores and clinical metadata

This function tests for associations between transcription factor (TF)
module scores and available clinical traits. Numeric traits are
correlated with each module (Pearson or Spearman) and shown as a
heatmap. If `plot_grid = TRUE`, categorical traits (text, factor and
TRUE/FALSE columns) are additionally tested with one-way ANOVA (Tukey
HSD post-hoc) and shown as boxplot grids. All plots are saved in the
`Results/` directory, and the tables with the results are returned.

## Usage

``` r
compute.metadata.association(
  tfs.modules,
  coldata,
  pval = 0.05,
  corr_method = "p",
  file.name,
  ncol = 5,
  plot_grid = F,
  p_adjust = "none"
)
```

## Arguments

- tfs.modules:

  A numeric matrix or data frame of TF module scores across samples.
  Typically the first element (`TFs module matrix`) of the output from
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md).
  Rows represent samples, columns represent TF modules.

- coldata:

  A data frame containing clinical traits (both categorical and
  numerical) for the same samples. Rows must be the same samples, in the
  same order, as the rows of `tfs.modules` (checked by name; an error is
  raised otherwise). Numeric columns are treated as continuous: convert
  number-coded categories (e.g. a stage coded 1 to 4) to factor to test
  them as groups.

- pval:

  A numeric threshold (default = 0.05) to determine statistical
  significance (p-value, or adjusted p-value, `<= pval`). Only
  associations with (adjusted) p-values below this threshold are
  annotated in the heatmap or drawn as boxplots.

- corr_method:

  Character string specifying the correlation method to use for
  continuous variables. Options are `"p"` for pearson and `"s"` for
  spearman. Default is `"p"`.

- file.name:

  Character. Base file name for saving plots of results.

- ncol:

  Integer. Number of columns in the grid of association boxplots.
  Default is 5.

- plot_grid:

  Logical; if TRUE, tests categorical traits with ANOVA and saves
  boxplot grids. Default is FALSE.

- p_adjust:

  Character. Method to correct the p-values for multiple testing,
  applied across the modules for each trait (any method of
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html), e.g.
  `"BH"` or `"bonferroni"`). Default is `"none"`.

## Value

Invisibly, a list with the results:

- `correlation`: matrix of module-trait correlations for the numeric
  traits (modules x traits).

- `p_value`: matrix of their p-values.

- `p_adjusted`: matrix of their p-values adjusted with `p_adjust` (equal
  to `p_value` with `"none"`).

- `anova`: data frame with the ANOVA of every module for every
  categorical trait (`trait`, `module`, `F`, `p`, `p_adjusted`); `NULL`
  if `plot_grid = FALSE`.

The elements are `NULL` when there is no trait of that type. The
following files are saved to the `Results/` directory:

- `TF.modules_metadata_<file.name>.pdf`: a heatmap of module-trait
  correlations for numeric traits (only associations with p-value \<=
  `pval` are annotated, with the correlation and its significance: \* \<
  0.05, \*\* \< 0.01, \*\*\* \< 0.001, \*\*\*\* \< 0.0001).

- `ANOVA_boxplot_summary_<file.name>_<trait>.svg` (if
  `plot_grid = TRUE`): one boxplot grid per categorical trait with the
  significant modules.

## Details

Traits that cannot be tested are skipped with a message: numeric traits
without variation, and categorical traits with a single level or with as
many levels as samples (e.g. a sample identifier). Categorical traits
are only analysed when `plot_grid = TRUE`; a message lists them
otherwise.

The size of the figures is computed from their content (number of
modules and traits for the heatmap, number of panels for the boxplot
grids), and the y-axis of each boxplot panel from its data.

## Examples

``` r

data("network.tuto")
data("traitdata.tuto")

# The clinical traits of the tutorial data are categorical: they are tested with plot_grid = TRUE
res <- compute.metadata.association(
  tfs.modules = network.tuto[[1]],
  coldata = traitdata.tuto,
  pval = 0.05,
  file.name = 'Tutorial',
  plot_grid = TRUE
)
#> No module significantly associated with trait 'Best.Confirmed.Overall.Response'
#> No module significantly associated with trait 'binaryResponse'
#> Error in ggplot2::ggsave(paste0("Results/ANOVA_boxplot_summary", if (!is.null(file.name)) paste0("_",     file.name), "_", gsub("[^A-Za-z0-9._ -]", "_", tr), ".svg"),     plot = p, device = grDevices::svg, width = max(5, 3.2 * min(ncol,         n_panels) + 0.5), height = 3.4 * ceiling(n_panels/ncol) +         1.4, limitsize = FALSE): Cannot find directory Results.
#> ℹ Please supply an existing directory or use `create.dir = TRUE`.
head(res$anova)
#> Error: object 'res' not found
```
