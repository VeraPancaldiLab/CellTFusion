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
  scores,
  coldata,
  pval = 0.05,
  corr_type = "spearman",
  file_name = NULL,
  ncol_grid = NULL,
  plot_grid = FALSE,
  padj = FALSE
)
```

## Arguments

- scores:

  A numeric matrix or data frame of TF module scores across samples.
  Typically the first element (`TFs module matrix`) of the output from
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md).
  Rows represent samples, columns represent TF modules.

- coldata:

  A data frame containing clinical traits (both categorical and
  numerical) for the same samples. Rows must be the same samples, in the
  same order, as the rows of `scores` (checked by name; an error is
  raised otherwise). Numeric columns are treated as continuous: convert
  number-coded categories (e.g. a stage coded 1 to 4) to factor to test
  them as groups.

- pval:

  A numeric threshold (default = 0.05) to determine statistical
  significance (p-value, or adjusted p-value, `<= pval`). Only
  associations with (adjusted) p-values below this threshold are
  annotated in the heatmap or drawn as boxplots.

- corr_type:

  Correlation used for the numeric traits: `"spearman"` (default) or
  `"pearson"`.

- file_name:

  Optional character suffix of the names of the saved plots (`NULL` by
  default: no suffix).

- ncol_grid:

  Integer. Number of columns in the grid of association boxplots. If
  `NULL` (default), it is set to the ceiling of the square root of the
  number of panels.

- plot_grid:

  Logical; if TRUE, tests categorical traits with ANOVA and saves
  boxplot grids. Default is FALSE.

- padj:

  Logical; if TRUE, the p-values are corrected for multiple testing with
  the Benjamini-Hochberg method, across the modules for each trait, and
  the adjusted p-values are compared with `pval`. Default is FALSE.

## Value

Invisibly, a list with the results:

- `correlation`: matrix of module-trait correlations for the numeric
  traits (modules x traits).

- `p_value`: matrix of their p-values.

- `p_adjusted`: matrix of their Benjamini-Hochberg adjusted p-values
  (equal to `p_value` with `padj = FALSE`).

- `anova`: data frame with the ANOVA of every module for every
  categorical trait (`trait`, `module`, `F`, `p`, `p_adjusted`); `NULL`
  if `plot_grid = FALSE`.

The elements are `NULL` when there is no trait of that type. The
following files are saved to the `Results/` directory:

- `TF.modules_metadata_<file_name>.pdf`: a heatmap of module-trait
  correlations for numeric traits (only associations with p-value \<=
  `pval` are annotated, with the correlation and its significance: \* \<
  0.05, \*\* \< 0.01, \*\*\* \< 0.001, \*\*\*\* \< 0.0001).

- `ANOVA_boxplot_summary_<file_name>_<trait>.svg` (if
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

# Categorical traits (response, Enrollment.IC) are tested with plot_grid = TRUE; os and censOS are numeric
res <- compute.metadata.association(
  scores = network.tuto[[1]],
  coldata = traitdata.tuto,
  pval = 0.05,
  file_name = 'Tutorial',
  plot_grid = TRUE
)
head(res$anova)
#>                             trait module     F     p p_adjusted
#> 1 Best.Confirmed.Overall.Response  black 4.404 0.037      0.037
#> 2 Best.Confirmed.Overall.Response   blue 0.164 0.686      0.686
#> 3 Best.Confirmed.Overall.Response  brown 0.954 0.330      0.330
#> 4 Best.Confirmed.Overall.Response  green 0.807 0.370      0.370
#> 5 Best.Confirmed.Overall.Response   pink 8.859 0.003      0.003
#> 6 Best.Confirmed.Overall.Response    red 1.008 0.317      0.317
```
