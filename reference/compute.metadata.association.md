# Compute associations between TF module scores and clinical metadata

This function tests for associations between transcription factor (TF)
module scores and available clinical traits. Numeric traits are
correlated with each module (Pearson or Spearman) and shown as a labeled
heatmap. If `plot_grid = TRUE`, categorical traits are additionally
tested with one-way ANOVA (Tukey HSD post-hoc) and shown as boxplot
grids. All plots are saved in the `Results/` directory.

## Usage

``` r
compute.metadata.association(
  tfs.modules,
  coldata,
  pval = 0.05,
  corr_method = "p",
  file.name,
  width = 20,
  height = 8,
  ncol = 5,
  y_min = 0,
  y_max = 0.5,
  plot_grid = F,
  width_grid = 18,
  height_grid = 10
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
  numerical) for the same samples. Row names should match those of
  `tfs.modules`.

- pval:

  A numeric threshold (default = 0.05) to determine statistical
  significance. Only associations with p-values below this threshold are
  considered significant in the heatmap.

- corr_method:

  Character string specifying the correlation method to use for
  continuous variables. Options are `"p"` for pearson and `"s"` for
  spearman. Default is `"p"`.

- file.name:

  Character. Base file name for saving plots of results.

- width:

  A numeric value indicating the width (in inches) of the output heatmap
  plot (default = 20).

- height:

  A numeric value indicating the height (in inches) of the output
  heatmap plot (default = 8).

- ncol:

  Integer. Number of columns in the grid of association boxplots.
  Default is 5.

- y_min:

  Numeric. Lower y-axis limit for grid boxplots. Default is 0.

- y_max:

  Numeric. Upper y-axis limit for grid boxplots. Default is 0.5.

- plot_grid:

  Logical; if TRUE, tests categorical traits with ANOVA and saves
  boxplot grids. Default is FALSE.

- width_grid:

  Numeric width (in inches) of the grid plot output. Default is 18.

- height_grid:

  Numeric height (in inches) of the grid plot output. Default is 10.

## Value

Called for its side effects. Saves to the `Results/` directory:

- `TF.modules_metadata_<file.name>.pdf`: a labeled heatmap of
  module-trait correlations for numeric traits (only associations with
  p-value \< `pval` are annotated).

- `ANOVA_boxplot_summary_<file.name>_<trait>.svg` (if
  `plot_grid = TRUE`): one boxplot grid per categorical trait with the
  significant modules.

## Examples

``` r

data("network.tuto")
data("traitdata.tuto")

compute.metadata.association(
  tfs.modules = network.tuto[[1]],
  coldata = traitdata.tuto,
  pval = 0.05,
  file.name = 'Tutorial',
  width = 15,
  height = 10
)
```
