# Compute associations between TF module scores and clinical metadata

This function tests for associations between transcription factor (TF)
module scores and available clinical traits. It uses Pearson correlation
for continuous (numeric) traits, and ANOVA for categorical traits.
Results are visualized as a labeled heatmap and violin plots. All plots
are saved in the `Results/` directory.

## Usage

``` r
# S3 method for class 'metadata.association'
compute(
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
  Typically the output from
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

  Character. Base file name for saving PDF plots of results.

- width:

  A numeric value indicating the width (in inches) of the output heatmap
  plot (default = 20).

- height:

  A numeric value indicating the height (in inches) of the output
  heatmap plot (default = 8).

- ncol:

  Integer. Number of columns in the grid of association boxplots.

- y_min:

  Numeric. Lower y-axis limit for grid boxplots.

- y_max:

  Numeric. Upper y-axis limit for grid boxplots.

- plot_grid:

  Logical; if TRUE, generates a grid of boxplot summaries.

- width_grid:

  Numeric width of the grid plot output.

- height_grid:

  Numeric height of the grid plot output.

## Value

This function saves the following to the `Results/` directory:

- A labeled heatmap showing Pearson correlations and ANOVA test
  p-values.

- Individual violin plots for significant categorical trait
  associations.

The function does not return an object to the R environment.

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
