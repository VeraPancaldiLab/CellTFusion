# Summarize TF module-trait associations as ANOVA boxplot grids

For each categorical trait in `coldata`, runs a one-way ANOVA for each
TF module, performs Tukey HSD post-hoc tests for significant modules,
and saves a multi-panel SVG boxplot grid to `Results/`.

## Usage

``` r
# S3 method for class 'metadata.association.boxplot_summary'
compute(
  tfs.modules,
  coldata,
  pval = 0.05,
  file.name,
  ncol = 5,
  y_min = 0,
  y_max = 0.5,
  width = 18,
  height = 10
)
```

## Arguments

- tfs.modules:

  A numeric matrix or data frame of TF module scores (samples x
  modules), typically from
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md).

- coldata:

  A data frame of sample metadata. Only character and factor columns are
  used as traits.

- pval:

  Numeric. ANOVA p-value threshold for significance. Default 0.05.

- file.name:

  Character. Base name appended to output SVG file names.

- ncol:

  Integer. Number of columns in the boxplot facet grid. Default 5.

- y_min:

  Numeric. Lower y-axis limit for boxplots. Default 0.

- y_max:

  Numeric. Upper y-axis limit for boxplots. Default 0.5.

- width:

  Numeric. Width of the SVG output in inches. Default 18.

- height:

  Numeric. Height of the SVG output in inches. Default 10.

## Value

Called for its side effect (saves SVG files); returns `NULL` invisibly
when no significant traits are found.
