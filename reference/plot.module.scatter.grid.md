# Save a grid of scatter plots for significant module-feature pairs

For each pair of columns from `matA` and `matB` whose p-value in `p_mat`
is below `pval`, draws a scatter plot with a regression line and
annotates it with the correlation coefficient. The full grid is saved as
an SVG file in `Results/`.

## Usage

``` r
# S3 method for class 'module.scatter.grid'
plot(
  matA,
  matB,
  cor_mat,
  p_mat,
  file_name,
  pval = 0.05,
  cor_type = "p",
  width = width,
  height = width,
  only_sig = TRUE,
  ncol = NULL
)
```

## Arguments

- matA:

  A numeric matrix (samples x modules) for the y-axis of each panel.

- matB:

  A numeric matrix (samples x features) for the x-axis of each panel.

- cor_mat:

  A numeric matrix of correlation coefficients between columns of `matA`
  (rows) and `matB` (columns).

- p_mat:

  A numeric matrix of p-values corresponding to `cor_mat`.

- file_name:

  Character. Base name for the output SVG file (saved to
  `Results/<file_name>_scatter_grid.svg`).

- pval:

  Numeric. P-value cutoff for displaying a pair. Default 0.05.

- cor_type:

  Character. Label used in axis text (e.g., `"p"` for Pearson). Default
  `"p"`.

- width:

  Numeric. Width of the SVG output in inches. Default same as `height`.

- height:

  Numeric. Height of the SVG output in inches. Default same as `width`.

- only_sig:

  Logical. If `TRUE` (default), only pairs with `p <= pval` are plotted.

- ncol:

  Integer or `NULL`. Number of columns in the plot grid. If `NULL`, set
  to `ceiling(sqrt(n_pairs))`.

## Value

Called for its side effect (saves SVG); returns `NULL` invisibly.
