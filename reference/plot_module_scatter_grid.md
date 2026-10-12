# Save a grid of scatter plots for significant module-feature pairs

For each pair of columns from `matA` and `matB` whose p-value in `p_mat`
is at or below `pval`, draws a scatter plot with a linear regression
line and annotates it with the correlation coefficient and p-value taken
from `cor_mat` and `p_mat`. The full grid is saved as an SVG file in
`Results/`. The line is a linear fit of the values given in `matA` and
`matB`. Each panel is 3.2 x 3.2 inches, so the size of the file follows
from the number of panels. The column names are the axis titles; names
longer than 30 characters are split over several lines at `"_"` or
`"."`, and the margins are set from the number of lines.

## Usage

``` r
plot_module_scatter_grid(
  matA,
  matB,
  cor_mat,
  p_mat,
  file_name,
  pval = 0.05,
  only_sig = TRUE,
  ncol = NULL,
  max_pairs = 100,
  axis_suffix = ""
)
```

## Arguments

- matA:

  A numeric matrix or data frame (samples x modules) for the y-axis of
  each panel.
  [`compute.modules.relationship()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md)
  passes the ranks for Spearman correlations, and the values corrected
  for batch when `batch` is used.

- matB:

  A numeric matrix or data frame (samples x features) for the x-axis of
  each panel.

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

- only_sig:

  Logical. If `TRUE` (default), only pairs with `p <= pval` are plotted.

- ncol:

  Integer or `NULL`. Number of columns in the plot grid. If `NULL`, set
  to `ceiling(sqrt(n_pairs))`.

- max_pairs:

  Integer. Maximum number of panels: if there are more pairs, only the
  strongest correlations are drawn (with a message). Default 100.

- axis_suffix:

  Character added to the axis titles after the column names (e.g.
  `" (rank)"`). Default "".

## Value

Called for its side effect (saves SVG); returns `NULL` invisibly.
