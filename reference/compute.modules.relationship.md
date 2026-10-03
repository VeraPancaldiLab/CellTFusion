# Compute modules relationship

Performs linear correlation between two matrices (e.g., TF modules and
external features such as deconvolution estimates or pathway activities)
across the same set of samples, visualizing only significant
associations.

## Usage

``` r
compute.modules.relationship(
  matA,
  matB,
  file_name,
  batch = NULL,
  width = 8,
  height = 8,
  par_mar = NULL,
  pval = 0.05,
  padj = F,
  cor_type = "p",
  return = F,
  vertical = F,
  plot = T,
  plot.grid = F,
  width.grid = 12,
  height.grid = 10,
  ncol.grid = NULL
)
```

## Arguments

- matA:

  A numeric matrix or data frame of features (samples x features).

- matB:

  A numeric matrix or data frame of features to correlate with (samples
  x features).

- file_name:

  A string indicating the base name (without extension) of the figure to
  be saved in the "Results/" folder.

- batch:

  Optional vector indicating batch assignment for samples. If it has two
  or more levels, Pearson partial correlations
  ([`ppcor::pcor.test()`](https://rdrr.io/pkg/ppcor/man/pcor.test.html))
  are computed controlling for batch and `cor_type` is ignored. Batch is
  treated as categorical (one dummy variable per level beyond the
  first).

- width:

  An integer indicating the width (in inches) of the output PDF figure.
  Default is 8.

- height:

  An integer indicating the height (in inches) of the output PDF figure.
  Default is 8.

- par_mar:

  A numeric vector of length 4 specifying the margin sizes (bottom,
  left, top, right) for the heatmap. If NULL (default), reasonable
  defaults are chosen based on plot orientation.

- pval:

  A numeric threshold for the p-value to define significance. Default is
  0.05.

- padj:

  Logical; if TRUE, applies Bonferroni correction for multiple testing.
  Default is FALSE.

- cor_type:

  Type of correlation to compute: "p" (Pearson), "s" (Spearman), or "k"
  (Kendall). Default is "p".

- return:

  Logical; if TRUE, the function returns a list containing the
  correlation matrix and a named list of significant feature names per
  module. Default is FALSE.

- vertical:

  Logical; if TRUE, modules are on the x-axis and the `matB` features on
  the y-axis. Otherwise (default), the `matB` features are on the x-axis
  and modules on the y-axis.

- plot:

  Logical; if TRUE, saves the heatmap plot as a PDF (only when
  `return = FALSE`). Default is TRUE.

- plot.grid:

  Logical; if TRUE, generates per-pair scatter grid plots for
  significant associations, saved to
  "Results/\<file_name\>\_scatter_grid.svg". Default is FALSE.

- width.grid:

  Numeric width of the scatter grid output (increased if needed to fit
  all panels). Default is 12.

- height.grid:

  Numeric height of the scatter grid output (increased if needed to fit
  all panels). Default is 10.

- ncol.grid:

  Integer number of columns used in scatter grid layout. If NULL
  (default), it is set to the ceiling of the square root of the number
  of panels.

## Value

If `return = TRUE`, returns a list with:

- A correlation matrix between modules and external features.

- A named list with significant features per module (after p-value or
  adjusted p-value thresholding).

If `return = FALSE` and `plot = TRUE`, the function saves a heatmap of
significant correlations to "Results/\<file_name\>.pdf" (the ".pdf"
extension is added if missing).

## Details

The function assumes that `matA` and `matB` share the same rownames
(i.e., samples in the same order); names are compared after
[`make.names()`](https://rdrr.io/r/base/make.names.html), so e.g.
"TCGA-XX" and "TCGA.XX" match. The correlation is computed using WGCNA's
[`cor()`](https://rdrr.io/r/stats/cor.html) and `corPvalueStudent()`
functions. Insignificant correlations (based on p-value or adjusted
p-value) are excluded from the visualization.

## Examples

``` r
data("network.tuto")
data("deconv_subgroups.tuto")
corr = compute.modules.relationship(network.tuto[[1]],
                                    deconv_subgroups.tuto[[1]],
                                    "Deconvolution-TFs_Modules",
                                    plot = FALSE,
                                    return = TRUE,
                                    pval = 0.01)

if (FALSE) { # \dontrun{
# TF modules vs PROGENy pathways (downloads the PROGENy model from OmniPath)
data("counts.norm.tuto")
pathways <- compute.pathway.activity(counts.norm.tuto)
compute.modules.relationship(network.tuto[[1]],
                             pathways,
                             "Pathways_Progeny-TFs_Modules",
                             width = 15)
} # }
```
