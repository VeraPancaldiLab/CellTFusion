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
  plot = T
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

  Logical; if TRUE, produces a vertical heatmap (traits on x-axis,
  modules on y-axis). Otherwise, a horizontal layout is used. Default is
  FALSE.

- plot:

  Logical; if TRUE, saves the heatmap plot as a PDF. Default is TRUE.

## Value

If `return = TRUE`, returns a list with:

- A correlation matrix between modules and external features.

- A named list with significant features per module (after p-value or
  adjusted p-value thresholding).

If `return = FALSE`, the function saves a heatmap of significant
correlations to "Results/file_name.pdf".

## Details

The function assumes that `matA` and `matB` share the same rownames
(i.e., samples in the same order). The correlation is computed using
WGCNA's [`cor()`](https://rdrr.io/r/stats/cor.html) and
`corPvalueStudent()` functions. Insignificant correlations (based on
p-value or adjusted p-value) are excluded from the visualization.

## Examples

``` r
data("counts.norm.tuto")
data("network.tuto")

pathways <- compute.pathway.activity(counts.norm.tuto)
#> Warning: 'OmnipathR::get_annotation_resources' is deprecated.
#> Use 'annotation_resources' instead.
#> See help("Deprecated")
#> Warning: incomplete final line found on 'https://omnipathdb.org/resources'
#> Warning: 'OmnipathR::import_omnipath_annotations' is deprecated.
#> Use 'annotations' instead.
#> See help("Deprecated")
compute.modules.relationship(network.tuto[[1]],
                             pathways,
                             "Pathways_Progeny-TFs_Modules",
                             width = 15)
#> agg_record_220886245 
#>                    2 

data("deconv_subgroups.tuto")
corr = compute.modules.relationship(network.tuto[[1]],
                                    deconv_subgroups.tuto[[1]],
                                    "Deconvolution-TFs_Modules",
                                    plot = FALSE,
                                    return = TRUE,
                                    pval = 0.01)
```
