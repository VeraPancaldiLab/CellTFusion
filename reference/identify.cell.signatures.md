# Identify cell presence scores across important features from the trained machine learning models.

This function extracts cell subgroup compositions from a machine
learning model, identifies top predictive features based on variable
importance, and analyzes cell-type presence patterns using clustering.
It computes feature importance based on clustering impact and generates
visualizations for cell feature presence and TF module scores.

## Usage

``` r
# S3 method for class 'cell.signatures'
identify(
  cell_groups,
  deconvolution_processed,
  TF_network,
  deconvolution,
  var_importance,
  n_top = 20,
  sign,
  file.name
)
```

## Arguments

- cell_groups:

  A list containing the cell groups returned by
  [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md).

- deconvolution_processed:

  A list containing the subgroupped deconvolution results returned by
  `compute.deconvolution.analysis()` from the multideconv R package.

- TF_network:

  A list containing the TF modules network returned by
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md).

- deconvolution:

  A data.frame or matrix of deconvolution features where columns
  represent cell types or subgroups.

- var_importance:

  A data.frame or tibble containing variable importance scores (e.g.,
  SHAP values) for the model features.

- n_top:

  Integer. Number of top features to select based on importance and
  direction. Defaults to 20.

- sign:

  Character. Direction of effect to select features: "Increase" or
  "Decrease".

- file.name:

  Character. Base file name for saving PDF plots of results.

## Value

No return value. The function saves several PDF plots into the "Results"
directory and prints progress information during execution.

## Details

The function performs the following major steps:

1.  Extracts deconvolution subgroups composition from the model.

2.  Selects the top features based on average variable importance and
    direction.

3.  Constructs a presence matrix indicating cell type presence across
    top features.

4.  Performs hierarchical clustering on the presence matrix using
    Jaccard distance.

5.  Estimates the optimal number of clusters via silhouette analysis.

6.  Evaluates feature importance based on impact on clustering quality
    via permutation bootstrapping.

7.  Generates heatmaps and bar plots summarizing cell presence scores
    and TF module scores.
