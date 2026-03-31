# Compute Weighted TF-coactivity Network Analysis (WTCNA)

Construct a weighted signed or unsigned network using TF activity to
cluster protein regulators into modules that share similar activity
patterns. Each TF module will have a sample-level score represented by
the eigenvalue of the module.

## Usage

``` r
# S3 method for class 'WTCNA'
compute(
  TFs.matrix,
  batch = FALSE,
  network.type = "signed",
  clustering.method = "ward.D2",
  minMod = 15,
  corr_mod = 0.9,
  cor_type = "p",
  verbose = F,
  file.name = NULL,
  softPower = NULL,
  return = T
)
```

## Arguments

- TFs.matrix:

  Matrix of TF activity (samples x TFs).

- batch:

  Logical; if TRUE, performs consensus WGCNA across cohorts provided as
  a list.

- network.type:

  Network type: "signed", "unsigned", "signed hybrid", or "distance".
  Default is "signed".

- clustering.method:

  Clustering method for hierarchical clustering. Default is "ward.D2".

- minMod:

  Minimum number of TFs per module. Default is 15.

- corr_mod:

  Correlation threshold (0–1) for merging similar modules. Default is
  0.9.

- cor_type:

  Correlation type for adjacency calculation: "p" (Pearson), "s"
  (Spearman). Default is "p".

- verbose:

  Boolen value to whether print or no the function messages

- file.name:

  Optional character suffix used when writing WTCNA outputs.

- softPower:

  Optional numeric value specifying the soft-thresholding power to be
  used when constructing

- return:

  Logical, whether to save output plots and module list to "Results/".
  Default is TRUE.

## Value

A named list with:

- `TFs module matrix`: Matrix of module eigengenes (samples x modules).

- `TFs colors`: Vector of module colors assigned to each TF.

- `TFs per module`: List of TF names in each module.

- `Proportion of variance`: Matrix of variance explained per module.

## References

Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted
correlation network analysis. BMC Bioinformatics, 9, 559.
https://doi.org/10.1186/1471-2105-9-559

## Examples

``` r
data("tfs.tuto")
network <- compute.WTCNA(tfs.tuto, corr_mod = 0.9, clustering.method = "ward.D2", return = FALSE)
#> Warning: executing %dopar% sequentially: no parallel backend registered
```
