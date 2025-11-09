# Identify hub TFs

Identifies hub TFs per module using values of module membership and
degree.

## Usage

``` r
identify_hub_TFs(datExpr, TF.network, MM_thresh = 0.8, degree_thresh = 0.9)
```

## Arguments

- datExpr:

  A matrix of TF activity (TFs as rows and samples as columns).

- TF.network:

  TF network obtained from compute.WTCNA().

- MM_thresh:

  Threshold for module membership (e.g., 0.8).

- degree_thresh:

  Quantile threshold for degree (e.g., 0.9 for top 10%).

  TFs with high module membership (r \> MM_thresh) and among the top
  percentage of genes by degree (above degree_thresh quantile) are
  considered hub TFs.

## Value

A list with two elements:

- hubGenes:

  Named list of hub TFs per module

- detailedData:

  Dataframe with module, degree, and membership info

## Examples

``` r
data("tfs.tuto")
data("network.tuto")

hub_tfs <- identify_hub_TFs(t(tfs), network, MM_thresh = 0.8, degree_thresh = 0.9)
#> Error: object 'network' not found
```
