# Identify hub TFs

Identifies hub TFs per module using values of module membership and
degree. TFs with high module membership (correlation with the module
eigengene \> `MM_thresh`) and a degree at or above the `degree_thresh`
quantile of their module are considered hub TFs. The degree of a TF is
its intramodular connectivity: the sum of its Pearson correlations with
the other TFs of the same module.

## Usage

``` r
identify_hub_TFs(datExpr, TF.network, MM_thresh = 0.8, degree_thresh = 0.9)
```

## Arguments

- datExpr:

  A matrix of TF activity (TFs as rows and samples as columns), with TFs
  in the same order as the module colors of `TF.network` (e.g.
  `t(TF.network$TFs_matrix)`).

- TF.network:

  TF network obtained from compute.WTCNA().

- MM_thresh:

  Threshold for module membership (e.g., 0.8).

- degree_thresh:

  Quantile threshold for degree (e.g., 0.9 for top 10%).

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

hub_tfs <- identify_hub_TFs(t(tfs.tuto), network.tuto, MM_thresh = 0.8, degree_thresh = 0.9)
```
