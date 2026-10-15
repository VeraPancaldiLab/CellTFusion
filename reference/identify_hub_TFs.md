# Identify hub TFs

Identifies hub TFs per module using values of module membership and
degree. TFs with high module membership (correlation of their activity
with the module score \> `mm_thresh`, with the correlation used to build
the network, `corr_type` of
[`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md))
and a degree at or above the `degree_thresh` quantile of their module
are considered hub TFs. The degree of a TF is its intramodular
connectivity, as in WGCNA: the sum of its network adjacencies with the
other TFs of the same module. The adjacencies are computed as in
[`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md),
with the same soft-thresholding power, network type and correlation (for
a `"signed"` network, \\((1 + r)/2)^{power}\\, so a negative correlation
adds almost nothing instead of being subtracted).

## Usage

``` r
identify_hub_TFs(network, mm_thresh = 0.8, degree_thresh = 0.9)
```

## Arguments

- network:

  TF network obtained from
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md).

- mm_thresh:

  Threshold for module membership. Default is 0.8.

- degree_thresh:

  Quantile threshold for degree within each module. Default is 0.9 (top
  10%).

## Value

A list with two elements:

- hubGenes:

  Named list of hub TFs per module. Modules without hub TFs are left out
  (a message lists them).

- detailedData:

  Data frame with one row per hub TF (TF names as row names) and the
  columns `Module`, `ModuleMembership` and `Degree`.

## Details

The TF activities are taken from `network` (element `TFs_matrix`) and
matched to the module scores by sample name. For a multi-cohort network
([`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)
with `batch`), the per-cohort TF activities are stacked to compute the
membership, and the degree is the mean of the per-cohort degrees (each
with the power of its cohort).

Networks built before the soft-thresholding power was stored in
`network` (element `Network parameters`) have their power chosen again
with the rule of
[`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)
(with a message), assuming a signed network and Pearson correlation.

## Examples

``` r

data("network.tuto")

hub_tfs <- identify_hub_TFs(network.tuto, mm_thresh = 0.8, degree_thresh = 0.9)
```
