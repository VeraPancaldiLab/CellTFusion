# Identify cell groups

Builds, for each TF module, a dendrogram of the deconvolution features
significantly correlated with it. The features are clustered on the
value of that correlation: the distance between two features is the
absolute difference between their correlations with the module (e.g.
0.62 and 0.58 are at distance 0.04). Two features are therefore similar
when they are related to the module to the same degree (similar
biological activity), not when they vary together across samples. The
dendrograms are cut into cell groups by
[`cell.groups.computation()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.computation.md).

## Usage

``` r
# S3 method for class 'cell.groups'
identify(features, clustering.method = "ward.D2", file_name = NULL, return = T)
```

## Arguments

- features:

  The list returned by
  [`compute.modules.relationship()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md).
  Its first two elements are used (accessed by position):

  - `[[1]]` (`correlation`): A matrix of correlations between TF modules
    (rows, named) and deconvolution features (columns).

  - `[[2]]` (`significant`): A list, named by TF module, of the features
    significantly correlated with each module (p-value \<= the `pval`
    given to
    [`compute.modules.relationship()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md)).

  The correlations of each module are taken from the row with its name.

- clustering.method:

  Linkage method of the hierarchical clustering
  ([`stats::hclust()`](https://rdrr.io/r/stats/hclust.html)) of the
  features on the distance described above. Default: "ward.D2".

- file_name:

  Optional character suffix for the names of the saved plots.

- return:

  Logical; whether to save one dendrogram plot per module to
  `Results/Dendrogram_cell_types_<module>_<file_name>.pdf`
  (`Results/Dendrogram_cell_types_<module>.pdf` without `file_name`).
  Each plot shows the correlation of every feature with the module next
  to its name, and its size follows the number of features. Default:
  TRUE.

## Value

A named list of dendrograms (`hclust` objects, labelled with the feature
names), one per TF module. Each one also stores the correlation of its
features with the module (element `correlation`, used by
[`calculate_dendrogram_cuts()`](https://verapancaldilab.github.io/CellTFusion/reference/calculate_dendrogram_cuts.md)).
Modules with fewer than two significant features are discarded; `NULL`
is returned if no module remains.

## Examples

``` r
data("network.tuto")
data("deconv_subgroups.tuto")
rel <- compute.modules.relationship(network.tuto[[1]], deconv_subgroups.tuto[[1]], plot = FALSE)
dendrograms <- identify.cell.groups(rel, return = FALSE)
names(dendrograms)
#>  [1] "brown"     "purple"    "blue"      "pink"      "magenta"   "red"      
#>  [7] "green"     "yellow"    "black"     "turquoise"
```
