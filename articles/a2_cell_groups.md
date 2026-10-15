# Cell groups and latent factors

``` r

library(CellTFusion)
#> Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
#> 'DelayedArray::makeNindexFromArrayViewport' when loading 'HDF5Array'
#> Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
#> 'DelayedArray::makeNindexFromArrayViewport' when loading 'SummarizedExperiment'
#> 
```

This tutorial covers the core steps of `CellTFusion`: building cell
groups from TF modules and deconvolution features, extracting latent
factors and characterizing cell niches. It uses two objects from the
[Feature
computation](https://verapancaldilab.github.io/CellTFusion/articles/a1_feature_computation.md)
tutorial:

- `network` — TF co-activity modules from
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md).
- `dt` — deconvolution subgroups from
  [`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html).

## **Cell groups**

[`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
correlates each TF module with the deconvolution features, clusters the
significantly correlated features of each module, and cuts the
dendrograms into **cell groups**. Features are clustered on the value of
their correlation with the module, so a cell group is a set of features
with a similar biological activity: features related to the same TF
module to the same degree, which do not need to vary together across
samples. Each cell group covers a single range of correlation values:
the dendrogram is cut with a dynamic tree cut, clusters made of
separated ranges are split, and pieces smaller than `min_cluster_size`
are joined to the neighbouring group with the closest correlation. Each
cell group gets a composite score per sample: the first canonical
component of a canonical correlation analysis (CCA) between the group’s
deconvolution features and the score of its TF module, i.e. the
combination of the features that correlates most with the module score.
The score is oriented to be aligned with the TF module: a positive score
means that the sample is on the active side of the module for that cell
group, a negative score that it is on the inactive side.

Key arguments:

- `network`, `dt` — the inputs above.
- `batch` — optional cohort vector (see [Multi-cohort
  analysis](https://verapancaldilab.github.io/CellTFusion/articles/a6_batch_analysis.md)).
- `pval` — p-value threshold to select the deconvolution features
  significantly correlated with each module.
- `file_name` / `return` — with `return = TRUE` (default), the colored
  dendrograms are saved to
  `Results/Dendrogram_color_clusters_<file_name>.pdf`, and the
  composition and the scores of the cell groups to
  `Results/Cell.groups.composition_<file_name>.csv` and
  `Results/Cell.groups.scores_<file_name>.csv`.

``` r

cell_groups <- construct_cell_groups(
  network           = network,
  dt                = dt,
  pval              = 0.05,
  clustering_method = "ward.D2",
  file_name         = "Tutorial"
)
```

One dendrogram is drawn per TF module, showing how its correlated
deconvolution features cluster into cell groups. The leaves are drawn in
the order of their correlation with the module as far as the tree
allows, and the groups are numbered from the lowest to the highest
correlation. The number written next to each group is the one in the
cell group name (`Dendrogram_<module>_group_<number>`); a group that the
tree splits into two places is drawn, and numbered, in both. Features
that are not assigned to any cluster are not part of any cell group:

![Dendrogram of deconvolution features clustered into cell groups for
one TF module](figures/cell_groups_dendrogram_02.png)

[`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
returns a named list:

| Element | Description |
|----|----|
| `Cell_groups` | Data frame of composite scores (samples x cell groups) |
| `Composition` | Deconvolution features in each cell group |
| `Weights` | CCA projection parameters of each group (`xcoef`, `train_means`, `train_sds`; only `xcoef` with `batch`), used to score new samples |

### Re-running with precomputed features

[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
accepts precomputed `dt`, `tfs` and `pathways`, for example from a
previous run (`res$Processed_deconvolution`, `res$TFs_matrix`,
`res$Pathways_scores`). The corresponding steps are then skipped, which
is useful to try different module or cell group parameters without
recomputing the features:

``` r

res_rerun <- CellTFusion(
  raw_counts = raw.counts,
  normalized = TRUE,
  dt         = res$Processed_deconvolution,
  tfs        = res$TFs_matrix,
  pathways   = res$Pathways_scores,
  corr_mod   = 0.8,    # e.g. merge modules more readily than the default 0.9
  min_module_size     = 10,
  pval       = 0.05,
  return     = TRUE
)
```

## **Latent factors**

[`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md)
decomposes the cell group scores into non-negative latent factors with
NMF. Each score is first split into a positive and a negative part
(`_pos` and `_neg` columns: samples on the active and on the inactive
side of the group’s TF module), so that both directions can contribute
to the factors. If `rank = NULL`, the number of factors is chosen at the
elbow of the reconstruction error across ranks 2 to 8, never more than
the number of cell groups (with fewer than 4 cell groups, `rank = 2` is
used, and `rank = 1` with a single cell group).

``` r

latent_spaces <- compute.latent_factors(
  scores    = cell_groups$Cell_groups,
  rank      = NULL,
  seed      = 123,
  file_name = "Tutorial",
  return    = TRUE
)

head(latent_spaces$Z)  # samples x factors
```

The result contains `Z` (factor scores, samples x factors — the input
for statistics and machine learning), `W` (weights of each cell group
part in each factor), `nmf_input`, the NMF model and the per-sample
factor proportions. Each factor (column of `Z`) sums to 1 over the
samples; the size of each factor is in `nmf_model$d`. The factor
proportions of each sample take these sizes into account. With
`return = TRUE`, the mixture of factors in each sample is plotted:

![Stacked bar plot of latent factor proportions per
sample](figures/nmf_patient_mixture.png)

## **Cell niches**

[`compute_cells_niches()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_cells_niches.md)
characterizes each latent factor by the cell types enriched among the
cell groups that contribute most to it: the top-weighted cell groups are
selected (`quantile_cutoff`) and cell types more frequent in them than
in all cell groups (`enrich_thresh`) are kept. A niche lists the cell
types associated with the factor; it does not yet say whether the factor
has more or fewer of them (the `_pos` and `_neg` parts of the cell
groups are added together).

``` r

cells_niches <- compute_cells_niches(
  latent_factors  = latent_spaces,
  dt              = dt,
  cell_groups     = cell_groups,
  enrich_thresh   = 1.5,
  quantile_cutoff = 0.7,
  return          = TRUE,
  file_name       = "Tutorial"
)
```

For each factor, a star network of its enriched cell types is saved
(edge width reflects the weight of the cell type in the factor):

![Star network of the cell types enriched in latent factor
1](figures/niche_network_factor1.png)

## **Next steps**

With the latent factors you can:

- characterize TME states with Hallmark GSEA and meta-programs — [TME
  state
  characterization](https://verapancaldilab.github.io/CellTFusion/articles/a3_tme_states.md);
- test associations with clinical variables — [Statistical
  analysis](https://verapancaldilab.github.io/CellTFusion/articles/a4_statistical_analysis.md);
- use them as features for prediction — [Machine learning
  workflows](https://verapancaldilab.github.io/CellTFusion/articles/a5_machine_learning.md).
