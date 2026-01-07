# Package index

## Main

Compute cell type deconvolution

- [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
  : Compute one-step CellTFusion
- [`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md)
  : Compute Transcription Factor (TF) activity
- [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)
  : Compute Weighted TF-coactivity Network Analysis (WTCNA)
- [`compute.pathway.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md)
  : Computes TF-modules pathway activities scores
- [`identify_hub_TFs()`](https://verapancaldilab.github.io/CellTFusion/reference/identify_hub_TFs.md)
  : Identify hub TFs
- [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
  : Construct cell groups based on TF networks and deconvolution
- [`compute.test.set()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.test.set.md)
  : Compute composite scores on test set based on previous cell groups

## Utils

Visualization

- [`compute.metadata.association()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.metadata.association.md)
  : Compute associations between TF module scores and clinical metadata
- [`compute.modules.relationship()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md)
  : Compute modules relationship
- [`compute.modules.enrichment()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.enrichment.md)
  : Compute TF module enrichment using directed target genes

## Analysis

Cell groups analysis

- [`cell.groups.fisher.test()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.fisher.test.md)
  : Fisher test using cell groups scores
- [`cell.groups.anova.test()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.anova.test.md)
  : One-way ANOVA test for multi-group comparisons
- [`cell.groups.wilcox.test()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.wilcox.test.md)
  : Wilcoxon rank-sum test for binary traits
- [`cell.groups.kruskal.test()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.kruskal.test.md)
  : Kruskal–Wallis test for multi-group comparisons
- [`cell.groups.ttest()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.ttest.md)
  : Student's t-test for cell group comparisons
- [`cell.groups.stat.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.stat.analysis.md)
  : Perform statistical analysis on cell group scores using a specified
  test
- [`identify(`*`<cell.groups>`*`)`](https://verapancaldilab.github.io/CellTFusion/reference/identify.cell.groups.md)
  : Identify cell groups
- [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md)
  : Full NMF pipeline for latent immune states (single cohort)
- [`identify(`*`<cell.signatures>`*`)`](https://verapancaldilab.github.io/CellTFusion/reference/identify.cell.signatures.md)
  : Identify cell presence scores across important features from the
  trained machine learning models.

## Machine learning

- [`prepare_CellTFusion_folds()`](https://verapancaldilab.github.io/CellTFusion/reference/prepare_CellTFusion_folds.md)
  : Prepare CellTFusion folds for cross-validation with training and
  test data

## Utils

Internal use (not exported functions)

- [`calculate_dendrogram_cuts()`](https://verapancaldilab.github.io/CellTFusion/reference/calculate_dendrogram_cuts.md)
  : Calculate dendrogram cuts
- [`cell.groups.computation()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.computation.md)
  : Compute cell group scores from deconvolution and TF module network
- [`classify.deconvolution()`](https://verapancaldilab.github.io/CellTFusion/reference/classify.deconvolution.md)
  : Classify samples by high or low deconvolution values in given cell
  groups
- [`compute.TF.network.classification()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TF.network.classification.md)
  : Compute TF Network Classification
- [`compute.composition.matrix()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.composition.matrix.md)
  : Compute a cell-type composition matrix from deconvolution subgroups
- [`compute_cell_groups_signatures()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_cell_groups_signatures.md)
  : Compute projected cell group scores on an independent cohort
- [`compute_composite_score()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_composite_score.md)
  : Compute composite score for cell groups
- [`correlation()`](https://verapancaldilab.github.io/CellTFusion/reference/correlation.md)
  : Perform pairwise correlation across all features
- [`create_tfs_modules()`](https://verapancaldilab.github.io/CellTFusion/reference/create_tfs_modules.md)
  : Create TFs modules
- [`extract_cells()`](https://verapancaldilab.github.io/CellTFusion/reference/extract_cells.md)
  : Extract cells from cell type groups
- [`extract_colors()`](https://verapancaldilab.github.io/CellTFusion/reference/extract_colors.md)
  : Extract colors
- [`extract_wilcox_significant()`](https://verapancaldilab.github.io/CellTFusion/reference/extract_wilcox_significant.md)
  : Extract significant features using Wilcoxon test
- [`find.maximum.iteration()`](https://verapancaldilab.github.io/CellTFusion/reference/find.maximum.iteration.md)
  : Find maximum iteration from subgroups
- [`mergeModules()`](https://verapancaldilab.github.io/CellTFusion/reference/mergeModules.md)
  : Merge TFs modules
- [`module_enrich()`](https://verapancaldilab.github.io/CellTFusion/reference/module_enrich.md)
  : Module enrichment
- [`remove.cell.groups.corr()`](https://verapancaldilab.github.io/CellTFusion/reference/remove.cell.groups.corr.md)
  : Remove highly correlated cell groups
- [`remove_equal()`](https://verapancaldilab.github.io/CellTFusion/reference/remove_equal.md)
  : Remove cell groups with equal composition
- [`remove_single_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/remove_single_groups.md)
  : Remove cell groups with only one feature

## Package Data

Example data

- [`raw.counts.tuto`](https://verapancaldilab.github.io/CellTFusion/reference/raw.counts.tuto.md)
  : Raw counts
- [`traitdata.tuto`](https://verapancaldilab.github.io/CellTFusion/reference/traitdata.tuto.md)
  : Clinical data
- [`tfs.tuto`](https://verapancaldilab.github.io/CellTFusion/reference/tfs.tuto.md)
  : TFs data
- [`counts.norm.tuto`](https://verapancaldilab.github.io/CellTFusion/reference/counts.norm.tuto.md)
  : Log(TPM+1) normalized counts
- [`network.tuto`](https://verapancaldilab.github.io/CellTFusion/reference/network.tuto.md)
  : TF Network
- [`deconv.tuto`](https://verapancaldilab.github.io/CellTFusion/reference/deconv.tuto.md)
  : Example Deconvolution Results
- [`deconv_subgroups.tuto`](https://verapancaldilab.github.io/CellTFusion/reference/deconv_subgroups.tuto.md)
  : Cell subgroups
