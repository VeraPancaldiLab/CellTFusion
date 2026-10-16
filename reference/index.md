# Package index

## Main

Full pipeline in one call

- [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
  : Compute one-step CellTFusion

## Feature computation

TF activity, TF co-activity modules and pathway activity

- [`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md)
  : Compute Transcription Factor (TF) activity
- [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)
  : Compute Weighted TF-coactivity Network Analysis (WTCNA)
- [`compute.pathway.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md)
  : Compute pathway activity scores
- [`identify_hub_TFs()`](https://verapancaldilab.github.io/CellTFusion/reference/identify_hub_TFs.md)
  : Identify hub TFs
- [`compute.modules.enrichment()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.enrichment.md)
  : Compute TF module enrichment using directed target genes
- [`compute.modules.relationship()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md)
  : Compute modules relationship

## Cell groups and latent factors

Cell group construction, NMF latent factors and cell niches

- [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
  : Construct cell groups based on TF networks and deconvolution
- [`identify(`*`<cell.groups>`*`)`](https://verapancaldilab.github.io/CellTFusion/reference/identify.cell.groups.md)
  : Identify cell groups
- [`cell.groups.computation()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.computation.md)
  : Compute cell group scores from deconvolution and TF module network
- [`compute_composite_score()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_composite_score.md)
  : Compute composite score for cell groups
- [`compute.composition.matrix()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.composition.matrix.md)
  : Compute a cell-type composition matrix from deconvolution subgroups
- [`extract_cells()`](https://verapancaldilab.github.io/CellTFusion/reference/extract_cells.md)
  : Extract cell types from deconvolution feature names
- [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md)
  : Compute latent factors from cell group scores using NMF
- [`compute_cells_niches()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_cells_niches.md)
  : Identify cell niches from NMF latent factors

## Statistical analysis

Associations with clinical variables and survival

- [`scores.stat.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.stat.analysis.md)
  : Perform statistical analysis on scores using a specified test
- [`scores.wilcox.test()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.wilcox.test.md)
  : Wilcoxon rank-sum test for binary traits
- [`scores.ttest()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.ttest.md)
  : T-test for cell group comparisons
- [`scores.kruskal.test()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.kruskal.test.md)
  : Kruskal-Wallis test for multi-group comparisons
- [`scores.anova.test()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.anova.test.md)
  : One-way ANOVA test for multi-group comparisons
- [`scores.fisher.test()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.fisher.test.md)
  : Fisher's exact test for score-trait association
- [`compute.metadata.association()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.metadata.association.md)
  : Compute associations between TF module scores and clinical metadata
- [`compute.survival.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.survival.analysis.md)
  : Kaplan-Meier survival analysis on clinical groups or CellTFusion
  features

## Machine learning and new cohorts

Leakage-aware cross-validation with pipeML and projection of independent
data

- [`prepare_celltfusion_folds()`](https://verapancaldilab.github.io/CellTFusion/reference/prepare_celltfusion_folds.md)
  : Prepare CellTFusion cross-validation folds for pipeML
- [`compute.test.set()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.test.set.md)
  : Compute composite scores on test set based on previous cell groups
- [`project_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/project_factors.md)
  : Project cell group scores onto trained NMF latent factors

## Package data

Example data

- [`raw.counts.tuto`](https://verapancaldilab.github.io/CellTFusion/reference/raw.counts.tuto.md)
  : Raw counts
- [`counts.norm.tuto`](https://verapancaldilab.github.io/CellTFusion/reference/counts.norm.tuto.md)
  : Log(TPM+1) normalized counts
- [`traitdata.tuto`](https://verapancaldilab.github.io/CellTFusion/reference/traitdata.tuto.md)
  : Clinical data
- [`deconv.tuto`](https://verapancaldilab.github.io/CellTFusion/reference/deconv.tuto.md)
  : Example Deconvolution Results
- [`deconv_subgroups.tuto`](https://verapancaldilab.github.io/CellTFusion/reference/deconv_subgroups.tuto.md)
  : Cell subgroups
- [`tfs.tuto`](https://verapancaldilab.github.io/CellTFusion/reference/tfs.tuto.md)
  : TFs data
- [`network.tuto`](https://verapancaldilab.github.io/CellTFusion/reference/network.tuto.md)
  : TF Network
