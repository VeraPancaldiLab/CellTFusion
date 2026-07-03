# Package index

## Main

Compute cell type deconvolution

- [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
  : Compute one-step CellTFusion
- [`compute(`*`<TFs.activity>`*`)`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md)
  : Compute Transcription Factor (TF) activity
- [`compute(`*`<WTCNA>`*`)`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)
  : Compute Weighted TF-coactivity Network Analysis (WTCNA)
- [`compute(`*`<pathway.activity>`*`)`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md)
  : Computes TF-modules pathway activities scores
- [`compute_factor_gsea()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_factor_gsea.md)
  : Run multivariate feature-based GSEA using limma and Hallmark gene
  sets
- [`identify_hub_TFs()`](https://verapancaldilab.github.io/CellTFusion/reference/identify_hub_TFs.md)
  : Identify hub TFs
- [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
  : Construct cell groups based on TF networks and deconvolution
- [`compute(`*`<test.set>`*`)`](https://verapancaldilab.github.io/CellTFusion/reference/compute.test.set.md)
  : Compute composite scores on test set based on previous cell groups

## Utils

Visualization

- [`compute(`*`<metadata.association>`*`)`](https://verapancaldilab.github.io/CellTFusion/reference/compute.metadata.association.md)
  : Compute associations between TF module scores and clinical metadata
- [`compute(`*`<modules.relationship>`*`)`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md)
  : Compute modules relationship
- [`compute(`*`<modules.enrichment>`*`)`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.enrichment.md)
  : Compute TF module enrichment using directed target genes

## Analysis

Cell groups analysis

- [`scores.fisher.test()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.fisher.test.md)
  : Fisher's exact test for score-trait association
- [`scores.anova.test()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.anova.test.md)
  : One-way ANOVA test for multi-group comparisons
- [`scores.wilcox.test()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.wilcox.test.md)
  : Wilcoxon rank-sum test for binary traits
- [`scores.kruskal.test()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.kruskal.test.md)
  : Kruskal-Wallis test for multi-group comparisons
- [`scores.ttest()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.ttest.md)
  : Student's t-test for cell group comparisons
- [`scores.stat.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.stat.analysis.md)
  : Perform statistical analysis on scores using a specified test
- [`compute(`*`<latent_factors>`*`)`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md)
  : Compute latent factors from cell group scores using NMF
- [`identify(`*`<cell.groups>`*`)`](https://verapancaldilab.github.io/CellTFusion/reference/identify.cell.groups.md)
  : Identify cell groups
- [`compute(`*`<survival.analysis>`*`)`](https://verapancaldilab.github.io/CellTFusion/reference/compute.survival.analysis.md)
  : Kaplan-Meier survival analysis on clinical groups or CellTFusion
  features

## Meta-programs and TME

Meta-program derivation and TME mapping

- [`derive_meta_programs()`](https://verapancaldilab.github.io/CellTFusion/reference/derive_meta_programs.md)
  : Derive TME meta-programs by clustering Hallmarks across NMF factors
- [`map_factors_to_metaprograms()`](https://verapancaldilab.github.io/CellTFusion/reference/map_factors_to_metaprograms.md)
  : Map study factors to TCGA meta-programs
- [`map_factors_to_TME()`](https://verapancaldilab.github.io/CellTFusion/reference/map_factors_to_TME.md)
  : Annotate NMF factors with Bagaev et al. (2021) MFP subtypes
- [`annotate_metaprograms_TME()`](https://verapancaldilab.github.io/CellTFusion/reference/annotate_metaprograms_TME.md)
  : Annotate meta-programs with Bagaev TME subtypes
- [`build_nes_matrix()`](https://verapancaldilab.github.io/CellTFusion/reference/build_nes_matrix.md)
  : Build a Hallmarks x factors NES matrix from GSEA results
- [`project_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/project_factors.md)
  : Project cell group scores onto trained NMF latent factors
- [`project_test_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/project_test_factors.md)
  : Project test-set samples onto training NMF factors

## Utils

Internal use (not exported functions)

- [`cell.groups.computation()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.computation.md)
  : Compute cell group scores from deconvolution and TF module network
- [`classify.deconvolution()`](https://verapancaldilab.github.io/CellTFusion/reference/classify.deconvolution.md)
  : Classify samples by high or low deconvolution values in given cell
  groups
- [`compute(`*`<TF.network.classification>`*`)`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TF.network.classification.md)
  : Compute TF Network Classification
- [`compute(`*`<composition.matrix>`*`)`](https://verapancaldilab.github.io/CellTFusion/reference/compute.composition.matrix.md)
  : Compute a cell-type composition matrix from deconvolution subgroups
- [`compute_composite_score()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_composite_score.md)
  : Compute composite score for cell groups
- [`create_tfs_modules()`](https://verapancaldilab.github.io/CellTFusion/reference/create_tfs_modules.md)
  : Create TFs modules
- [`extract_cells()`](https://verapancaldilab.github.io/CellTFusion/reference/extract_cells.md)
  : Extract cells from cell type groups
- [`extract_colors()`](https://verapancaldilab.github.io/CellTFusion/reference/extract_colors.md)
  : Extract colors
- [`extract_wilcox_significant()`](https://verapancaldilab.github.io/CellTFusion/reference/extract_wilcox_significant.md)
  : Extract significant features using Wilcoxon test

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
