# CellTFusion

This tutorial demonstrates how to use the `CellTFusion` package and
explains the main functions of the pipeline. We walk through each step
of the pipeline, explain key parameters, and illustrate how to interpret
and save results. For more information about the methods we invite to
read the main paper of the tool.

``` r

library(CellTFusion)
#> 
#> 
```

Load the pre-packaged example data included in CellTFusion:

``` r

raw.counts = CellTFusion::raw.counts.tuto ## Corresponds to bulk RNAseq counts
traitdata = CellTFusion::traitdata.tuto ## Corresponds to the clinical data per samples
```

CellTFusion requires two primary inputs: cell-type deconvolution
proportions and transcription factor (TF) activity scores. To obtain
these features, we will leverage publicly available tools designed for
this purpose.

### **Cell-type deconvolution**

To compute cell-type proportions, we will use the `multideconv` R
package. This tool integrates multiple deconvolution algorithms and
signature matrices to estimate cell-type abundances from bulk RNA-seq
data. For more details, visit the [multideconv GitHub
repository](https://github.com/VeraPancaldiLab/multideconv).

``` r

deconv = multideconv::compute.deconvolution(raw.counts, methods = c("Quantiseq", "Epidish"), normalized = T, return = F)
#> Performing TPM normalization ................................................................................
#> Converting input to matrix.
#> Running deconvolution using the following methods...............................................................
#> 
#> * Quantiseq
#> * Epidish
#> 
#> Running Quantiseq...............................................................
#> 
#> >>> Running quantiseq
#> 
#> Running quanTIseq deconvolution module
#> Gene expression normalization and re-annotation (arrays: FALSE)
#> Removing 17 noisy genes
#> Removing 15 genes with high expression in tumors
#> Signature genes found in data set: 128/138 (92.75%)
#> Mixture deconvolution (method: lsei)
#> Deconvolution successful!
#> 
#> The following method-signature combinations are going to be calculated...............................................................
#> 
#> Methods
#> * Epidish
#> 
#> Signatures
#> * CBSX-HNSCC-scRNAseq
#> * CBSX-Melanoma-scRNAseq
#> * CBSX-NSCLC-PBMCs-scRNAseq
#> * CBSX-Vanderbilt-scRNAseq
#> * CBSX-Zilionis-scRNAseq
#> * LM22
#> * TIL10
#> 
#> Running Epidish...............................................................
#> 
#> 
#> Running Epidish...............................................................
#> 
#> 
#> Running Epidish...............................................................
#> 
#> 
#> Running Epidish...............................................................
#> 
#> 
#> Running Epidish...............................................................
#> 
#> 
#> Running Epidish...............................................................
#> 
#> 
#> Running Epidish...............................................................
#> 
#> Deconvolution cache files removed.
#> Preprocessing deconvolution features...............................................................
#> 
#> Checking consistency in deconvolution cell fractions across patients...............................................................
#> 
#> No extra cell types provided. Only the following cell types will be considered:
#>  B.cells
#> B.naive.cells
#> B.memory.cells
#> Macrophages.cells
#> Macrophages.M0
#> Macrophages.M1
#> Macrophages.M2
#> Monocytes
#> Neutrophils
#> NK.cells
#> NK.activated
#> NK.resting
#> NKT.cells
#> CD4.cells
#> CD4.memory.activated
#> CD4.memory.resting
#> CD4.naive
#> CD8.cells
#> CD4.regulatory
#> CD4.non.regulatory
#> T.cells.helper
#> T.cells.gamma.delta
#> Dendritic.cells
#> Dendritic.activated.cells
#> Dendritic.resting.cells
#> Cancer
#> Endothelial
#> Eosinophils
#> Plasma
#> Myocytes
#> Fibroblast
#> Mast.cells
#> Mast.activated.cells
#> Mast.resting.cells
#> CAF
#> uncharacterized_cell 
#> 
#> If you want to consider other cell types (e.g. from a custom signature) which are not included in the package by default (see README), please provide them in the cells_extra argument.
#> 
#> Total sum across samples of combination Quantiseq_TIL10 is 1
#> Total sum across samples of combination Epidish_CBSX.HNSCC.scRNAseq is 1
#> Total sum across samples of combination Epidish_CBSX.Melanoma.scRNAseq is 1
#> Total sum across samples of combination Epidish_CBSX.NSCLC.PBMCs.scRNAseq is 1
#> Total sum across samples of combination Epidish_CBSX.Vanderbilt.scRNAseq is 1
#> Total sum across samples of combination Epidish_TIL10 is 1
#> Total sum across samples of combination Epidish_CBSX.Zilionis.scRNAseq is 1
#> Total sum across samples of combination Epidish_LM22 is 1
```

``` r

head(deconv[,1:5])
#>                 Quantiseq_TIL10_B.cells Epidish_CBSX.HNSCC.scRNAseq_B.cells
#> SAM7f0d9cc7f001             0.050210988                         0.025626151
#> SAM4305ab968b90             0.008603929                         0.000000000
#> SAMcf018fee2acd             0.042857608                         0.080914523
#> SAMcc4675f394a1             0.025216439                         0.001434758
#> SAM49f9b2e57aa5             0.022668958                         0.137032856
#> SAM2e7aa8fa0ab3             0.012022510                         0.282602048
#>                 Epidish_CBSX.Melanoma.scRNAseq_B.cells
#> SAM7f0d9cc7f001                            0.066233918
#> SAM4305ab968b90                            0.008756584
#> SAMcf018fee2acd                            0.009133098
#> SAMcc4675f394a1                            0.000000000
#> SAM49f9b2e57aa5                            0.036205984
#> SAM2e7aa8fa0ab3                            0.005377105
#>                 Epidish_CBSX.NSCLC.PBMCs.scRNAseq_B.cells
#> SAM7f0d9cc7f001                                 0.9270601
#> SAM4305ab968b90                                 0.7335563
#> SAMcf018fee2acd                                 0.8704203
#> SAMcc4675f394a1                                 0.7806137
#> SAM49f9b2e57aa5                                 0.7666496
#> SAM2e7aa8fa0ab3                                 0.9091358
#>                 Epidish_CBSX.Vanderbilt.scRNAseq_B.cells
#> SAM7f0d9cc7f001                                        0
#> SAM4305ab968b90                                        0
#> SAMcf018fee2acd                                        0
#> SAMcc4675f394a1                                        0
#> SAM49f9b2e57aa5                                        0
#> SAM2e7aa8fa0ab3                                        0
```

### **TF activity inference**

To infer transcription factor activity, we will use the `viper` package
(Alvarez et al. ([2016](#ref-Alvarez2016))) in combination with the
`CollecTRI` regulon database (Müller-Dott et al.
([2023](#ref-10.1093/nar/gkad841))). This approach estimates TF activity
based on the expression levels of their downstream target genes.

``` r

#We first normalize by log2(TPM + 1)
counts.norm = data.frame(ADImpute::NormalizeTPM(raw.counts, log = T))
universe = decoupleR::get_collectri(organism = "human", split_complexes = F)
tfs = compute.TFs.activity(counts.norm, universe = universe)
```

``` r

head(tfs[,1:5])
```

Once features have been calculated we can start analyzing our data using
CellTFusion

Before constructing the cell groups, we first reduce the TF activity
matrix into modules of TFs that exhibit similar activity patterns across
samples. For this, we apply a method inspired by the Weighted Gene
Correlation Network Analysis (WGCNA) approach (Langfelder and Horvath
([2008](#ref-langfelder2008wgcna))), adapted here for TF activity, which
we refer to as Weighted TF Correlation Network Analysis (WTCNA).

``` r

network = compute.WTCNA(tfs, corr_mod = 0.8, clustering.method = "ward.D2", return = T) 
```

To explore how clinical variables are associated with the TF modules, we
can run the following analysis. This function will generate and save the
corresponding plots in the `Results/` directory:

``` r

compute.metadata.association(network[[1]], traitdata, pval = 0.05, file.name = "Tutorial", width = 10) 
```

To functionally characterize the TF modules, we provide a utility to
estimate pathway activities using a multivariate linear model (MLM)
implemented via the `decoupleR` package (Badia-i-Mompel et al.
([2022](#ref-10.1093/bioadv/vbac016))). By default, the `PROGENy`
pathway database (Schubert et al. ([2018](#ref-Schubert2018))) is used;
however, users can supply a custom gene set, in which case Gene Set
Variation Analysis (GSVA) will be applied instead.

``` r

pathways = compute.pathway.activity(counts.norm)
```

To explore the relationship between features, users can run the
following function. This will generate and save a labeled heatmap in the
`Results/` directory, displaying the correlation between the two feature
sets along with significance annotations.

``` r

compute.modules.relationship(network[[1]], pathways, "Pathways_Progeny-TFs_Modules", width = 15)
```

If pathway analysis alone does not provide sufficient insight into the
biological relevance of each module, users can perform
Over-Representation Analysis (ORA) using the Reactome database. The
function below first identifies hub transcription factors (TFs) for each
module based on module membership (MM) and node degree thresholds. It
then performs enrichment analysis using the target genes of these hub
TFs. Resulting dot plots are automatically saved in the `Results/`
directory.

``` r

hub_tfs = identify_hub_TFs(t(tfs), network, MM_thresh = 0.8, degree_thresh = 0.9)
compute.modules.enrichment(counts.norm, hub_tfs)
```

Once the TF modules are defined and functionally characterized, we
further refine the cell-type deconvolution results to reduce
dimensionality and enhance statistical power in downstream association
analyses. To achieve this, we use the `compute.deconvolution.analysis()`
function from the `multideconv` R package. This function identifies
groups of cell types with similar abundance patterns across samples and
returns a reduced feature matrix.

``` r

dt = multideconv::compute.deconvolution.analysis(deconv, corr = 0.7, seed = 123) 
dt = multideconv::deconvolution_dictionary(dt, pathways) ## Apply dictionary of deconvolution (To be added in compute.deconvolution.analysis() soon)
#Association of TFs modules with deconvolution 
compute.modules.relationship(network[[1]], dt[[1]], "Deconvolution-TFs_Modules", vertical = T, height = 30, width = 10, pval = 0.05)
```

### **Construct cell groups**

Cell groups construction is the main part of CellTFusion. Using the
previously calculated features, we identify clusters of cells that share
similar biological activity patterns, reflected by the TF modules scores
across samples.

Supervised analysis

In this approach, we incorporate clinical traits to guide the
identification of cell groups associated with specific phenotypes (e.g.,
treatment response).

``` r

cell_groups = construct_cell_groups(counts.norm, tfs, deconv, network, dt, traitdata, pval = 0.05, 
                                    trait = "Best.Confirmed.Overall.Response", positive = "CR")
```

Unsupervised analysis

Alternatively, we can identify cell groups solely based on molecular
features without using clinical annotations.

``` r

cell_groups = construct_cell_groups(counts.norm, tfs, deconv, network, dt, pval = 0.05)
```

### **Analize results**

To evaluate the association between scores and a clinical trait, you can
perform statistical tests such as Fisher’s exact test, ANOVA, or any
other supported method (Wilcoxon, Kruskal–Wallis, t-test) using the
unified function
[`scores.stat.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.stat.analysis.md).
It accepts cell group scores, NMF latent factors
(`compute.latent.factors()` output), or any samples × features matrix.

``` r

res_fisher <- scores.stat.analysis(scores = cell_groups, coldata = traitdata,
                                   trait = "Best.Confirmed.Overall.Response", method = "fisher", pval = 0.05)

res_anova <- scores.stat.analysis(scores = cell_groups, coldata = traitdata,
                                  trait = "Best.Confirmed.Overall.Response", method = "anova", pval = 0.05)

res_wilcox <- scores.stat.analysis(scores = cell_groups, coldata = traitdata,
                                   trait = "Best.Confirmed.Overall.Response", method = "wilcox", pval = 0.05)

res_kruskal <- scores.stat.analysis(scores = cell_groups, coldata = traitdata,
                                    trait = "Best.Confirmed.Overall.Response", method = "kruskal", pval = 0.05)

# Also works directly with NMF latent factors:
nmf <- compute.latent.factors(cell_groups)
res_nmf <- scores.stat.analysis(scores = nmf, coldata = traitdata,
                                trait = "Best.Confirmed.Overall.Response", method = "anova", pval = 0.05)
```

### **One-step CellTFusion**

For convenience, we offer an all-in-one function that automates the
entire CellTFusion workflow—calculating features, performing
intermediate analyses, and ultimately returning the cell group scores.
To run this streamlined pipeline, simply execute:

``` r

source("R/CellTFusion.R") ## To be added in the package soon
library(dplyr)
res <- CellTFusion(
  raw.counts = raw.counts,
  normalized = TRUE,
  coldata = traitdata, # Optional metadata
  trait = "Best.Confirmed.Overall.Response",  # Optional for supervised analysis
  trait.positive = "CR",                      # Define positive class for trait
  deconv_methods = c("Quantiseq", "Epidish"), # Choose from Quantiseq, Epidish, DeconRNASeq, DWLS, CibersortX
  file_name = "TestRun",
  corr = 0.7,
  pval = 0.05,
  high_corr_groups = 0.85,
  return = FALSE
)
```

### **Machine Learning models**

Below are examples demonstrating how the computed cell groups can be
leveraged to train machine learning models for predicting specific
clinical traits.

We use the R package `pipeML` for model training and prediction. For
more details, visit the [pipeML GitHub
repository](https://github.com/VeraPancaldiLab/pipeML).

To start, we split our dataset into training and testing sets as
follows:

``` r

index = caret::createDataPartition(traitdata[,"Best.Confirmed.Overall.Response"], times = 1, p = 0.8, list = FALSE) 

# Train cohort
traitData_train = traitdata[index, ]
raw.counts_train = raw.counts[,index]

# Test cohort
traitData_test = traitdata[-index, ]
raw.counts_test = raw.counts[,-index]
```

We performed CellTFusion in the training set

``` r

res_training <- CellTFusion(
  raw.counts = raw.counts_train,
  normalized = TRUE,
  coldata = traitData_train, # Optional metadata
  trait = "Best.Confirmed.Overall.Response",  # Optional supervised analysis
  trait.positive = "PD",                      # Define positive class for trait
  deconv_methods = c("Quantiseq", "Epidish"), # Choose from Quantiseq, Epidish, DeconRNASeq, DWLS, CibersortX
  file_name = "TestRun",
  corr = 0.7,
  pval = 0.05,
  high_corr_groups = 0.85,
  return = FALSE
)
```

Next, we train the machine learning model using the cell group scores:

``` r

library(caret)
library(pipeML)
res = pipeML::compute_features.training.ML(features_train = res_training$Cell_groups[[1]], 
                                           target_var = traitData_train$Best.Confirmed.Overall.Response, 
                                           trait.positive = "PD", 
                                           metric = "AUROC", 
                                           task_type = "classification",
                                           stack = F, k_folds = 2, 
                                           n_rep = 2,  
                                           ncores = 2,
                                           return = F)
```

We replicate the trained cell groups into an independent dataset

``` r

## Compute deconvolution in the independent set
deconv_test = multideconv::compute.deconvolution(raw.counts_test, methods = c("Quantiseq", "Epidish"), normalized = T, return = F)

## Replicate cell groups into independent set
testing_set = compute.test.set(res_training$Processed_deconvolution, res_training$Cell_groups, names(res_training$Cell_groups[[2]]), deconv_test)
```

Predicting clinical trait on the test set

``` r

pred = pipeML::compute_prediction(res, testing_set, traitData_test$Best.Confirmed.Overall.Response, "PD", stack = F)
head(pred$Metrics[,1:5])
pred$AUC$AUROC
```

### **Custom k-fold cross validation with CellTFusion**

In certain use cases it is essential to carefully design how the
training and test splits are created to prevent data leakage. This is
especially important for feature construction methods that depend on
sample-level correlations or deconvolution results. If such features are
computed on the full dataset before splitting, information from the test
set can inadvertently influence the training set, leading to inflated
performance estimates.

For this purpose, `pipeML` allows the user to define a custom fold
construction function. This makes it possible to recompute features at
every fold while ensuring that test samples are never included in the
feature learning process.

For this purpose, `CellTFusion` provides the function
`prepare_CellTFusion_folds()`, which is fully compatible with the
`pipeML` package for machine learning model training and evaluation.

The `prepare_CellTFusion_folds()` function handles cross-validation by
recomputing `CellTFusion` features inside each fold. This ensures that:

- Test samples are never included in the feature learning process.

- Hyperparameter tuning can be performed safely without introducing
  bias.

- Parallelization (via foreach and doParallel) allows multiple folds to
  be processed simultaneously, reducing runtime.

When integrated with `pipeML`, the function can be passed as a custom
fold construction function (e.g., `fold_construction_fun` =
`prepare_CellTFusion_folds`), ensuring that the cross-validation setup
mimics real-world application to unseen data.

Here is an example of how to integrate it into the machine learning
workflow:

``` r

universe <- decoupleR::get_collectri(organism = 'human', split_complexes = FALSE)
paths <- decoupleR::get_progeny(organism = 'human', top = 500)

res_groups = compute_features.training.ML(t(raw.counts_train), 
                                          traitData_train$Best.Confirmed.Overall.Response, 
                                          trait.positive = "PD",
                                          metric = "AUROC", 
                                          stack = F, 
                                          k_folds = 3, 
                                          n_rep = 5, 
                                          LODO = F, 
                                          file_name = "Test", 
                                          ncores =  2,                                                                                        
                                          return = T,
                                          fold_construction_fun = prepare_CellTFusion_folds, 
                                          fold_construction_args_fixed = list(deconv = deconv, 
                                                                              universe = universe, 
                                                                              paths = paths, 
                                                                              ncores = 2,
                                                                              normalized = TRUE, 
                                                                              coldata = traitData_train, 
                                                                              trait = "Best.Confirmed.Overall.Response", 
                                                                              trait.positive = "PD"),
                                          fold_construction_args_tunable = list(min_targets_size = c(5, 10, 15, 20), 
                                                                                minMod = c(5, 10, 15, 20), 
                                                                                corr_mod = c(0.7, 0.8, 0.9), 
                                                                                corr = c(0.7, 0.8, 0.9), 
                                                                                high_corr_groups = 0.9))
```

In this setup:

- fold_construction_args_fixed defines parameters that are kept constant
  during the fold construction.

- fold_construction_args_tunable defines the hyperparameters that will
  be tuned within the cross-validation process.

👉 For more details on what these arguments mean and how to set them up,
please refer to the
[pipeML](https://verapancaldilab.github.io/pipeML/articles/pipeML.html)
tutorial.

## References

Alvarez, Mariano J., Yao Shen, Federico M. Giorgi, et al. 2016.
“Functional Characterization of Somatic Mutations in Cancer Using
Network-Based Inference of Protein Activity.” *Nature Genetics* 48 (8):
838–47. <https://doi.org/10.1038/ng.3593>.

Badia-i-Mompel, Pau, Jesús Vélez Santiago, Jana Braunger, et al. 2022.
“decoupleR: Ensemble of Computational Methods to Infer Biological
Activities from Omics Data.” *Bioinformatics Advances* 2 (1): vbac016.
<https://doi.org/10.1093/bioadv/vbac016>.

Langfelder, Peter, and Steve Horvath. 2008. “WGCNA: An r Package for
Weighted Correlation Network Analysis.” *BMC Bioinformatics* 9 (1): 559.
<https://doi.org/10.1186/1471-2105-9-559>.

Müller-Dott, Sophia, Eirini Tsirvouli, Miguel Vazquez, et al. 2023.
“Expanding the Coverage of Regulons from High-Confidence Prior Knowledge
for Accurate Estimation of Transcription Factor Activities.” *Nucleic
Acids Research* 51 (20): 10934–49.
<https://doi.org/10.1093/nar/gkad841>.

Schubert, Michael, Bertram Klinger, Martina Klünemann, et al. 2018.
“Perturbation-Response Genes Reveal Signaling Footprints in Cancer Gene
Expression.” *Nature Communications* 9 (1): 20.
<https://doi.org/10.1038/s41467-017-02391-6>.
