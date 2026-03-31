# CellTFusion

This tutorial demonstrates how to use the `CellTFusion` package and
explains the main functions of the pipeline. We walk through each step
of the pipeline, explain key parameters, and illustrate how to interpret
and save results. For more information about the methods we invite to
read the main paper of the tool.

``` r
library(CellTFusion)
#> 
#> Warning: replacing previous import 'dplyr::select' by 'AnnotationDbi::select'
#> when loading 'CellTFusion'
#> Warning: replacing previous import 'AnnotationDbi::select' by 'dplyr::select'
#> when loading 'CellTFusion'
#> 
#> Warning: replacing previous import 'dendextend::cutree' by 'stats::cutree' when
#> loading 'CellTFusion'
```

Load the pre-packaged example data included in CellTFusion:

``` r
raw.counts = CellTFusion::raw.counts.tuto ## Corresponds to bulk RNAseq counts
traitdata = CellTFusion::traitdata.tuto ## Corresponds to the clinical data per sample
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
#> * BPRNACan
#> * BPRNACan3DProMet
#> * BPRNACanProMet
#> * CBSX-HNSCC-scRNAseq
#> * CBSX-Melanoma-scRNAseq
#> * CBSX-NSCLC-PBMCs-scRNAseq
#> * CBSX-NSCLC-scRNAseq
#> * CCLE-TIL10
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
#> 
#> Running Epidish...............................................................
#> 
#> 
#> Running Epidish...............................................................
#> 
#> 
#> Running Epidish...............................................................
#> 
#> Preprocessing deconvolution features...............................................................
#> 
#> Checking consistency in deconvolution cell fractions across patients...............................................................
#> 
#> 
#> 
#> Total sum across samples of combination Quantiseq is not 1! Remember these are proportions and the total should be 1
#> Samples which sum with combination Quantiseq is not 1:
#> 
#>  SAM7f0d9cc7f001
#> SAMcf018fee2acd
#> SAM49f9b2e57aa5
#> SAMdf3e42c8672a
#> SAMd027124354ce
#> SAMe7bf6c015192
#> SAM6dd7ad1d797d
#> SAMc692536a795a
#> SAM9a2cf3c06fb3
#> SAM557dde1b9f3e
#> SAM23aa15d4a0b0
#> SAMb963dda93cfd
#> SAMbcbc7957c264
#> SAM7fb6987514a4
#> SAM63405b04ab2d
#> SAM18bc1078bc15
#> SAMd1bd63734394
#> SAMe9ae8beb82fa
#> SAMbe83eae4026e
#> SAMe5bc41772bc9
#> SAM23095936e611
#> SAMdb3f50c9129c
#> SAMbf1a3ae828e6
#> SAM52e3fa3ad574
#> SAMd4c0837b0997
#> SAM9cafb905b36a
#> SAM032c642382a7
#> SAM0ce9c983b20f
#> SAM2f228939632f
#> SAM297c0301e861
#> SAM1fa6bcb7fc48
#> SAM6d2ae0c39b96
#> SAMeff2ce356ccb
#> SAM110501d0eedb
#> SAM2e9ac0b1b250
#> SAMc0ef41aa6c8b
#> SAMd636e3461955
#> SAMd98bac0a070f
#> SAM8e8ef2368dfa
#> SAMb4c7a001537d
#> SAMfd947610629d
#> SAM943df5cf15df
#> SAM39eb94fa504d
#> SAM62fb1388c871
#> SAM2dc3f04e45e9
#> SAM4501e41e4751
#> SAM5b57e47fdcb3
#> SAMae4da274eded
#> SAMd35318127278
#> SAM75142fcab9df
#> SAM166a419a4e5a
#> SAM025b45c27e05
#> SAM560f23d6a3ad
#> SAMf3a9bce50099
#> SAMdab9ca8fb5de
#> SAM34430ef08e5b
#> SAM4b0175e8db6e
#> SAM5e3bae090b8c
#> SAM3cb94b0d5297
#> SAMcb132b0cdd2c
#> SAMe97af0feefdf
#> SAMee3844cc0b9f
#> SAMdad5c29dc105
#> SAMa424c75831b4
#> SAM45c8e6412c66
#> SAM73663ee4a96e
#> SAM0bdb3428bd13
#> SAM3330c03fdf00
#> SAM0a7c2091dd56
#> SAM99a46b9eec27
#> SAM553c3c35b847
#> SAMc2a1820d4e6b
#> SAMe1eb5d988760
#> SAMe3210d3632b4
#> SAMb15ad09d6e24
#> SAM7893196e0e89
#> SAM961d04c42bd9
#> SAMb0d11db9aa79
#> SAMa90d73f8d891
#> SAM415f36ad349e
#> SAMad83c9c53537
#> SAM3ee5dcd894f0
#> SAMd7d57ee3a863
#> SAMeb29625f76a5
#> SAM563d6233dfa2
#> SAMd135d5867fe3
#> SAM065890737112
#> SAMb470eb8f04be
#> SAM675a12a09c15
#> SAM1e9c4d1d39ae
#> SAM4b7ea015fd9e
#> SAM9306c5c92444
#> SAM63b2189c36d7
#> SAM18be5b395318
#> SAM6cbc10abddb0
#> SAMdee1011782cd
#> SAM203dcf14f927
#> SAMe41b1e773582
#> SAM978a587b207e
#> SAM5234688806a7
#> SAM2c9586161ce6
#> SAM76a431ba6ce1
#> SAMd3bd67996035
#> SAMd3601288319e
#> SAMba1a34b5a060
#> SAM18a4dabbc557
#> SAMfed609955db9
#> SAMf2aae1443f67
#> SAM2bba8cb35e48
#> SAM5c139c5c1c4f
#> SAM6780ed436b55
#> SAM85f0a3ac1c45
#> SAM49d48750e294
#> SAM7aa01fc49a80
#> SAMa321770ac31c
#> SAM3894ac3956a5
#> SAM1c0ecfb3eb63
#> SAMbcb07ba81cee
#> SAM9aa6a095a9d6
#> SAM957378bd907f
#> SAM2624229effe8
#> SAMd2492b2a31bb
#> SAM670649e105b5
#> SAM5fe7a81a39dd
#> SAM29da928587ad
#> SAM0a0f2bac4b20
#> SAM8e469834acc1
#> SAM181b638b8248
#> SAM3b15b4c6311d
#> SAM1c8b086175ca
#> SAMeaa477a5384b
#> SAM1dda30f1c5be
#> SAM9daccafc18db
#> SAMbc8dc3a7b54e
#> SAM7829a341b9f3
#> SAM0f956e757453
#> SAMb2e4a082541a
#> SAM771445e92421
#> SAM59fda9035d1d
#> SAM6f2a102a99df
#> SAMe56c96c51190
#> SAM6662f5181f87
#> SAM18b9351e265a
#> SAM7fb7a13c096b
#> SAM87a8e18eb45b
#> SAM9725303dce0c
#> SAM6ff654a20f98
#> SAM94859b440b1d
#> SAMaabf4afe4213
#> SAMce39dd79b441
#> SAM6964a6d7b967
#> SAMf20b827dca51
#> SAM31d9176e11fb
#> SAM80c6183220e6
#> SAM1f83ebd6be9b
#> SAMe7e4f7c076a7
#> SAMc6eff056c89a
#> SAM5cfa1699bdb7
#> SAMda4d892fddc8
#> SAMe3d4266775a9 
#> 
#> Total sum across samples of combination Epidish_BPRNACan_ is 1
#> Total sum across samples of combination Epidish_BPRNACanProMet is 1
#> Total sum across samples of combination Epidish_BPRNACan3DProMet is 1
#> Total sum across samples of combination Epidish_CBSX.HNSCC.scRNAseq is 1
#> Total sum across samples of combination Epidish_CBSX.Melanoma.scRNAseq is 1
#> Total sum across samples of combination Epidish_CBSX.NSCLC.PBMCs.scRNAseq is 1
#> 
#> Total sum across samples of combination Epidish_CCLE.TIL10 is not 1! Remember these are proportions and the total should be 1
#> Samples which sum with combination Epidish_CCLE.TIL10 is not 1:
#> 
#>  SAM7f0d9cc7f001
#> SAMcf018fee2acd
#> SAMe7bf6c015192
#> SAM6dd7ad1d797d
#> SAMc692536a795a
#> SAM9a2cf3c06fb3
#> SAM557dde1b9f3e
#> SAM23aa15d4a0b0
#> SAM7fb6987514a4
#> SAMd1bd63734394
#> SAMbe83eae4026e
#> SAMe5bc41772bc9
#> SAM7114d99032ec
#> SAMbf1a3ae828e6
#> SAM52e3fa3ad574
#> SAM9cafb905b36a
#> SAM032c642382a7
#> SAM0ce9c983b20f
#> SAM2f228939632f
#> SAM1fa6bcb7fc48
#> SAM110501d0eedb
#> SAM2e9ac0b1b250
#> SAMd636e3461955
#> SAMd98bac0a070f
#> SAM8e8ef2368dfa
#> SAMb4c7a001537d
#> SAMfd947610629d
#> SAM943df5cf15df
#> SAM39eb94fa504d
#> SAM62fb1388c871
#> SAMae4da274eded
#> SAM75142fcab9df
#> SAM166a419a4e5a
#> SAM025b45c27e05
#> SAM560f23d6a3ad
#> SAMf3a9bce50099
#> SAM34430ef08e5b
#> SAM4b0175e8db6e
#> SAM5e3bae090b8c
#> SAM3cb94b0d5297
#> SAMe97af0feefdf
#> SAMa424c75831b4
#> SAM73663ee4a96e
#> SAM3330c03fdf00
#> SAM99a46b9eec27
#> SAM553c3c35b847
#> SAMc2a1820d4e6b
#> SAMe1eb5d988760
#> SAMb15ad09d6e24
#> SAM961d04c42bd9
#> SAMb0d11db9aa79
#> SAM19fec8f3b3bd
#> SAM415f36ad349e
#> SAMad83c9c53537
#> SAMd135d5867fe3
#> SAM065890737112
#> SAMb470eb8f04be
#> SAM4b7ea015fd9e
#> SAM18be5b395318
#> SAM0d855cff64e6
#> SAM6cbc10abddb0
#> SAM203dcf14f927
#> SAM76a431ba6ce1
#> SAMd3bd67996035
#> SAMd3601288319e
#> SAMba1a34b5a060
#> SAMfed609955db9
#> SAM5c139c5c1c4f
#> SAM49d48750e294
#> SAM7aa01fc49a80
#> SAM1c0ecfb3eb63
#> SAM9aa6a095a9d6
#> SAM957378bd907f
#> SAM670649e105b5
#> SAM5fe7a81a39dd
#> SAM181b638b8248
#> SAM3b15b4c6311d
#> SAM1dda30f1c5be
#> SAM7829a341b9f3
#> SAMb2e4a082541a
#> SAM771445e92421
#> SAMe56c96c51190
#> SAM18b9351e265a
#> SAM9725303dce0c
#> SAM6ff654a20f98
#> SAMce39dd79b441
#> SAMf20b827dca51
#> SAM31d9176e11fb
#> SAM80c6183220e6
#> SAM572f19794c96
#> SAMe7e4f7c076a7
#> SAMc6eff056c89a
#> SAMda4d892fddc8 
#> 
#> 
#> Total sum across samples of combination Epidish_TIL10 is not 1! Remember these are proportions and the total should be 1
#> Samples which sum with combination Epidish_TIL10 is not 1:
#> 
#>  SAM7f0d9cc7f001
#> SAMcf018fee2acd
#> SAMdf3e42c8672a
#> SAMe7bf6c015192
#> SAM6dd7ad1d797d
#> SAMc692536a795a
#> SAM557dde1b9f3e
#> SAM23aa15d4a0b0
#> SAMb963dda93cfd
#> SAM7fb6987514a4
#> SAM63405b04ab2d
#> SAM18bc1078bc15
#> SAMd1bd63734394
#> SAMbe83eae4026e
#> SAM23095936e611
#> SAM7114d99032ec
#> SAMbf1a3ae828e6
#> SAM52e3fa3ad574
#> SAMd4c0837b0997
#> SAM9cafb905b36a
#> SAM032c642382a7
#> SAM0ce9c983b20f
#> SAM2f228939632f
#> SAM297c0301e861
#> SAM1fa6bcb7fc48
#> SAM110501d0eedb
#> SAM2e9ac0b1b250
#> SAMd636e3461955
#> SAMb4c7a001537d
#> SAM39eb94fa504d
#> SAM62fb1388c871
#> SAM5b57e47fdcb3
#> SAMae4da274eded
#> SAMd35318127278
#> SAM75142fcab9df
#> SAM560f23d6a3ad
#> SAMf3a9bce50099
#> SAM34430ef08e5b
#> SAM5e3bae090b8c
#> SAM3cb94b0d5297
#> SAMcb132b0cdd2c
#> SAMe97af0feefdf
#> SAMdad5c29dc105
#> SAMa424c75831b4
#> SAM45c8e6412c66
#> SAM73663ee4a96e
#> SAM99a46b9eec27
#> SAM553c3c35b847
#> SAMc2a1820d4e6b
#> SAMe1eb5d988760
#> SAMe3210d3632b4
#> SAMb15ad09d6e24
#> SAM961d04c42bd9
#> SAMb0d11db9aa79
#> SAM415f36ad349e
#> SAMad83c9c53537
#> SAM3ee5dcd894f0
#> SAM8a1b0e02ee42
#> SAMd7d57ee3a863
#> SAMeb29625f76a5
#> SAMd135d5867fe3
#> SAM065890737112
#> SAMb470eb8f04be
#> SAM4b7ea015fd9e
#> SAM9306c5c92444
#> SAM18be5b395318
#> SAM6cbc10abddb0
#> SAMdee1011782cd
#> SAM203dcf14f927
#> SAMe41b1e773582
#> SAM76a431ba6ce1
#> SAMd3bd67996035
#> SAMd3601288319e
#> SAMba1a34b5a060
#> SAM18a4dabbc557
#> SAMfed609955db9
#> SAM5c139c5c1c4f
#> SAM6780ed436b55
#> SAM85f0a3ac1c45
#> SAM7aa01fc49a80
#> SAM1c0ecfb3eb63
#> SAMbcb07ba81cee
#> SAM9aa6a095a9d6
#> SAM670649e105b5
#> SAM5fe7a81a39dd
#> SAM29da928587ad
#> SAM181b638b8248
#> SAM3b15b4c6311d
#> SAMeaa477a5384b
#> SAM1dda30f1c5be
#> SAM7829a341b9f3
#> SAMb2e4a082541a
#> SAM771445e92421
#> SAMe56c96c51190
#> SAM18b9351e265a
#> SAM87a8e18eb45b
#> SAM9725303dce0c
#> SAMce39dd79b441
#> SAMf20b827dca51
#> SAM31d9176e11fb
#> SAM80c6183220e6
#> SAM572f19794c96
#> SAMe7e4f7c076a7
#> SAMc6eff056c89a
#> SAMda4d892fddc8 
#> 
#> 
#> Total sum across samples of combination Epidish_LM22 is not 1! Remember these are proportions and the total should be 1
#> Samples which sum with combination Epidish_LM22 is not 1:
#> 
#>  SAMd35318127278
#> SAM560f23d6a3ad
#> SAM73663ee4a96e
#> SAM5234688806a7
#> SAMae02629a97f7
#> SAM87a8e18eb45b 
#> 
#> Total sum across samples of combination Epidish_CBSX.Melanoma.scRNAseq is 1
#> Warning in compute.deconvolution.preprocessing(data.frame(all_deconvolution_table)): 
#> Please verify your matrix
```

``` r
head(deconv[,1:5])
#>                 Quantiseq_TIL10_B.cells Epidish_BPRNACan_B.cells
#> SAM7f0d9cc7f001             0.050210988             3.563950e-02
#> SAM4305ab968b90             0.008603929             1.073465e-04
#> SAMcf018fee2acd             0.042857608             1.708482e-02
#> SAMcc4675f394a1             0.025216439             9.074189e-05
#> SAM49f9b2e57aa5             0.022668958             1.900745e-02
#> SAM2e7aa8fa0ab3             0.012022510             4.130289e-03
#>                 Epidish_BPRNACan3DProMet_B.cells Epidish_BPRNACanProMet_B.cells
#> SAM7f0d9cc7f001                      0.051634688                    0.061927146
#> SAM4305ab968b90                      0.007359397                    0.013455460
#> SAMcf018fee2acd                      0.024649003                    0.029724100
#> SAMcc4675f394a1                      0.000000000                    0.002167176
#> SAM49f9b2e57aa5                      0.025668504                    0.032688442
#> SAM2e7aa8fa0ab3                      0.008716558                    0.011090919
#>                 Epidish_CBSX.HNSCC.scRNAseq_B.cells
#> SAM7f0d9cc7f001                         0.025626151
#> SAM4305ab968b90                         0.000000000
#> SAMcf018fee2acd                         0.080914523
#> SAMcc4675f394a1                         0.001434758
#> SAM49f9b2e57aa5                         0.137032856
#> SAM2e7aa8fa0ab3                         0.282602048
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

To evaluate the association between cell group scores and a clinical
trait, you can perform statistical tests such as Fisher’s exact test,
ANOVA, or any other supported method (Wilcoxon, Kruskal–Wallis, t-test)
using the unified function
[`cell.groups.stat.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.stat.analysis.md).

``` r
res_fisher <- cell.groups.stat.analysis(cell.groups = cell_groups, coldata = traitdata, 
                                        trait = "Best.Confirmed.Overall.Response", method = "fisher", pval = 0.05)

res_anova <- cell.groups.stat.analysis(cell.groups = cell_groups, coldata = traitdata,
                                       trait = "Best.Confirmed.Overall.Response", method = "anova", pval = 0.05)

res_wilcox <- cell.groups.stat.analysis(cell.groups = cell_groups, coldata = traitdata, 
                                        trait = "Best.Confirmed.Overall.Response", method = "wilcox", pval = 0.05)

res_krustal <- cell.groups.stat.analysis(cell.groups = cell_groups, coldata = traitdata,
                                         trait = "Best.Confirmed.Overall.Response", method = "krustal", pval = 0.05)
```

### **One-step CellTFusion**

For convenience, we offer an all-in-one function that automates the
entire CellTFusion workflow—calculating features, performing
intermediate analyses, and ultimately returning the cell group scores.
To run this streamlined pipeline, simply execute:

``` r
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
[`prepare_CellTFusion_folds()`](https://verapancaldilab.github.io/CellTFusion/reference/prepare_CellTFusion_folds.md),
which is fully compatible with the `pipeML` package for machine learning
model training and evaluation.

The
[`prepare_CellTFusion_folds()`](https://verapancaldilab.github.io/CellTFusion/reference/prepare_CellTFusion_folds.md)
function handles cross-validation by recomputing `CellTFusion` features
inside each fold. This ensures that:

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

Alvarez, Mariano J., Yao Shen, Federico M. Giorgi, Alexander Lachmann,
B. Belinda Ding, B. Hilda Ye, and Andrea Califano. 2016. “Functional
Characterization of Somatic Mutations in Cancer Using Network-Based
Inference of Protein Activity.” *Nature Genetics* 48 (8): 838–47.
<https://doi.org/10.1038/ng.3593>.

Badia-i-Mompel, Pau, Jesús Vélez Santiago, Jana Braunger, Celina Geiss,
Daniel Dimitrov, Sophia Müller-Dott, Petr Taus, et al. 2022. “decoupleR:
Ensemble of Computational Methods to Infer Biological Activities from
Omics Data.” *Bioinformatics Advances* 2 (1): vbac016.
<https://doi.org/10.1093/bioadv/vbac016>.

Langfelder, Peter, and Steve Horvath. 2008. “WGCNA: An r Package for
Weighted Correlation Network Analysis.” *BMC Bioinformatics* 9 (1): 559.
<https://doi.org/10.1186/1471-2105-9-559>.

Müller-Dott, Sophia, Eirini Tsirvouli, Miguel Vazquez, Ricardo O Ramirez
Flores, Pau Badia-i-Mompel, Robin Fallegger, Dénes Türei, Astrid
Lægreid, and Julio Saez-Rodriguez. 2023. “Expanding the Coverage of
Regulons from High-Confidence Prior Knowledge for Accurate Estimation of
Transcription Factor Activities.” *Nucleic Acids Research* 51 (20):
10934–49. <https://doi.org/10.1093/nar/gkad841>.

Schubert, Michael, Bertram Klinger, Martina Klünemann, Anja Sieber,
Florian Uhlitz, Sascha Sauer, Mathew J. Garnett, Nils Blüthgen, and
Julio Saez-Rodriguez. 2018. “Perturbation-Response Genes Reveal
Signaling Footprints in Cancer Gene Expression.” *Nature Communications*
9 (1): 20. <https://doi.org/10.1038/s41467-017-02391-6>.
