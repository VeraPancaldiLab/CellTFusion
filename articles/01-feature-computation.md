# Feature Computation

``` r

library(CellTFusion)
#> 
#> 
```

CellTFusion requires two primary inputs: **cell-type deconvolution
proportions** and **transcription factor (TF) activity scores**. This
article walks through how to compute both, followed by TF module
construction and pathway activity scoring.

Load the pre-packaged example data:

``` r

raw.counts  <- CellTFusion::raw.counts.tuto
traitdata   <- CellTFusion::traitdata.tuto
```

### Cell-type deconvolution

To compute cell-type proportions, we use the `multideconv` R package,
which integrates multiple deconvolution algorithms and signature
matrices to estimate cell-type abundances from bulk RNA-seq data. For
more details, visit the [multideconv GitHub
repository](https://github.com/VeraPancaldiLab/multideconv).

``` r

deconv <- multideconv::compute.deconvolution(
  raw.counts,
  methods    = c("Quantiseq", "Epidish"),
  normalized = TRUE,
  return     = FALSE
)
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
head(deconv[, 1:5])
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

### TF activity inference

To infer transcription factor activity, we use the `viper` package
([Alvarez et al. 2016](#ref-Alvarez2016)) combined with the `CollecTRI`
regulon database ([Müller-Dott et al. 2023](#ref-10.1093/nar/gkad841)).
This approach estimates TF activity from the expression levels of
downstream target genes.

``` r

# Normalize by log2(TPM + 1)
counts.norm <- data.frame(ADImpute::NormalizeTPM(raw.counts, log = TRUE))

universe <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)
tfs <- compute.TFs.activity(counts.norm, universe = universe)
head(tfs[, 1:5])
```

### TF module construction (WTCNA)

Before constructing cell groups, we reduce the TF activity matrix into
modules of TFs that exhibit similar activity patterns across samples.
This is done using a Weighted TF Correlation Network Analysis (WTCNA),
adapted from WGCNA ([Langfelder and Horvath
2008](#ref-langfelder2008wgcna)).

``` r

network <- compute.WTCNA(tfs, corr_mod = 0.8, clustering.method = "ward.D2", return = TRUE)
```

To explore how clinical variables associate with TF modules, use
[`compute.metadata.association()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.metadata.association.md).
Plots are saved in the `Results/` directory:

``` r

compute.metadata.association(network[[1]], traitdata, pval = 0.05, file.name = "Tutorial", width = 10)
```

### Pathway activity scoring

To functionally characterize TF modules, we estimate pathway activities
using a multivariate linear model (MLM) via `decoupleR` ([Badia-i-Mompel
et al. 2022](#ref-10.1093/bioadv/vbac016)) with the `PROGENy` database
([Schubert et al. 2018](#ref-Schubert2018)) by default. A custom gene
set can also be supplied, in which case GSVA is applied instead.

``` r

pathways <- compute.pathway.activity(counts.norm)
```

To visualize the relationship between TF modules and pathways, use
[`compute.modules.relationship()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md).
This saves a labeled correlation heatmap to `Results/`:

``` r

compute.modules.relationship(network[[1]], pathways, "Pathways_Progeny-TFs_Modules", width = 15)
```

### Hub TF enrichment analysis

If pathway scoring alone is insufficient, you can perform
Over-Representation Analysis (ORA) on the Reactome database using hub
TFs. The function identifies hub TFs per module (based on module
membership and node degree) and runs enrichment on their target genes:

``` r

hub_tfs <- identify_hub_TFs(t(tfs), network, MM_thresh = 0.8, degree_thresh = 0.9)
compute.modules.enrichment(counts.norm, hub_tfs)
```

### Deconvolution dimensionality reduction

To reduce dimensionality in the deconvolution output and improve
statistical power downstream, we group cell types with similar abundance
patterns using `multideconv`:

``` r

dt <- multideconv::compute.deconvolution.analysis(deconv, corr = 0.7, seed = 123)
dt <- multideconv::deconvolution_dictionary(dt, pathways)

# Correlate deconvolution groups with TF modules
compute.modules.relationship(
  network[[1]], dt[[1]], "Deconvolution-TFs_Modules",
  vertical = TRUE, height = 30, width = 10, pval = 0.05
)
```

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
