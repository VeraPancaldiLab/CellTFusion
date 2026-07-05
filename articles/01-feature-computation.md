# Feature Computation

``` r

library(CellTFusion)
#> 
#> 
```

CellTFusion combines two families of features computed from bulk RNAseq
data: **cell-type deconvolution proportions** and **transcription factor
(TF) activity scores**, further summarized into TF co-activity modules
and pathway activities. This article walks through each function
individually — what it does, its key arguments, and the plots/objects it
produces. All of these steps are run automatically, in the right order,
by the
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
wrapper (see the package
[README](https://verapancaldilab.github.io/CellTFusion/index.md)); read
this article if you want to run, inspect, or customize each stage on its
own.

Load the pre-packaged example data:

``` r

raw.counts <- CellTFusion::raw.counts.tuto
traitdata  <- CellTFusion::traitdata.tuto
```

## 1. Normalization

Before TF activity or pathway scoring, raw counts are normalized to
log2(TPM + 1). Deconvolution methods use TPM without the log transform,
and this step is handled internally by `compute.deconvolution()`.

``` r

counts.norm <- data.frame(ADImpute::NormalizeTPM(raw.counts, log = TRUE))
```

## 2. Cell-type deconvolution

Cell-type proportions are estimated with `compute.deconvolution()` from
the companion
[`multideconv`](https://github.com/VeraPancaldiLab/multideconv) package,
which wraps multiple deconvolution algorithms (Quantiseq, Epidish,
DeconRNASeq, DWLS, CIBERSORTx) and averages/combines their outputs into
a single feature matrix (samples x cell-type/method combinations).

Key arguments:

- `raw.counts` — raw counts matrix (genes x samples).
- `methods` — deconvolution algorithms to run.
- `normalized` — set to `TRUE` if `raw.counts` is already normalized
  (the function otherwise normalizes internally to TPM).
- `credentials.mail` / `credentials.token` — required only if `"CBSX"`
  (CIBERSORTx) is included in `methods`.
- `return` — if `TRUE`, saves the resulting matrix to `Results/`.

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

## 3. TF activity inference

TF activity per sample is inferred with
[`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md)
from normalized expression using a VIPER/consensus-scoring approach (via
`decoupleR`) ([Alvarez et al. 2016](#ref-Alvarez2016)), combined with a
regulon (TF-target network).

Key arguments:

- `RNA.counts` — normalized expression matrix (genes x samples),
  typically log2(TPM+1).
- `TF.collection` — source of the TF-target regulon:
  - `"CollecTRI"` (default) ([Müller-Dott et al.
    2023](#ref-10.1093/nar/gkad841)) and `"Dorothea"` are curated,
    literature-derived regulons fetched automatically from OmnipathR.
  - `"ARACNE"` uses a **data-driven, cohort-specific** regulatory
    network reverse-engineered from expression data with the ARACNe
    algorithm ([Margolin et al. 2006](#ref-Margolin2006)) (Algorithm for
    the Reconstruction of Accurate Cellular Networks), instead of a
    generic literature-derived network. ARACNe infers TF-target edges
    from mutual information between gene expression profiles across the
    cohort, followed by removal of indirect interactions via the Data
    Processing Inequality. This is useful when you have a large cohort
    of samples from the same cancer type/context and want TF-target
    relationships inferred directly from its co-expression structure
    rather than relying on curated interactions that may not hold in
    that tissue — at the cost of requiring a sizeable, homogeneous
    cohort to estimate mutual information reliably. It requires a
    pre-computed 3-column network file (`regulator`, `target`,
    `mutual information`) located at
    `input/ARACNE/<cancer.type>/network/network.txt`; pass the
    corresponding label via `cancer.type`.
- `min_targets_size` — minimum number of target genes required per
  regulon (TFs with fewer targets are dropped). Default is 5; increasing
  it keeps only better-supported TFs.
- `universe` — optional user-supplied regulon (data frame of TF-target
  interactions), bypassing the automatic OmnipathR/ARACNE fetch.
- `cores` — number of cores used for the activity inference.
- `scale` — if `TRUE` (default), z-scores TF activity across samples.

``` r

tfs <- compute.TFs.activity(
  RNA.counts       = counts.norm,
  TF.collection    = "CollecTRI",
  min_targets_size = 5,
  cores            = 3,
  return           = TRUE
)
head(tfs[, 1:5])
```

You can also supply a pre-built regulon directly via `universe`:

``` r

universe <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)
tfs <- compute.TFs.activity(counts.norm, universe = universe)
```

Or use a cohort-specific ARACNe network:

``` r

tfs_aracne <- compute.TFs.activity(
  RNA.counts    = counts.norm,
  TF.collection = "ARACNE",
  cancer.type   = "skcm"   # matches input/ARACNE/skcm/network/network.txt
)
```

## 4. TF co-activity modules

TF activity scores (potentially hundreds of TFs) are reduced into a
small number of co-activity **modules** using
[`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md),
an implementation of Weighted TF Co-activity Network Analysis (WTCNA),
an adaptation of WGCNA ([Langfelder and Horvath
2008](#ref-langfelder2008wgcna)) applied to TFs instead of genes. Each
module gets one sample-level score (its eigengene).

Key arguments:

- `TFs.matrix` — TF activity matrix (or, if `batch = TRUE`, a **list**
  of per-cohort matrices) from
  [`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md).
- `batch` — if `TRUE`, runs consensus WGCNA across cohorts instead of a
  single-matrix analysis (see the [Batch/Multi-cohort
  Analysis](https://verapancaldilab.github.io/CellTFusion/articles/07-batch-analysis.md)
  article for details).
- `network.type` — `"signed"` (default), `"unsigned"`,
  `"signed hybrid"`, or `"distance"`.
- `minMod` — minimum number of TFs required to form a module.
- `corr_mod` — correlation threshold for merging similar modules.
- `cor_type` — `"p"` (Pearson) or `"s"` (Spearman).

``` r

network <- compute.WTCNA(
  TFs.matrix        = tfs,
  batch             = FALSE,
  network.type      = "signed",
  clustering.method = "ward.D2",
  minMod            = 15,
  corr_mod          = 0.9,
  return            = TRUE
)
```

When `return = TRUE`, this saves diagnostic plots to `Results/`. First,
the scale-free topology fit used to pick the WGCNA soft-thresholding
power:

![Scale-free topology model fit (signed R^2) versus soft-thresholding
power for the WGCNA network](figures/soft_threshold.png)

Then the TF clustering dendrogram, colored by module, before and after
merging highly correlated modules:

![TF clustering dendrogram with module colors before
merging](figures/gene_dendrogram_before_merging.png)![TF clustering
dendrogram with module colors after
merging](figures/gene_dendrogram_after_merging.png)

`network[[1]]` (`"TFs module matrix"`) holds the module eigengene scores
(samples x modules) used downstream by
[`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md).

## 5. Pathway activity

Pathway activities are estimated with
[`compute.pathway.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md),
using a multivariate linear model (MLM) via `decoupleR` ([Badia-i-Mompel
et al. 2022](#ref-10.1093/bioadv/vbac016)), with PROGENy ([Schubert et
al. 2018](#ref-Schubert2018)) as the default pathway database. If a
custom `gene_sets` list is supplied, GSVA scoring is computed in
addition to PROGENy.

Key arguments:

- `RNA.tpm` — normalized expression matrix (genes x samples).
- `gene_sets` — optional list of custom gene sets; if provided, GSVA
  scores are also returned.
- `paths` — optional custom PROGENy-style pathway-gene weight table; if
  `NULL`, the default human PROGENy model is used.

``` r

pathways <- compute.pathway.activity(
  RNA.tpm = counts.norm,
  return  = TRUE
)
```

The result is a samples x pathways score matrix (14 PROGENy pathways by
default). Clustering samples and pathways on this matrix gives a quick
overview of which pathways co-vary across the cohort:

![Heatmap of PROGENy pathway activity scores clustered by sample and
pathway](figures/pathway_heatmap.png)

## 6. Deconvolution feature reduction

Deconvolution features from different methods/signatures are often
highly correlated. `compute.deconvolution.analysis()` (from
`multideconv`) groups correlated cell-type/method columns into
**subgroups**, reducing redundancy and improving statistical power in
later steps.

Key arguments:

- `deconvolution` — output of `compute.deconvolution()`.
- `corr` — minimum correlation to group features into the same subgroup.
- `corr_type` — `"spearman"` (default) or `"pearson"`.
- `batch` — optional batch/cohort vector; if provided, correlations are
  computed as partial correlations controlling for batch (see
  [Batch/Multi-cohort
  Analysis](https://verapancaldilab.github.io/CellTFusion/articles/07-batch-analysis.md)).

``` r

dt <- multideconv::compute.deconvolution.analysis(
  deconvolution = deconv,
  corr          = 0.7,
  seed          = 123,
  return        = FALSE
)
```

## 7. Cell group construction

[`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
combines the TF module network (`network`) and reduced deconvolution
structure (`dt`) into composite **cell groups**: clusters of cell-type
features whose abundance correlates with specific TF module activity. It
correlates cell-type subgroups against TF modules, hierarchically
clusters the significant pairs, and cuts the resulting dendrogram into
cell groups.

``` r

cell_groups <- construct_cell_groups(
  network = network,
  dt      = dt,
  pval    = 0.05
)
```

One dendrogram is produced per TF module, showing how its correlated
cell-type features cluster into cell groups (colors indicate module
membership of the underlying TF driving each branch):

![Dendrogram of cell-type features clustered into cell groups per TF
module](figures/cell_groups_dendrogram.png)

See the [Cell Group
Construction](https://verapancaldilab.github.io/CellTFusion/articles/02-cell-groups.md)
article for full details, including the difference between unsupervised
and supervised construction.

## 8. Latent factor extraction

Cell group scores are further compressed into a small set of
non-negative matrix factorization (NMF) **latent factors**, a compact
representation of the TME landscape used for downstream statistics and
machine learning, via
[`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md).
See the [Cell Group
Construction](https://verapancaldilab.github.io/CellTFusion/articles/02-cell-groups.md)
article for a full walkthrough.

## 9. Cell niche characterization

[`compute_cells_niches()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_cells_niches.md)
characterizes each latent factor by identifying which cell types are
enriched among the cell groups that contribute most to it. See the [Cell
Group
Construction](https://verapancaldilab.github.io/CellTFusion/articles/02-cell-groups.md)
article for a full walkthrough.

## 10. Functional and meta-program annotation

[`compute_factor_gsea()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_factor_gsea.md)
and
[`map_factors_to_metaprograms()`](https://verapancaldilab.github.io/CellTFusion/reference/map_factors_to_metaprograms.md)
functionally annotate latent factors via Hallmark GSEA and map them onto
reference cancer meta-programs derived from TCGA. See the [TME State
Characterisation](https://verapancaldilab.github.io/CellTFusion/articles/03-tme-states.md)
article for a detailed walkthrough, including the biological rationale,
references, and example plots behind meta-program mapping.

## Putting it together

Every function above is called, in this order, by
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md).
Once you are comfortable with the individual building blocks, you
normally do not need to run them manually — simply call:

``` r

res <- CellTFusion(
  raw.counts     = raw.counts,
  normalized     = TRUE,
  coldata        = traitdata,
  task           = "unsupervised",
  deconv_methods = c("Quantiseq", "Epidish"),
  TF.collection  = "CollecTRI",
  cancer_type    = "skcm",
  corr           = 0.7,
  corr_mod       = 0.9,
  pval           = 0.05,
  file_name      = "Tutorial",
  return         = TRUE
)
```

`res` then contains every intermediate object described above
(`$Deconvolution`, `$TFs_matrix`, `$TF_network`, `$Pathways_scores`,
`$Processed_deconvolution`, `$Cell_groups`, `$Latent_spaces`,
`$Cells_niches`, `$TME_states`, `$Metaprograms_reference`) — see the
[README](https://verapancaldilab.github.io/CellTFusion/index.md) for the
full output structure.

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

Margolin, Adam A., Ilya Nemenman, Katia Basso, et al. 2006. “ARACNE: An
Algorithm for the Reconstruction of Gene Regulatory Networks in a
Mammalian Cellular Context.” *BMC Bioinformatics* 7 (Suppl 1): S7.
<https://doi.org/10.1186/1471-2105-7-S1-S7>.

Müller-Dott, Sophia, Eirini Tsirvouli, Miguel Vazquez, et al. 2023.
“Expanding the Coverage of Regulons from High-Confidence Prior Knowledge
for Accurate Estimation of Transcription Factor Activities.” *Nucleic
Acids Research* 51 (20): 10934–49.
<https://doi.org/10.1093/nar/gkad841>.

Schubert, Michael, Bertram Klinger, Martina Klünemann, et al. 2018.
“Perturbation-Response Genes Reveal Signaling Footprints in Cancer Gene
Expression.” *Nature Communications* 9 (1): 20.
<https://doi.org/10.1038/s41467-017-02391-6>.
