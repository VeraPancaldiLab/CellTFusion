# Compute modules relationship

Correlates every column of one matrix with every column of another
(e.g., TF modules and external features such as deconvolution estimates
or pathway activities) across the same set of samples, and visualizes
the significant associations. It is meant for continuous variables only:
all the columns of both matrices must be numeric. For categorical
variables (e.g. clinical groups), use
[`compute.metadata.association()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.metadata.association.md).

## Usage

``` r
compute.modules.relationship(
  matA,
  matB,
  file_name = NULL,
  batch = NULL,
  pval = 0.05,
  padj = FALSE,
  cor_type = "p",
  vertical = F,
  plot = T,
  plot.grid = F,
  ncol.grid = NULL
)
```

## Arguments

- matA:

  A numeric matrix or data frame of continuous features (samples x
  features), e.g. TF module scores.

- matB:

  A numeric matrix or data frame of continuous features to correlate
  with (samples x features).

- file_name:

  A string indicating the base name (without extension) of the figures
  to be saved in the "Results/" folder. If `NULL` (default), no figure
  is saved.

- batch:

  Optional vector indicating batch assignment for samples (no missing
  values). If it has two or more levels, partial correlations
  ([`ppcor::pcor.test()`](https://rdrr.io/pkg/ppcor/man/pcor.test.html),
  Pearson or Spearman according to `cor_type`) are computed controlling
  for batch. Batch is treated as categorical (one dummy variable per
  level beyond the first).

- pval:

  A numeric threshold for the p-value to define significance (p-value,
  or adjusted p-value with `padj = TRUE`, \<= `pval`). Default is 0.05.

- padj:

  Logical; if TRUE, the p-values are corrected for multiple testing with
  the Benjamini-Hochberg method, across the columns of `matB` for each
  column of `matA`, and the adjusted p-values are used to define
  significance (in `significant` and in the figures). Default is FALSE.

- cor_type:

  Type of correlation to compute: "p" (Pearson) or "s" (Spearman).
  Default is "p". It also sets the type of the partial correlations
  computed with `batch`.

- vertical:

  Logical; if TRUE, modules are on the x-axis and the `matB` features on
  the y-axis. Otherwise (default), the `matB` features are on the x-axis
  and modules on the y-axis.

- plot:

  Logical; if TRUE, saves the heatmap plot as a PDF (`file_name` is
  needed). Default is TRUE.

- plot.grid:

  Logical; if TRUE, generates per-pair scatter grid plots for
  significant associations, saved to
  "Results/\<file_name\>\_scatter_grid.svg" (the 100 strongest pairs at
  most). With Spearman, the ranks of the values are shown, as the
  correlation is computed on them. With `batch`, the values shown are
  corrected for batch, as the correlations are. Default is FALSE.

- ncol.grid:

  Integer number of columns used in scatter grid layout. If NULL
  (default), it is set to the ceiling of the square root of the number
  of panels.

## Value

Invisibly, a list with:

- `correlation`: the correlation matrix (columns of `matA` as rows,
  columns of `matB` as columns).

- `significant`: a named list with, for each column of `matA`, the
  columns of `matB` significantly correlated with it (p-value, or
  adjusted p-value with `padj = TRUE`, \<= `pval`).

- `p_value`: the matrix of p-values.

- `p_adjusted`: the matrix of Benjamini-Hochberg adjusted p-values
  (equal to `p_value` with `padj = FALSE`).

If `plot = TRUE` and `file_name` is given, a heatmap of the correlations
is saved to "Results/\<file_name\>.pdf" (the ".pdf" extension is added
if missing). Its size is computed from the number of features. With up
to 40 features per axis, the significant cells show the correlation and
its significance (\* \< 0.05, \*\* \< 0.01, \*\*\* \< 0.001, \*\*\*\* \<
0.0001); with more, they are marked with an asterisk.

## Details

`matA` and `matB` must have the same samples in the same order (same
rownames); names are compared after
[`make.names()`](https://rdrr.io/r/base/make.names.html), so e.g.
"TCGA-XX" and "TCGA.XX" match. Correlations and p-values use, for each
pair of columns, the samples with a value in both (missing values are
allowed); a pair with fewer than three such samples, or without
variation among them, gets `NA`. Columns without variation are skipped
with a message. Column names must be unique within each matrix.

## Examples

``` r
data("network.tuto")
data("deconv_subgroups.tuto")
corr = compute.modules.relationship(network.tuto[[1]],
                                    deconv_subgroups.tuto[[1]],
                                    plot = FALSE,
                                    pval = 0.01)
corr$correlation[1:3, 1:3]
#>        Quantiseq_B.cells DeconRNASeq_CBSX.Melanoma.scRNAseq_B.cells
#> brown         0.02792722                                -0.04395709
#> purple       -0.06263770                                -0.03061520
#> blue         -0.06526158                                -0.05794571
#>        DeconRNASeq_CBSX.NSCLC.PBMCs.scRNAseq_B.cells
#> brown                                     0.05627460
#> purple                                   -0.06805893
#> blue                                      0.04245497
corr$significant[[1]]
#>  [1] "B.cells_Subgroup.4"                                  
#>  [2] "DeconRNASeq_CBSX.Melanoma.scRNAseq_Macrophages.cells"
#>  [3] "Macrophages.cells_Subgroup.1"                        
#>  [4] "Macrophages.M0_Subgroup.1"                           
#>  [5] "Macrophages.M0_Subgroup.2"                           
#>  [6] "Quantiseq_Macrophages.M1"                            
#>  [7] "Macrophages.M1_Subgroup.1"                           
#>  [8] "Macrophages.M1_Subgroup.2"                           
#>  [9] "Macrophages.M1_Subgroup.4"                           
#> [10] "Macrophages.M1_Subgroup.5"                           
#> [11] "Quantiseq_Macrophages.M2"                            
#> [12] "CBSX_CCLE.TIL10_Macrophages.M2"                      
#> [13] "DeconRNASeq_LM22_Macrophages.M2"                     
#> [14] "Macrophages.M2_Subgroup.2"                           
#> [15] "Macrophages.M2_Subgroup.3"                           
#> [16] "Macrophages.M2_Subgroup.4"                           
#> [17] "Macrophages.M2_Subgroup.5"                           
#> [18] "Monocytes_Subgroup.1"                                
#> [19] "Monocytes_Subgroup.4"                                
#> [20] "CBSX_BPRNACanProMet_Neutrophils"                     
#> [21] "CBSX_TIL10_Neutrophils"                              
#> [22] "Neutrophils_Subgroup.1"                              
#> [23] "Neutrophils_Subgroup.7"                              
#> [24] "Neutrophils_Subgroup.8"                              
#> [25] "Quantiseq_NK.cells"                                  
#> [26] "DeconRNASeq_CBSX.Melanoma.scRNAseq_NK.cells"         
#> [27] "CBSX_TIL10_NK.cells"                                 
#> [28] "NK.cells_Subgroup.1"                                 
#> [29] "DeconRNASeq_LM22_NK.resting"                         
#> [30] "DeconRNASeq_CBSX.NSCLC.PBMCs.scRNAseq_NKT.cells"     
#> [31] "Epidish_CBSX.NSCLC.PBMCs.scRNAseq_NKT.cells"         
#> [32] "DeconRNASeq_BPRNACan_CD4.cells"                      
#> [33] "DeconRNASeq_CBSX.Melanoma.scRNAseq_CD4.cells"        
#> [34] "DeconRNASeq_CBSX.NSCLC.PBMCs.scRNAseq_CD4.cells"     
#> [35] "CD4.memory.resting_Subgroup.1"                       
#> [36] "CD4.naive_Subgroup.1"                                
#> [37] "DeconRNASeq_BPRNACan_CD8.cells"                      
#> [38] "DeconRNASeq_CBSX.HNSCC.scRNAseq_CD8.cells"           
#> [39] "DeconRNASeq_CBSX.Melanoma.scRNAseq_CD8.cells"        
#> [40] "DeconRNASeq_CBSX.NSCLC.PBMCs.scRNAseq_CD8.cells"     
#> [41] "DeconRNASeq_LM22_T.cells.helper"                     
#> [42] "DeconRNASeq_CCLE.TIL10_Cancer"                       
#> [43] "DeconRNASeq_CBSX.Melanoma.scRNAseq_Cancer"           
#> [44] "Epidish_CBSX.Melanoma.scRNAseq_Cancer"               
#> [45] "Cancer_Subgroup.1"                                   
#> [46] "Cancer_Subgroup.3"                                   
#> [47] "Cancer_Subgroup.4"                                   
#> [48] "DeconRNASeq_CBSX.Melanoma.scRNAseq_Endothelial"      
#> [49] "Endothelial_Subgroup.2"                              
#> [50] "DeconRNASeq_LM22_Eosinophils"                        
#> [51] "DeconRNASeq_LM22_Plasma.cells"                       
#> [52] "Plasma_Subgroup.1"                                   
#> [53] "DeconRNASeq_CBSX.HNSCC.scRNAseq_Myocytes"            
#> [54] "DeconRNASeq_CBSX.HNSCC.scRNAseq_CAF"                 
#> [55] "DeconRNASeq_CBSX.Melanoma.scRNAseq_CAF"              
#> [56] "CAF_Subgroup.1"                                      

if (FALSE) { # \dontrun{
# TF modules vs PROGENy pathways (downloads the PROGENy model from OmniPath)
data("counts.norm.tuto")
pathways <- compute.pathway.activity(counts.norm.tuto)
compute.modules.relationship(network.tuto[[1]],
                             pathways,
                             "Pathways_Progeny-TFs_Modules")
} # }
```
