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
  corr_type = "spearman",
  vertical = FALSE,
  return = TRUE,
  plot_grid = FALSE,
  ncol_grid = NULL
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
  Pearson or Spearman according to `corr_type`) are computed controlling
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

- corr_type:

  Type of correlation to compute: "spearman" (default) or "pearson". It
  also sets the type of the partial correlations computed with `batch`.

- vertical:

  Logical; if TRUE, modules are on the x-axis and the `matB` features on
  the y-axis. Otherwise (default), the `matB` features are on the x-axis
  and modules on the y-axis.

- return:

  Logical; if TRUE, saves the heatmap plot as a PDF (`file_name` is
  needed). Default is TRUE.

- plot_grid:

  Logical; if TRUE, generates per-pair scatter grid plots for
  significant associations, saved to
  "Results/\<file_name\>\_scatter_grid.svg" (the 100 strongest pairs at
  most). With Spearman, the ranks of the values are shown, as the
  correlation is computed on them. With `batch`, the values shown are
  corrected for batch, as the correlations are. Default is FALSE.

- ncol_grid:

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

If `return = TRUE` and `file_name` is given, a heatmap of the
correlations is saved to "Results/\<file_name\>.pdf" (the ".pdf"
extension is added if missing). Its size is computed from the number of
features. With up to 40 features per axis, the significant cells show
the correlation and its significance (\* \< 0.05, \*\* \< 0.01, \*\*\*
\< 0.001, \*\*\*\* \< 0.0001); with more, they are marked with an
asterisk.

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
                                    return = FALSE,
                                    pval = 0.01)
corr$correlation[1:3, 1:3]
#>        Quantiseq_B.cells DeconRNASeq_CBSX.Melanoma.scRNAseq_B.cells
#> pink          -0.3092277                                 0.05039945
#> yellow         0.3771651                                -0.12406919
#> red            0.1006260                                -0.23509345
#>        DeconRNASeq_CBSX.NSCLC.PBMCs.scRNAseq_B.cells
#> pink                                     -0.05170327
#> yellow                                    0.27144427
#> red                                       0.16008769
corr$significant[[1]]
#>  [1] "Quantiseq_B.cells"                              
#>  [2] "DWLS_CCLE.TIL10_B.cells"                        
#>  [3] "CBSX_CBSX.NSCLC.PBMCs.scRNAseq_B.cells"         
#>  [4] "CBSX_CCLE.TIL10_B.cells"                        
#>  [5] "B.cells_Subgroup.6"                             
#>  [6] "DeconRNASeq_LM22_B.naive.cells"                 
#>  [7] "Epidish_LM22_B.memory.cells"                    
#>  [8] "CBSX_LM22_B.memory.cells"                       
#>  [9] "Macrophages.cells_Subgroup.2"                   
#> [10] "Macrophages.M0_Subgroup.2"                      
#> [11] "DeconRNASeq_LM22_Macrophages.M1"                
#> [12] "Macrophages.M1_Subgroup.2"                      
#> [13] "Macrophages.M1_Subgroup.3"                      
#> [14] "Quantiseq_Macrophages.M2"                       
#> [15] "CBSX_CCLE.TIL10_Macrophages.M2"                 
#> [16] "Macrophages.M2_Subgroup.5"                      
#> [17] "DeconRNASeq_CBSX.NSCLC.PBMCs.scRNAseq_Monocytes"
#> [18] "CBSX_CCLE.TIL10_Monocytes"                      
#> [19] "Monocytes_Subgroup.2"                           
#> [20] "Monocytes_Subgroup.5"                           
#> [21] "Monocytes_Subgroup.6"                           
#> [22] "Neutrophils_Subgroup.1"                         
#> [23] "Quantiseq_NK.cells"                             
#> [24] "Epidish_BPRNACanProMet_NK.cells"                
#> [25] "DeconRNASeq_CBSX.Melanoma.scRNAseq_NK.cells"    
#> [26] "DeconRNASeq_CBSX.NSCLC.PBMCs.scRNAseq_NK.cells" 
#> [27] "Epidish_CCLE.TIL10_NK.cells"                    
#> [28] "DWLS_CCLE.TIL10_NK.cells"                       
#> [29] "CBSX_BPRNACanProMet_NK.cells"                   
#> [30] "NK.activated_Subgroup.1"                        
#> [31] "DeconRNASeq_CBSX.NSCLC.PBMCs.scRNAseq_NKT.cells"
#> [32] "Epidish_CCLE.TIL10_CD4.cells"                   
#> [33] "DWLS_CCLE.TIL10_CD4.cells"                      
#> [34] "CBSX_CCLE.TIL10_CD4.cells"                      
#> [35] "CD4.memory.resting_Subgroup.1"                  
#> [36] "DeconRNASeq_LM22_CD4.naive"                     
#> [37] "DeconRNASeq_BPRNACan_CD8.cells"                 
#> [38] "DeconRNASeq_CBSX.HNSCC.scRNAseq_CD8.cells"      
#> [39] "DeconRNASeq_CCLE.TIL10_CD8.cells"               
#> [40] "Epidish_LM22_CD8.cells"                         
#> [41] "DWLS_TIL10_CD8.cells"                           
#> [42] "DeconRNASeq_LM22_T.cells.helper"                
#> [43] "CBSX_LM22_T.cells.helper"                       
#> [44] "T.cells.helper_Subgroup.1"                      
#> [45] "Dendritic.activated.cells_Subgroup.1"           
#> [46] "Dendritic.resting.cells_Subgroup.1"             
#> [47] "DeconRNASeq_CCLE.TIL10_Cancer"                  
#> [48] "Cancer_Subgroup.1"                              
#> [49] "Cancer_Subgroup.3"                              
#> [50] "Endothelial_Subgroup.2"                         
#> [51] "DeconRNASeq_CBSX.HNSCC.scRNAseq_Myocytes"       
#> [52] "DeconRNASeq_LM22_Mast.activated.cells"          
#> [53] "Mast.activated.cells_Subgroup.1"                
#> [54] "DeconRNASeq_LM22_Mast.resting.cells"            
#> [55] "Mast.resting.cells_Subgroup.1"                  
#> [56] "DeconRNASeq_CBSX.HNSCC.scRNAseq_CAF"            
#> [57] "CAF_Subgroup.1"                                 

if (FALSE) { # \dontrun{
# TF modules vs PROGENy pathways (downloads the PROGENy model from OmniPath)
data("counts.norm.tuto")
pathways <- compute.pathway.activity(counts.norm.tuto)
compute.modules.relationship(network.tuto[[1]],
                             pathways,
                             "Pathways_Progeny-TFs_Modules")
} # }
```
