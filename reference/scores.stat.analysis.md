# Perform statistical analysis on scores using a specified test

Unified interface for testing associations between any score matrix and
a clinical/experimental trait. Accepts cell group scores, NMF latent
factors (output of
[`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md)),
or any samples x features matrix.

## Usage

``` r
scores.stat.analysis(
  scores,
  coldata,
  trait,
  test = c("fisher", "wilcox", "anova", "kruskal", "ttest"),
  pval = 0.05,
  padj = FALSE
)
```

## Arguments

- scores:

  A list, NMF output from
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md),
  or a score matrix (samples x features). When a list, the first element
  must be the score matrix; optional second and third elements (e.g. the
  output of
  [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md))
  are subsetted to the significant features and returned with their
  names.

- coldata:

  A data frame (or tibble) containing clinical or experimental metadata
  for samples. Must include the column specified in `trait`. Its rows
  are matched to the samples of the scores by row name; without sample
  names, they must be in the order of the scores. Samples with a missing
  trait are left out.

- trait:

  Character. The name of the column in `coldata` representing the
  clinical or experimental trait to test against (e.g., response,
  subtype, etc.). It must be categorical.

- test:

  Character. Statistical test to perform. One of:

  - `"fisher"` - Fisher's exact test (scores binarised at median)

  - `"wilcox"` - Wilcoxon rank-sum test (non-parametric, binary traits)

  - `"anova"` - One-way ANOVA (parametric, two or more groups)

  - `"kruskal"` - Kruskal-Wallis test (non-parametric, two or more
    groups)

  - `"ttest"` - Welch's t-test (parametric, binary traits)

  Only one can be used per call; defaults to the first one (`"fisher"`).

- pval:

  Numeric. Threshold for significance: a feature is significant if its
  (adjusted) p-value is `<= pval` (default: 0.05).

- padj:

  Logical; if TRUE, the p-values are corrected for multiple testing
  across the features with the Benjamini-Hochberg method, and the
  adjusted p-values are compared with the threshold. Default FALSE.

## Value

A list whose first element is the score matrix restricted to the
significant features (no column if none is significant). Optional second
and third elements (from the input list, e.g. the composition and
weights of the cell groups) are subsetted accordingly and keep their
names. The last element, `Statistics`, is a data frame with the test of
every feature (`feature`, `test`, `statistic_name`, `statistic`,
`p_value`, `p_adjusted`, `significant`).

## Details

The function automatically calls the corresponding statistical test
function based on the `test` argument:

- [`scores.fisher.test()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.fisher.test.md)

- [`scores.wilcox.test()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.wilcox.test.md)

- [`scores.anova.test()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.anova.test.md)

- [`scores.kruskal.test()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.kruskal.test.md)

- [`scores.ttest()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.ttest.md)

Each test produces both a statistical result and visual outputs (PDF
plots) stored in the `"Results/"` folder for the significant features.
These visualizations include the relevant test results (p-values)
annotated on the plots.

## Examples

``` r
data("network.tuto")
data("deconv_subgroups.tuto")
data("traitdata.tuto")
cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Dendrogram plots saved to: Results/Dendrogram_color_clusters.pdf
#> 44 of 897 features were not assigned to any cluster by the tree cut and are not part of any cell group (pink: 2, yellow: 8, red: 2, blue: 6, brown: 4, turquoise: 12, black: 4, green: 6)
sig <- scores.stat.analysis(cell_groups, traitdata.tuto, trait = "Enrollment.IC",
                            test = "kruskal", pval = 0.05, padj = TRUE)
#> 83 significant feature(s) after Kruskal-Wallis test (adjusted p-value (BH) <= 0.05): Dendrogram_pink_group_1, Dendrogram_pink_group_2, Dendrogram_pink_group_3, Dendrogram_pink_group_4, Dendrogram_pink_group_5, Dendrogram_pink_group_6, Dendrogram_pink_group_7, Dendrogram_pink_group_8, Dendrogram_pink_group_10, Dendrogram_pink_group_11, ...
head(sig$Statistics)
#>                   feature                test statistic_name statistic
#> 1 Dendrogram_pink_group_1 Kruskal-Wallis test              H 10.091884
#> 2 Dendrogram_pink_group_2 Kruskal-Wallis test              H  8.075064
#> 3 Dendrogram_pink_group_3 Kruskal-Wallis test              H 29.615321
#> 4 Dendrogram_pink_group_4 Kruskal-Wallis test              H 15.049547
#> 5 Dendrogram_pink_group_5 Kruskal-Wallis test              H  8.563492
#> 6 Dendrogram_pink_group_6 Kruskal-Wallis test              H 34.623612
#>        p_value   p_adjusted significant
#> 1 6.435395e-03 1.625577e-02        TRUE
#> 2 1.764095e-02 3.867440e-02        TRUE
#> 3 3.707787e-07 1.585079e-05        TRUE
#> 4 5.395509e-04 2.469786e-03        TRUE
#> 5 1.381852e-02 3.109166e-02        TRUE
#> 6 3.030946e-08 1.727639e-06        TRUE
names(sig$Composition)
#>  [1] "Dendrogram_pink_group_1"       "Dendrogram_pink_group_2"      
#>  [3] "Dendrogram_pink_group_3"       "Dendrogram_pink_group_4"      
#>  [5] "Dendrogram_pink_group_5"       "Dendrogram_pink_group_6"      
#>  [7] "Dendrogram_pink_group_7"       "Dendrogram_pink_group_8"      
#>  [9] "Dendrogram_pink_group_10"      "Dendrogram_pink_group_11"     
#> [11] "Dendrogram_yellow_group_5"     "Dendrogram_yellow_group_7"    
#> [13] "Dendrogram_yellow_group_10"    "Dendrogram_yellow_group_11"   
#> [15] "Dendrogram_yellow_group_13"    "Dendrogram_yellow_group_14"   
#> [17] "Dendrogram_yellow_group_15"    "Dendrogram_yellow_group_16"   
#> [19] "Dendrogram_yellow_group_17"    "Dendrogram_red_group_4"       
#> [21] "Dendrogram_red_group_6"        "Dendrogram_red_group_7"       
#> [23] "Dendrogram_red_group_10"       "Dendrogram_red_group_12"      
#> [25] "Dendrogram_red_group_13"       "Dendrogram_red_group_14"      
#> [27] "Dendrogram_red_group_16"       "Dendrogram_red_group_18"      
#> [29] "Dendrogram_red_group_19"       "Dendrogram_blue_group_1"      
#> [31] "Dendrogram_blue_group_2"       "Dendrogram_blue_group_4"      
#> [33] "Dendrogram_blue_group_7"       "Dendrogram_blue_group_11"     
#> [35] "Dendrogram_blue_group_15"      "Dendrogram_blue_group_17"     
#> [37] "Dendrogram_brown_group_2"      "Dendrogram_brown_group_6"     
#> [39] "Dendrogram_brown_group_7"      "Dendrogram_brown_group_11"    
#> [41] "Dendrogram_brown_group_13"     "Dendrogram_brown_group_14"    
#> [43] "Dendrogram_brown_group_16"     "Dendrogram_turquoise_group_1" 
#> [45] "Dendrogram_turquoise_group_2"  "Dendrogram_turquoise_group_3" 
#> [47] "Dendrogram_turquoise_group_4"  "Dendrogram_turquoise_group_7" 
#> [49] "Dendrogram_turquoise_group_8"  "Dendrogram_turquoise_group_9" 
#> [51] "Dendrogram_turquoise_group_11" "Dendrogram_turquoise_group_14"
#> [53] "Dendrogram_turquoise_group_16" "Dendrogram_turquoise_group_19"
#> [55] "Dendrogram_turquoise_group_21" "Dendrogram_turquoise_group_22"
#> [57] "Dendrogram_black_group_1"      "Dendrogram_black_group_3"     
#> [59] "Dendrogram_black_group_4"      "Dendrogram_black_group_8"     
#> [61] "Dendrogram_black_group_10"     "Dendrogram_black_group_12"    
#> [63] "Dendrogram_black_group_13"     "Dendrogram_black_group_14"    
#> [65] "Dendrogram_black_group_18"     "Dendrogram_black_group_21"    
#> [67] "Dendrogram_black_group_22"     "Dendrogram_black_group_23"    
#> [69] "Dendrogram_black_group_24"     "Dendrogram_green_group_1"     
#> [71] "Dendrogram_green_group_2"      "Dendrogram_green_group_3"     
#> [73] "Dendrogram_green_group_5"      "Dendrogram_green_group_9"     
#> [75] "Dendrogram_green_group_10"     "Dendrogram_green_group_12"    
#> [77] "Dendrogram_green_group_14"     "Dendrogram_green_group_15"    
#> [79] "Dendrogram_green_group_22"     "Dendrogram_green_group_23"    
#> [81] "Dendrogram_green_group_25"     "Dendrogram_green_group_26"    
#> [83] "Dendrogram_green_group_27"    
```
