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
  method = c("fisher", "wilcox", "anova", "kruskal", "ttest"),
  pval = 0.05,
  p_adjust = "none"
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

- method:

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

- p_adjust:

  Character. Method to correct the p-values for multiple testing across
  the features (see
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html), e.g.
  `"BH"`). Default `"none"`.

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
function based on the `method` argument:

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
#> 20 of 764 features were not assigned to any cluster by the tree cut and are not part of any cell group (blue: 6, pink: 2, red: 2, yellow: 8, black: 2)
sig <- scores.stat.analysis(cell_groups, traitdata.tuto, trait = "Enrollment.IC",
                            method = "kruskal", pval = 0.05, p_adjust = "BH")
#> 73 significant feature(s) after Kruskal-Wallis test (adjusted p-value (BH) <= 0.05): Dendrogram_brown_group_3, Dendrogram_brown_group_4, Dendrogram_brown_group_5, Dendrogram_brown_group_6, Dendrogram_brown_group_7, Dendrogram_brown_group_11, Dendrogram_purple_group_2, Dendrogram_blue_group_3, Dendrogram_blue_group_9, Dendrogram_blue_group_11, ...
head(sig$Statistics)
#>                    feature                test statistic_name statistic
#> 1 Dendrogram_brown_group_1 Kruskal-Wallis test              H  2.695207
#> 2 Dendrogram_brown_group_2 Kruskal-Wallis test              H  2.612820
#> 3 Dendrogram_brown_group_3 Kruskal-Wallis test              H 10.356225
#> 4 Dendrogram_brown_group_4 Kruskal-Wallis test              H 12.480027
#> 5 Dendrogram_brown_group_5 Kruskal-Wallis test              H 26.813980
#> 6 Dendrogram_brown_group_6 Kruskal-Wallis test              H 31.124228
#>        p_value   p_adjusted significant
#> 1 2.598623e-01 3.433043e-01       FALSE
#> 2 2.707904e-01 3.544893e-01       FALSE
#> 3 5.638640e-03 1.249176e-02        TRUE
#> 4 1.949829e-03 5.105006e-03        TRUE
#> 5 1.504590e-06 1.083305e-05        TRUE
#> 6 1.743652e-07 3.138574e-06        TRUE
names(sig$Composition)
#>  [1] "Dendrogram_brown_group_3"     "Dendrogram_brown_group_4"    
#>  [3] "Dendrogram_brown_group_5"     "Dendrogram_brown_group_6"    
#>  [5] "Dendrogram_brown_group_7"     "Dendrogram_brown_group_11"   
#>  [7] "Dendrogram_purple_group_2"    "Dendrogram_blue_group_3"     
#>  [9] "Dendrogram_blue_group_9"      "Dendrogram_blue_group_11"    
#> [11] "Dendrogram_blue_group_13"     "Dendrogram_blue_group_15"    
#> [13] "Dendrogram_pink_group_1"      "Dendrogram_pink_group_3"     
#> [15] "Dendrogram_pink_group_9"      "Dendrogram_pink_group_10"    
#> [17] "Dendrogram_pink_group_11"     "Dendrogram_pink_group_12"    
#> [19] "Dendrogram_pink_group_13"     "Dendrogram_pink_group_15"    
#> [21] "Dendrogram_pink_group_16"     "Dendrogram_pink_group_17"    
#> [23] "Dendrogram_pink_group_19"     "Dendrogram_pink_group_20"    
#> [25] "Dendrogram_pink_group_21"     "Dendrogram_magenta_group_1"  
#> [27] "Dendrogram_magenta_group_2"   "Dendrogram_magenta_group_3"  
#> [29] "Dendrogram_magenta_group_4"   "Dendrogram_magenta_group_6"  
#> [31] "Dendrogram_magenta_group_7"   "Dendrogram_magenta_group_8"  
#> [33] "Dendrogram_magenta_group_9"   "Dendrogram_magenta_group_10" 
#> [35] "Dendrogram_magenta_group_11"  "Dendrogram_magenta_group_12" 
#> [37] "Dendrogram_magenta_group_14"  "Dendrogram_magenta_group_15" 
#> [39] "Dendrogram_red_group_1"       "Dendrogram_red_group_2"      
#> [41] "Dendrogram_red_group_6"       "Dendrogram_red_group_8"      
#> [43] "Dendrogram_red_group_13"      "Dendrogram_red_group_14"     
#> [45] "Dendrogram_red_group_15"      "Dendrogram_red_group_21"     
#> [47] "Dendrogram_red_group_22"      "Dendrogram_red_group_23"     
#> [49] "Dendrogram_green_group_2"     "Dendrogram_green_group_4"    
#> [51] "Dendrogram_green_group_7"     "Dendrogram_green_group_8"    
#> [53] "Dendrogram_green_group_10"    "Dendrogram_yellow_group_1"   
#> [55] "Dendrogram_yellow_group_2"    "Dendrogram_yellow_group_3"   
#> [57] "Dendrogram_yellow_group_4"    "Dendrogram_yellow_group_6"   
#> [59] "Dendrogram_yellow_group_10"   "Dendrogram_yellow_group_11"  
#> [61] "Dendrogram_yellow_group_13"   "Dendrogram_yellow_group_14"  
#> [63] "Dendrogram_yellow_group_16"   "Dendrogram_yellow_group_17"  
#> [65] "Dendrogram_black_group_2"     "Dendrogram_black_group_3"    
#> [67] "Dendrogram_black_group_6"     "Dendrogram_black_group_9"    
#> [69] "Dendrogram_black_group_15"    "Dendrogram_turquoise_group_1"
#> [71] "Dendrogram_turquoise_group_2" "Dendrogram_turquoise_group_3"
#> [73] "Dendrogram_turquoise_group_8"
```
