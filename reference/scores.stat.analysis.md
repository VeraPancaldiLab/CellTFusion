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
  pval = 0.05
)
```

## Arguments

- scores:

  A list, NMF output from
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md),
  or a score matrix (samples x features). When a list, the first element
  must be the score matrix; optional second and third elements are
  passed through to the returned result unchanged.

- coldata:

  A data frame containing clinical or experimental metadata for samples.
  Must include the column specified in `trait`.

- trait:

  Character. The name of the column in `coldata` representing the
  clinical or experimental trait to test against (e.g., response,
  subtype, etc.).

- method:

  Character. Statistical test to perform. One of:

  - `"fisher"` - Fisher's exact test (scores binarised at median)

  - `"wilcox"` - Wilcoxon rank-sum test (non-parametric, binary traits)

  - `"anova"` - One-way ANOVA (parametric, \>2 groups)

  - `"kruskal"` - Kruskal-Wallis test (non-parametric, \>2 groups)

  - `"ttest"` - Student's t-test (parametric, binary traits)

  Defaults to all available options, but only one can be used per call.

- pval:

  Numeric. P-value threshold for significance (default: 0.05).

## Value

A list of significant features, where the first element contains the
subset of the original score matrix for significant features. Optional
second and third elements (from the input list) are subsetted
accordingly. Returns `NULL` if no significant features are found.

## Details

The function automatically calls the corresponding statistical test
function based on the `method` argument:

- [`scores.fisher.test()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.fisher.test.md)

- [`scores.wilcox.test()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.wilcox.test.md)

- [`scores.anova.test()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.anova.test.md)

- [`scores.kruskal.test()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.kruskal.test.md)

- [`scores.ttest()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.ttest.md)

Each test produces both a statistical result and visual outputs (PDF
plots) stored in the `"Results/"` folder. These visualizations include
the relevant test results (p-values) annotated on the plots.

## Examples

``` r
if (FALSE) { # \dontrun{
# Cell group scores
sig <- scores.stat.analysis(cell.groups, coldata, trait = "response",
                            method = "kruskal", pval = 0.05)

# NMF latent factors
nmf <- compute.latent.factors(cell.groups)
sig <- scores.stat.analysis(nmf, coldata, trait = "response", method = "anova")
} # }
```
