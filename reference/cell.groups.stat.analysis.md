# Perform statistical analysis on cell group scores using a specified test

This function provides a unified interface to perform one of several
statistical tests (Fisher’s exact, Wilcoxon rank-sum, ANOVA,
Kruskal–Wallis, or Student’s t-test) on cell group scores in relation to
a given clinical or experimental trait. The choice of test is specified
by the `method` argument. Each test identifies significant associations
and saves a corresponding visualization for each significant feature in
the "Results/" directory.

## Usage

``` r
cell.groups.stat.analysis(
  cell.groups,
  coldata,
  trait,
  method = c("fisher", "wilcox", "anova", "kruskal", "ttest"),
  pval = 0.05
)
```

## Arguments

- cell.groups:

  A list where:

  - The first element is a data frame or matrix containing cell group
    scores (rows = samples, columns = cell groups).

  - The second and third elements contain corresponding metadata or
    annotations for these groups (e.g., names, features, etc.).

- coldata:

  A data frame containing clinical or experimental metadata for samples.
  Must include the column specified in `trait`.

- trait:

  Character. The name of the column in `coldata` representing the
  clinical or experimental trait to test against (e.g., response,
  subtype, etc.).

- method:

  Character. Statistical test to perform. One of:

  - `"fisher"` — Fisher’s exact test (for categorical data)

  - `"wilcox"` — Wilcoxon rank-sum test (non-parametric, binary traits)

  - `"anova"` — One-way ANOVA (parametric, \>2 groups)

  - `"kruskal"` — Kruskal–Wallis test (non-parametric, \>2 groups)

  - `"ttest"` — Student’s t-test (parametric, binary traits)

  Defaults to all available options, but only one can be used per call.

- pval:

  Numeric. P-value threshold for significance (default: 0.05).

## Value

A list of significant cell groups, where:

- The first element contains the subset of the original score matrix for
  significant features.

- The second and third elements contain associated metadata or feature
  annotations.

Returns `NULL` if no significant features are found.

## Details

The function automatically calls the corresponding statistical test
function based on the `method` argument:

- [`cell.groups.fisher.test()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.fisher.test.md)

- [`cell.groups.wilcox.test()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.wilcox.test.md)

- [`cell.groups.anova.test()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.anova.test.md)

- [`cell.groups.kruskal.test()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.kruskal.test.md)

- [`cell.groups.ttest()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.ttest.md)

Each test produces both a statistical result and visual outputs (PDF
plots) stored in the `"Results/"` folder. These visualizations include
the relevant test results (p-values) annotated on the plots.

## Examples

``` r
if (FALSE) { # \dontrun{
# Example usage:
sig.groups <- cell.groups.stat.analysis(cell.groups, coldata, trait = "response",
                                        method = "kruskal", pval = 0.05)
} # }
```
