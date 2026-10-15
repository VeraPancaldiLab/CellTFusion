# Kaplan-Meier survival analysis on clinical groups or CellTFusion features

Performs Kaplan-Meier survival analysis and log-rank testing, either on
a predefined grouping variable (e.g. a clinical risk group or a
supervised cell group split), or automatically for each column of a
feature matrix (e.g. latent factor scores from
[`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md),
TF module scores from
[`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md),
or cell group scores from
[`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)).
In the latter case, samples are split into "High"/"Low" groups per
feature using a quantile cutoff, and only features with a significant
log-rank test are returned.

## Usage

``` r
compute.survival.analysis(
  coldata,
  time_column,
  event_column,
  file_name = NULL,
  scores = NULL,
  pval = 0.05,
  quantile_cutoff = 0.5,
  group_column = NULL,
  padj = FALSE
)
```

## Arguments

- coldata:

  A data frame (or tibble) of clinical/survival metadata, with samples
  as rows.

- time_column:

  Character. Column name in `coldata` with the survival/follow-up time.

- event_column:

  Character. Column name in `coldata` with the event indicator (1 =
  event occurred, 0 = censored).

- file_name:

  Optional character. Suffix used when saving Kaplan-Meier plots to
  `Results/`.

- scores:

  Optional. A samples x features numeric matrix or data frame (e.g.
  `latent_spaces$Z`). If provided (and `group_column` is `NULL`), each
  feature is tested individually. Rows must be the same samples, in the
  same order, as the rows of `coldata` (checked by name). Mutually
  exclusive with `group_column`.

- pval:

  Numeric. Log-rank test threshold used to keep a feature as significant
  when `scores` is used: (adjusted) p-value `<= pval`. Default is 0.05.

- quantile_cutoff:

  Numeric between 0 and 1. Quantile cutoff used to split each feature
  into High/Low groups when `scores` is used (High = at or above the
  cutoff). Default is 0.5 (median split).

- group_column:

  Optional character. Column name in `coldata` defining a predefined
  categorical grouping (e.g. cell group membership, cluster, or trait
  class). Samples without a value are left out (with a message).
  Mutually exclusive with `scores`.

- padj:

  Logical; if TRUE, the log-rank p-values are corrected for multiple
  testing across the features with the Benjamini-Hochberg method, when
  `scores` is used. Default FALSE.

## Value

If `group_column` is provided, a list with:

- km_plot:

  The
  [`survminer::ggsurvplot`](https://rdrr.io/pkg/survminer/man/ggsurvplot.html)
  object.

- median_PFS:

  Median survival time per group
  ([`survminer::surv_median()`](https://rdrr.io/pkg/survminer/man/surv_median.html)).

- p_value:

  Log-rank test p-value across groups.

- median_follow_up:

  Median follow-up time (reverse Kaplan-Meier).

If `scores` is provided, a named list with one `Surv() ~ feature`
formula per significant feature, and a last element `Statistics`: a data
frame with the test of every feature (`feature`, `cutoff`, `n_high`,
`n_low`, `p_value`, `p_adjusted`, `significant`). When no feature is
significant, the list only contains `Statistics` (with a message).

A Kaplan-Meier plot (with risk table) is saved as an SVG file to
`Results/SurvPlot_<group-or-feature>_<file_name>.svg`
(`Results/SurvPlot_<group-or-feature>.svg` without `file_name`;
characters not allowed in file names, e.g. "/", replaced by "\_"):
always for `group_column`, and one per significant feature for `scores`.

## Details

Requires the `survival`, `survminer`, and `gridExtra` packages (see
`Suggests`).

## Examples

``` r
if (FALSE) { # \dontrun{
# Overall survival of the tutorial data (os: months, censOS: 1 = death)
data("traitdata.tuto")

# Predefined clinical groups
compute.survival.analysis(
  coldata       = traitdata.tuto,
  time_column   = "os",
  event_column  = "censOS",
  group_column  = "Best.Confirmed.Overall.Response",
  file_name     = "Tutorial"
)

# Automatic screening of latent factor scores (res: output of CellTFusion() on the tutorial data)
compute.survival.analysis(
  coldata       = traitdata.tuto,
  time_column   = "os",
  event_column  = "censOS",
  scores        = res$Latent_spaces$Z,
  pval          = 0.05,
  quantile_cutoff = 0.5,
  file_name     = "Tutorial"
)
} # }
```
