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
  survival.data,
  PFS,
  PFS_event,
  file_name = NULL,
  features = NULL,
  p.value = 0.05,
  thres = 0.5,
  group_column = NULL,
  p_adjust = "none"
)
```

## Arguments

- survival.data:

  A data frame (or tibble) of clinical/survival metadata, with samples
  as rows.

- PFS:

  Character. Column name in `survival.data` with the survival/follow-up
  time.

- PFS_event:

  Character. Column name in `survival.data` with the event indicator (1
  = event occurred, 0 = censored).

- file_name:

  Optional character. Suffix used when saving Kaplan-Meier plots to
  `Results/`.

- features:

  Optional. A samples x features numeric matrix or data frame (e.g.
  `latent_spaces$Z`). If provided (and `group_column` is `NULL`), each
  feature is tested individually. Rows must be the same samples, in the
  same order, as the rows of `survival.data` (checked by name). Mutually
  exclusive with `group_column`.

- p.value:

  Numeric. Log-rank test threshold used to keep a feature as significant
  when `features` is used: (adjusted) p-value `<= p.value`. Default is
  0.05.

- thres:

  Numeric between 0 and 1. Quantile cutoff used to split each feature
  into High/Low groups when `features` is used (High = at or above the
  cutoff). Default is 0.5 (median split).

- group_column:

  Optional character. Column name in `survival.data` defining a
  predefined categorical grouping (e.g. cell group membership, cluster,
  or trait class). Samples without a value are left out (with a
  message). Mutually exclusive with `features`.

- p_adjust:

  Character. Method to correct the log-rank p-values for multiple
  testing across the features (see
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html), e.g.
  `"BH"`), when `features` is used. Default `"none"`.

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

If `features` is provided, a named list with one `Surv() ~ feature`
formula per significant feature, and a last element `Statistics`: a data
frame with the test of every feature (`feature`, `cutoff`, `n_high`,
`n_low`, `p_value`, `p_adjusted`, `significant`). When no feature is
significant, the list only contains `Statistics` (with a message).

A Kaplan-Meier plot (with risk table) is saved as an SVG file to
`Results/SurvPlot_<group-or-feature>_<file_name>.svg`
(`Results/SurvPlot_<group-or-feature>.svg` without `file_name`;
characters not allowed in file names, e.g. "/", replaced by "\_"):
always for `group_column`, and one per significant feature for
`features`.

## Details

Requires the `survival`, `survminer`, and `gridExtra` packages (see
`Suggests`).

## Examples

``` r
if (FALSE) { # \dontrun{
# traitdata: clinical data with the survival columns (here "PFS" and "PFS_event");
# the tutorial data (traitdata.tuto) does not include them

# Predefined clinical/cell groups
compute.survival.analysis(
  survival.data = traitdata,
  PFS           = "PFS",
  PFS_event     = "PFS_event",
  group_column  = "Best.Confirmed.Overall.Response",
  file_name     = "Tutorial"
)

# Automatic screening of latent factor scores
compute.survival.analysis(
  survival.data = traitdata,
  PFS           = "PFS",
  PFS_event     = "PFS_event",
  features      = res$Latent_spaces$Z,
  p.value       = 0.05,
  thres         = 0.5,
  file_name     = "Tutorial"
)
} # }
```
