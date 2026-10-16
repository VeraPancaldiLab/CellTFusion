# Prepare CellTFusion cross-validation folds for pipeML

Fold construction function for `pipeML::compute_features.training.ML()`
(argument `fold_construction_fun`) that computes `CellTFusion` cell
group scores within each cross-validation fold, without information
leakage: cell groups are learned with
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
on the training samples of each fold only, and the test samples of the
fold are scored on them with
[`compute.test.set()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.test.set.md).
The features of the models are the cell group scores; the latent factors
are a downstream characterization: they are not computed in the folds,
and are computed in final mode (unless `compute_factors = FALSE` is
passed) for the characterization of the final model; new samples can be
projected onto them with
[`project_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/project_factors.md).
It works for classification and survival tasks.

## Usage

``` r
prepare_celltfusion_folds(
  data,
  folds = NULL,
  bestune = NULL,
  deconv,
  raw_counts = NULL,
  batch = NULL,
  do_parallel = FALSE,
  workers = NULL,
  return = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  Data frame provided by `pipeML`: samples as rows and genes as columns
  (the `features_train` given to `compute_features.training.ML()`, i.e.
  the transposed count matrix), plus the outcome columns (`target` for
  classification, `time` and `event` for survival).

- folds:

  Named list with the training rows of each fold, provided by `pipeML`.

- bestune:

  Provided by `pipeML` when building the features of the final model;
  `NULL` while running the folds.

- deconv:

  Deconvolution matrix of all the samples (samples x features, e.g. from
  [`multideconv::compute.deconvolution()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.html)),
  with rows in the same order as `data`. Deconvolution is computed
  independently for each sample, so computing it once for all samples
  does not leak information between folds.

- raw_counts:

  Optional count matrix (genes x samples; TPM values if
  `normalized = FALSE` is passed to
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md))
  with the samples of `data`. Its columns are matched to the samples of
  `deconv` by name (an error is raised if they are not the same
  samples); without column names, they must be in the same order as
  `data`. `pipeML` converts the column names of `features_train` with
  [`make.names()`](https://rdrr.io/r/base/make.names.html) (e.g.
  `"HLA-A"` becomes `"HLA.A"`), so gene symbols in `data` may no longer
  match TF regulons and gene sets; passing `raw_counts` keeps the
  original gene symbols. If `NULL` (default), `t(data)` is used.

- batch:

  Optional vector with the cohort of each sample, in the same order as
  `data` (e.g. `coldata$Cohort`). If given,
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
  runs in multi-cohort mode (see
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md))
  and the test samples of each cohort are projected separately, z-scored
  on their own means and SDs. Default `NULL` (single cohort).

- do_parallel:

  Logical; whether to compute the folds in parallel. Default is FALSE
  (one after the other).

- workers:

  (Optional) Number of parallel workers when `do_parallel = TRUE`. If
  `NULL` (default), the number of available cores minus one is used.

- return:

  Logical; in final mode, whether
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
  saves its outputs and plots in `Results/`. Outputs are never saved for
  the folds. Default TRUE.

- verbose:

  Logical; whether
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
  prints its progress messages. Default TRUE.

- ...:

  Other arguments passed to
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
  (e.g. `normalized`, `TF.collection`, `min_targets_size`,
  `min_module_size`, `corr_mod`, `corr`, `pval`, `cancer_type`). `dt`,
  `tfs` and `pathways` are not allowed, as they must be computed within
  each fold.

## Value

- In fold mode (`bestune = NULL`): each fold is saved to
  `Results/fold_<fold name>.rds`, which `pipeML` reads back, as a list
  with `train_data` (cell group scores of the training samples plus the
  outcome columns), `test_data` (cell group scores of the test samples,
  plus `time` and `event` for survival), `obs_test` (observed outcome of
  the test samples), `rowIndex` (rows of the test samples) and
  `fold_name`. The list of folds is returned invisibly. If
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
  fails in a fold (e.g. no cell groups are found on its training
  samples, or a cell group cannot be scored in the test samples), that
  fold is skipped with a warning and is not saved; an error is raised
  only if all the folds fail.

- In final mode (`bestune` provided): a list with (1) the cell group
  scores of all the samples plus the outcome columns, (2) the
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
  result on all the samples, available as `Custom_output` in the
  `pipeML` result and used to compute the features of new samples with
  [`compute.test.set()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.test.set.md),
  and (3) `bestune`.

## Details

Pass the fixed arguments (`deconv`, and optionally `raw_counts`, `batch`
and any
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
argument) through `fold_construction_args_fixed`. Tunable arguments
(`fold_construction_args_tunable`) are not supported. With `batch`, a
warning names the folds and cohorts with fewer than 10 test samples
(their test samples are z-scored on their own means and SDs, which is
unreliable with few samples; leave-one-cohort-out folds avoid it). With
`do_parallel = TRUE`, the first fold is computed before the others start
in parallel, so that the TF and pathway collections cached in `Results/`
are downloaded only once; the package must be installed for the parallel
workers.

## Examples

``` r
if (FALSE) { # \dontrun{
raw_counts <- CellTFusion::raw.counts.tuto
traitdata  <- CellTFusion::traitdata.tuto
deconv     <- CellTFusion::deconv.tuto  # deconvolution of all the samples

# Classification
ml_res <- pipeML::compute_features.training.ML(
  features_train = t(raw_counts),
  task_type      = "classification",
  target_var     = traitdata$Best.Confirmed.Overall.Response,
  trait.positive = "PD",
  metric         = "AUROC",
  k_folds        = 5,
  n_rep          = 1,
  fold_construction_fun        = prepare_celltfusion_folds,
  fold_construction_args_fixed = list(deconv = deconv, raw_counts = raw_counts)
)

# Features of new samples: scored on the cell groups learned on all the training samples
final <- ml_res$Custom_output
features_test <- compute.test.set(final$Processed_deconvolution, final$Cell_groups,
                                  deconv_test = deconv_test)

# Survival (time and event of each sample)
surv_res <- pipeML::compute_features.training.ML(
  features_train = t(raw_counts),
  task_type      = "survival",
  time_var       = time,
  event_var      = event,
  k_folds        = 5,
  n_rep          = 1,
  fold_construction_fun        = prepare_celltfusion_folds,
  fold_construction_args_fixed = list(deconv = deconv, raw_counts = raw_counts)
)
} # }
```
