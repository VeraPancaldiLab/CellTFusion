# Prepare CellTFusion cross-validation folds for pipeML

Fold construction function for `pipeML::compute_features.training.ML()`
(argument `fold_construction_fun`) that computes `CellTFusion` latent
factors within each cross-validation fold, without information leakage:
cell groups and latent factors are learned with
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
on the training samples of each fold only, and the test samples of the
fold are projected onto them with
[`project_test_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/project_test_factors.md).
It works for classification and survival tasks.

## Usage

``` r
prepare_celltfusion_folds(
  data,
  folds = NULL,
  bestune = NULL,
  deconv,
  raw.counts = NULL,
  coldata = NULL,
  batch = FALSE,
  batch_id = NULL,
  ncores = 1,
  return = FALSE,
  verbose = FALSE,
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

- raw.counts:

  Optional count matrix (genes x samples) with the samples in the same
  order as `data`. `pipeML` converts the column names of
  `features_train` with
  [`make.names()`](https://rdrr.io/r/base/make.names.html) (e.g.
  `"HLA-A"` becomes `"HLA.A"`), so gene symbols in `data` may no longer
  match TF regulons and gene sets; passing `raw.counts` keeps the
  original gene symbols. If `NULL` (default), `t(data)` is used.

- coldata:

  Optional data frame of sample metadata with rows in the same order as
  `data`. Required when `batch = TRUE`.

- batch:

  Logical; whether to run
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
  in multi-cohort mode (see
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)).
  The test samples of each cohort are then projected separately,
  centered on their own means. Default FALSE.

- batch_id:

  Column of `coldata` with the cohort of each sample. Required when
  `batch = TRUE`.

- ncores:

  Integer. Number of folds computed in parallel. Default 1 (sequential).

- return:

  Logical; in final mode, whether
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
  saves its outputs and plots in `Results/`. Outputs are never saved for
  the folds. Default FALSE.

- verbose:

  Logical; whether
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
  prints its progress messages. Default FALSE.

- ...:

  Other arguments passed to
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
  (e.g. `normalized`, `TF.collection`, `min_targets_size`, `minMod`,
  `corr_mod`, `corr`, `pval`, `cancer_type`). `dt`, `tfs` and `pathways`
  are not allowed, as they must be computed within each fold.

## Value

- In fold mode (`bestune = NULL`): each fold is saved to
  `Results/fold_<fold name>.rds`, which `pipeML` reads back, as a list
  with `train_data` (latent factors of the training samples plus the
  outcome columns), `test_data` (latent factors of the test samples,
  plus `time` and `event` for survival), `obs_test` (observed outcome of
  the test samples), `rowIndex` (rows of the test samples) and
  `fold_name`. The list of folds is returned invisibly.

- In final mode (`bestune` provided): a list with (1) the latent factors
  of all the samples plus the outcome columns, (2) the
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
  result on all the samples, available as `Custom_output` in the
  `pipeML` result and used to compute the features of new samples with
  [`project_test_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/project_test_factors.md),
  and (3) `bestune`.

## Details

Pass the fixed arguments (`deconv`, and optionally `raw.counts`,
`coldata`, `batch`, `batch_id` and any
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
argument) through `fold_construction_args_fixed`. Tunable arguments
(`fold_construction_args_tunable`) are not supported. With `ncores > 1`,
the first fold is computed before the others start in parallel, so that
the TF and pathway collections cached in `Results/` are downloaded only
once; the package must be installed for the parallel workers.

## Examples

``` r
if (FALSE) { # \dontrun{
raw.counts <- CellTFusion::raw.counts.tuto
traitdata  <- CellTFusion::traitdata.tuto
deconv     <- CellTFusion::deconv.tuto  # deconvolution of all the samples

# Classification
ml_res <- pipeML::compute_features.training.ML(
  features_train = t(raw.counts),
  task_type      = "classification",
  target_var     = traitdata$Best.Confirmed.Overall.Response,
  trait.positive = "PD",
  metric         = "AUROC",
  k_folds        = 5,
  n_rep          = 1,
  fold_construction_fun        = prepare_celltfusion_folds,
  fold_construction_args_fixed = list(deconv = deconv, raw.counts = raw.counts)
)

# Features of new samples: projected onto the model trained on all the training samples
features_test <- project_test_factors(ml_res$Custom_output, deconv_test)

# Survival (time and event of each sample)
surv_res <- pipeML::compute_features.training.ML(
  features_train = t(raw.counts),
  task_type      = "survival",
  time_var       = time,
  event_var      = event,
  k_folds        = 5,
  n_rep          = 1,
  fold_construction_fun        = prepare_celltfusion_folds,
  fold_construction_args_fixed = list(deconv = deconv, raw.counts = raw.counts)
)
} # }
```
