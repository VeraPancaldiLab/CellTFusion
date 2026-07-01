# Machine Learning

``` r

library(CellTFusion)
#> 
#> 
```

Cell group scores computed by CellTFusion can be used directly as input
features for machine learning models to predict clinical traits. This
article demonstrates how to train and evaluate models using the
[`pipeML`](https://github.com/VeraPancaldiLab/pipeML) R package.

> **Installation:** `pipeML` is not on CRAN. Install it from GitHub
> before running the examples in this article:
>
> ``` r
>
> remotes::install_github("VeraPancaldiLab/pipeML")
> ```

## Train/test split

Start by splitting the dataset into training and testing cohorts:

``` r

raw.counts <- CellTFusion::raw.counts.tuto
traitdata  <- CellTFusion::traitdata.tuto

index <- caret::createDataPartition(
  traitdata[, "Best.Confirmed.Overall.Response"],
  times = 1, p = 0.8, list = FALSE
)

traitData_train  <- traitdata[index, ]
raw.counts_train <- raw.counts[, index]

traitData_test   <- traitdata[-index, ]
raw.counts_test  <- raw.counts[, -index]
```

## Run CellTFusion on the training set

``` r

res_training <- CellTFusion(
  raw.counts     = raw.counts_train,
  normalized     = TRUE,
  coldata        = traitData_train,
  trait          = "Best.Confirmed.Overall.Response",
  trait.positive = "PD",
  deconv_methods = c("Quantiseq", "Epidish"),
  file_name      = "TestRun",
  corr           = 0.7,
  pval           = 0.05,
  high_corr_groups = 0.85,
  return         = FALSE
)
```

## Train a machine learning model

Use `pipeML::compute_features.training.ML()` to train classifiers on the
cell group scores:

``` r

library(pipeML)

res <- pipeML::compute_features.training.ML(
  features_train = res_training$Cell_groups[[1]],
  target_var     = traitData_train$Best.Confirmed.Overall.Response,
  trait.positive = "PD",
  metric         = "AUROC",
  task_type      = "classification",
  stack          = FALSE,
  k_folds        = 2,
  n_rep          = 2,
  ncores         = 2,
  return         = FALSE
)
```

## Project cell groups to the test set

To apply the trained cell groups to an independent dataset, replicate
the deconvolution and cell group structure:

``` r

deconv_test <- multideconv::compute.deconvolution(
  raw.counts_test,
  methods    = c("Quantiseq", "Epidish"),
  normalized = TRUE,
  return     = FALSE
)

testing_set <- compute.test.set(
  res_training$Processed_deconvolution,
  res_training$Cell_groups,
  names(res_training$Cell_groups[[2]]),
  deconv_test
)
```

## Predict on the test set

``` r

pred <- pipeML::compute_prediction(
  res, testing_set,
  traitData_test$Best.Confirmed.Overall.Response,
  "PD", stack = FALSE
)

head(pred$Metrics[, 1:5])
pred$AUC$AUROC
```

## Custom k-fold cross-validation

In certain use cases it is essential to carefully design train/test
splits to prevent data leakage. This is especially important for feature
construction methods that depend on sample-level correlations or
deconvolution results — if such features are computed on the full
dataset before splitting, information from the test set can
inadvertently influence training.

`CellTFusion` provides `prepare_CellTFusion_folds()`, which is fully
compatible with `pipeML` and recomputes CellTFusion features inside each
fold. This ensures that:

- Test samples are never included in the feature learning process.
- Hyperparameter tuning can be performed safely without bias.
- Parallelization (via `foreach`/`doParallel`) reduces runtime.

``` r

universe <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)
paths    <- decoupleR::get_progeny(organism = "human", top = 500)

res_groups <- pipeML::compute_features.training.ML(
  t(raw.counts_train),
  traitData_train$Best.Confirmed.Overall.Response,
  trait.positive = "PD",
  metric         = "AUROC",
  stack          = FALSE,
  k_folds        = 3,
  n_rep          = 5,
  LODO           = FALSE,
  file_name      = "Test",
  ncores         = 2,
  return         = TRUE,
  fold_construction_fun = prepare_CellTFusion_folds,
  fold_construction_args_fixed = list(
    deconv     = deconv,
    universe   = universe,
    paths      = paths,
    ncores     = 2,
    normalized = TRUE,
    coldata    = traitData_train,
    trait      = "Best.Confirmed.Overall.Response",
    trait.positive = "PD"
  ),
  fold_construction_args_tunable = list(
    min_targets_size = c(5, 10, 15, 20),
    minMod           = c(5, 10, 15, 20),
    corr_mod         = c(0.7, 0.8, 0.9),
    corr             = c(0.7, 0.8, 0.9),
    high_corr_groups = 0.9
  )
)
```

`fold_construction_args_fixed` holds parameters that stay constant
across folds; `fold_construction_args_tunable` defines hyperparameters
to tune during cross-validation.

👉 For full details on `pipeML` arguments, see the [pipeML
tutorial](https://verapancaldilab.github.io/pipeML/articles/pipeML.html).
