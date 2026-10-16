# CellTFusion

`CellTFusion` integrates immune cell type deconvolution with
transcription factor (TF)–gene regulatory networks to characterize
immune cell states in the tumor microenvironment (TME) from bulk RNA-seq
data.

Starting from a count matrix, `CellTFusion` builds **cell groups** —
sets of deconvolution features related to the same TF co-activity module
to the same degree (a similar correlation with the module activity) —
and summarizes them into **latent factors**, a compact representation of
the TME that can be annotated as TME states, tested against clinical
variables and used as features for machine learning.

## Installation

To avoid GitHub API rate limit issues, set up a Personal Access Token
(PAT) before installing:

``` r

# install.packages(c("usethis", "gitcreds"))
usethis::create_github_token()
gitcreds::gitcreds_set()
```

Install `CellTFusion` from GitHub:

``` r

# install.packages("pak")
pak::pkg_install("VeraPancaldiLab/CellTFusion")
```

## Quick start

The
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
wrapper runs the whole pipeline in one call, using the example data
shipped with the package. Intermediate results and plots are saved in a
`Results/` folder in the working directory.

``` r

library(CellTFusion)

res <- CellTFusion(
  raw_counts     = CellTFusion::raw.counts.tuto,
  normalized     = TRUE,
  deconv_methods = c("Quantiseq", "Epidish"),
  file_name      = "Tutorial"
)

head(res$Latent_spaces$Z)  # latent factor scores (samples x factors)
```

## Pipeline

| Step | Function | Tutorial |
|----|----|----|
| Cell type deconvolution | [`multideconv::compute.deconvolution()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.html) | [Feature computation](https://verapancaldilab.github.io/CellTFusion/articles/a1_feature_computation.md) |
| TF activity inference | [`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md) | [Feature computation](https://verapancaldilab.github.io/CellTFusion/articles/a1_feature_computation.md) |
| TF co-activity modules | [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md) | [Feature computation](https://verapancaldilab.github.io/CellTFusion/articles/a1_feature_computation.md) |
| Pathway activity | [`compute.pathway.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md) | [Feature computation](https://verapancaldilab.github.io/CellTFusion/articles/a1_feature_computation.md) |
| Cell groups | [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md) | [Cell groups and latent factors](https://verapancaldilab.github.io/CellTFusion/articles/a2_cell_groups.md) |
| Latent factors | [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md) | [Cell groups and latent factors](https://verapancaldilab.github.io/CellTFusion/articles/a2_cell_groups.md) |
| Cell niches | [`compute_cells_niches()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_cells_niches.md) | [Cell groups and latent factors](https://verapancaldilab.github.io/CellTFusion/articles/a2_cell_groups.md) |
| Clinical associations | [`scores.stat.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.stat.analysis.md) | [Statistical analysis](https://verapancaldilab.github.io/CellTFusion/articles/a4_statistical_analysis.md) |
| Survival analysis | [`compute.survival.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.survival.analysis.md) | [Statistical analysis](https://verapancaldilab.github.io/CellTFusion/articles/a4_statistical_analysis.md) |
| Cross-validation with `pipeML` | [`prepare_celltfusion_folds()`](https://verapancaldilab.github.io/CellTFusion/reference/prepare_celltfusion_folds.md) | [Machine learning workflows](https://verapancaldilab.github.io/CellTFusion/articles/a5_machine_learning.md) |
| Projection of new cohorts | [`compute.test.set()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.test.set.md) (cell groups), [`project_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/project_factors.md) (latent factors) | [Machine learning workflows](https://verapancaldilab.github.io/CellTFusion/articles/a5_machine_learning.md) |
| Multi-cohort analysis | `CellTFusion(batch = ...)` | [Multi-cohort analysis](https://verapancaldilab.github.io/CellTFusion/articles/a6_batch_analysis.md) |

## Tutorials

Step-by-step tutorials are available in the **Articles** section of the
navigation bar:

- **[Feature
  computation](https://verapancaldilab.github.io/CellTFusion/articles/a1_feature_computation.md)**
  — cell type deconvolution, TF activity, TF co-activity modules and
  pathway activity
- **[Cell groups and latent
  factors](https://verapancaldilab.github.io/CellTFusion/articles/a2_cell_groups.md)**
  — build cell groups, extract latent factors and characterize cell
  niches
- **[Statistical
  analysis](https://verapancaldilab.github.io/CellTFusion/articles/a4_statistical_analysis.md)**
  — associations with clinical variables and survival
- **[Machine learning
  workflows](https://verapancaldilab.github.io/CellTFusion/articles/a5_machine_learning.md)**
  — use cell group scores as features and score independent cohorts
- **[Multi-cohort
  analysis](https://verapancaldilab.github.io/CellTFusion/articles/a6_batch_analysis.md)**
  — correct for cohort effects with `batch`

## Shiny app

`CellTFusion` includes an interactive app to run the pipeline on the
example data or on your own data:

``` r

shiny::runApp(system.file("shiny", package = "CellTFusion"))
```

## Citation

If you use `CellTFusion` in a scientific publication, please cite:

> Hurtado, M., & Pancaldi, V. (2026). *CellTFusion: A transcriptional
> regulatory network framework for the identification of functional
> multicellular states from bulk RNA-seq data.* bioRxiv.
> <https://doi.org/10.64898/2026.06.30.735682>
