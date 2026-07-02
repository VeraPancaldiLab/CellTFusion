# One-step Pipeline

``` r

library(CellTFusion)
#> 
#> 
```

[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
is an all-in-one wrapper that runs the full pipeline in a single call:

1.  Log-TPM normalisation
2.  Cell-type deconvolution (`multideconv`)
3.  TF activity inference (`compute.TFs.activity`)
4.  TF module construction (`compute.WTCNA`)
5.  Pathway activity scoring (`compute.pathway.activity`)
6.  Deconvolution dimensionality reduction
7.  Cell group identification (`construct_cell_groups`)
8.  Latent factor extraction (`compute.latent_factors`)
9.  Cell niche derivation (`compute_cells_niches`)
10. Hallmark GSEA per factor (`compute_factor_gsea`)
11. Meta-program mapping (`map_factors_to_metaprograms`)

For a step-by-step explanation of each stage see: [Feature
Computation](https://verapancaldilab.github.io/CellTFusion/articles/01-feature-computation.md)
· [Cell
Groups](https://verapancaldilab.github.io/CellTFusion/articles/02-cell-groups.md)
· [TME
States](https://verapancaldilab.github.io/CellTFusion/articles/03-tme-states.md)
· [Statistical
Analysis](https://verapancaldilab.github.io/CellTFusion/articles/04-analysis.md)

------------------------------------------------------------------------

## Unsupervised mode

Use `task = "unsupervised"` when no clinical contrast is available. This
is the most common entry point.

``` r

raw.counts <- CellTFusion::raw.counts.tuto
traitdata  <- CellTFusion::traitdata.tuto

res <- CellTFusion(
  raw.counts     = raw.counts,
  normalized     = TRUE,
  coldata        = traitdata,
  task           = "unsupervised",
  deconv_methods = c("Quantiseq", "Epidish"),
  TF.collection  = "CollecTRI",
  cancer_type    = "skcm",
  corr           = 0.7,
  corr_mod       = 0.9,
  pval           = 0.05,
  file_name      = "Tutorial",
  return         = TRUE
)
```

## Supervised mode

Use `task = "supervised"` when you have a clinical contrast of interest.
TF activity is computed on differentially expressed genes (DEGs) between
the two groups, focusing the analysis on trait-relevant TFs.

Requires:

- `contrast` — column name in `coldata` defining the comparison groups
- `ref_level` — the reference level (baseline) for the contrast

> **Note:** supervised mode requires raw (non-normalised) counts since
> it runs differential expression analysis internally.

``` r

res <- CellTFusion(
  raw.counts     = raw.counts,
  normalized     = FALSE,
  coldata        = traitdata,
  task           = "supervised",
  contrast       = "Best.Confirmed.Overall.Response",
  ref_level      = "PD",
  deconv_methods = c("Quantiseq", "Epidish"),
  TF.collection  = "CollecTRI",
  cancer_type    = "skcm",
  corr           = 0.7,
  corr_mod       = 0.9,
  pval           = 0.05,
  file_name      = "Tutorial_supervised",
  return         = TRUE
)
```

## Multi-cohort (batch) mode

When samples come from multiple cohorts, set `batch = TRUE` and specify
the cohort column in `coldata`. TF activity is computed per cohort and
consensus WTCNA is applied.

``` r

res <- CellTFusion(
  raw.counts     = raw.counts,
  normalized     = FALSE,
  deconv         = deconv,
  coldata        = traitdata,
  task           = "unsupervised",
  batch          = TRUE,
  batch_id       = "Cohort",
  TF.collection  = "ARACNE",
  cancer_type    = "skcm",
  corr           = 0.7,
  pval           = 0.05,
  file_name      = "Tutorial_batch",
  return         = TRUE
)
```

## Output structure

[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
returns a named list with 10 elements:

| Element | Description |
|----|----|
| `$Deconvolution` | Raw deconvolution proportions matrix |
| `$TFs_matrix` | TF activity matrix (samples x TFs) |
| `$TF_network` | WTCNA network object |
| `$Pathways_scores` | PROGENy pathway activity scores |
| `$Processed_deconvolution` | Reduced deconvolution subgroups |
| `$Cell_groups` | Cell group scores, compositions, and loadings |
| `$Latent_spaces` | NMF latent factor object; `$Z` holds the scores matrix |
| `$Cells_niches` | Cell niche enrichment per latent factor |
| `$TME_states` | Factor-to-meta-program mapping table |
| `$Metaprograms_reference` | Reference NES matrix used for meta-program mapping |

``` r

# Latent factor scores (samples x factors) — primary input for downstream analyses
head(res$Latent_spaces$Z)

# Cell group composition per group
res$Cell_groups[[2]]

# TME state annotation per factor
head(res$TME_states)
```

## Key parameters

| Parameter | Default | Description |
|----|----|----|
| `task` | `"unsupervised"` | `"supervised"` or `"unsupervised"` |
| `contrast` | `NULL` | Column in `coldata` for supervised DEG contrast |
| `ref_level` | `NULL` | Reference level for the contrast |
| `normalized` | `TRUE` | Whether counts are already log-normalised |
| `deconv_methods` | all 5 | Deconvolution algorithms to run |
| `TF.collection` | `"CollecTRI"` | TF regulon: `"CollecTRI"`, `"Dorothea"`, `"ARACNE"` |
| `cancer_type` | `NULL` | For meta-program mapping: `"skcm"`, `"blca"`, `"luad"` |
| `corr` | `0.7` | Correlation threshold for deconvolution subgrouping |
| `corr_mod` | `0.9` | Correlation threshold for TF module merging |
| `pval` | `0.05` | P-value threshold for cell group construction |
| `batch` | `FALSE` | Enable multi-cohort mode |
| `batch_id` | `NULL` | Column in `coldata` identifying cohorts |
| `return` | `TRUE` | Return the result object |
