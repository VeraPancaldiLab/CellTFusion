# TME State Characterisation

``` r

library(CellTFusion)
#> 
#> 
```

Once latent factors have been extracted (see [Cell Group
Construction](https://verapancaldilab.github.io/CellTFusion/articles/02-cell-groups.md)),
each factor can be functionally annotated by linking it to known
biological programs. This article covers:

1.  **Hallmark GSEA** — associate each latent factor with MSigDB
    Hallmark gene sets
2.  **Meta-program mapping** — compare factors to cancer-type-specific
    reference programs derived from TCGA
3.  **TME subtype annotation** — annotate factors with established TME
    immune subtypes

These steps require:

- `counts.norm` — log-normalized expression matrix (genes × samples)
- `latent_spaces` — output of
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md);
  the `$Z` element is the samples × factors matrix

------------------------------------------------------------------------

## Hallmark GSEA per latent factor

[`compute_factor_gsea()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_factor_gsea.md)
fits a multivariate `limma` model using the latent factor scores as
covariates, ranks genes by moderated t-statistic per factor, and runs a
pre-ranked Hallmark gene set enrichment analysis via `fgsea`
([Korotkevich et al. 2021](#ref-Korotkevich2021)) against the MSigDB
Hallmark collection ([Liberzon et al. 2015](#ref-Liberzon2015)). The
Hallmark collection summarizes ~50 curated gene sets representing
well-defined, non-redundant biological processes (e.g. EMT, interferon
response, hypoxia, inflammatory response), making it a natural,
interpretable reference to functionally annotate a data-driven latent
factor. A dot plot of the top enriched pathways is saved per factor.

``` r

gsea_results <- compute_factor_gsea(
  RNA.tpm     = counts.norm,        # genes × samples expression matrix
  features_df = latent_spaces$Z,    # samples × factors (from compute.latent_factors)
  plot_dot    = TRUE,
  top_n       = 10,
  file_name   = "Tutorial",
  width       = 8,
  height      = 10
)
```

The result is a list with two elements:

| Element         | Description                                            |
|-----------------|--------------------------------------------------------|
| `$DE_results`   | Named list of limma DEG tables, one per factor         |
| `$GSEA_results` | Named list of fgsea result data frames, one per factor |

------------------------------------------------------------------------

## Map factors to cancer-type meta-programs

### Background: what are meta-programs?

A **meta-program (MP)** is a recurrent module of co-expressed genes
reflecting a specific transcriptional cell state (e.g. cell cycle,
hypoxia, epithelial-mesenchymal transition, interferon response, stress)
that recurs *across* tumors, patients, and even cancer types, rather
than being specific to one dataset. This concept was systematically
characterized by Gavish et al. ([Gavish et al. 2023](#ref-Gavish2023))
(“Hallmarks of transcriptional intratumour heterogeneity across a
thousand tumours”, *Nature*, 2023), who derived a curated, pan-cancer
set of intratumour heterogeneity meta-programs from single-cell RNA-seq
data spanning many cancer types. CellTFusion builds cancer-type-specific
reference meta-programs from bulk TCGA RNA-seq data using the same
underlying logic — recurrent, biologically interpretable transcriptional
programs — so that latent factors identified in an independent study can
be related back to a well-established vocabulary of TME states instead
of being described only in dataset-specific terms.

### How `map_factors_to_metaprograms()` works

1.  **Reference construction (offline, per cancer type):** for a given
    TCGA cohort (e.g. `"skcm"` for melanoma, `"blca"` for bladder
    cancer), a Hallmark GSEA is run per meta-program to build a
    reference Hallmarks x meta-programs normalized enrichment score
    (NES) matrix. This reference is shipped with the package as
    pre-built `.RData` objects, one per supported cancer type.
2.  **Study-side profile:** the Hallmark NES profile of each *study*
    latent factor, computed in the previous step by
    [`compute_factor_gsea()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_factor_gsea.md),
    is assembled into a Hallmarks x factors matrix
    ([`build_nes_matrix()`](https://verapancaldilab.github.io/CellTFusion/reference/build_nes_matrix.md)).
3.  **Matching:** for each study factor, the function computes, for
    every reference meta-program, the mean NES of the Hallmark gene sets
    that define that meta-program. The meta-program with the highest
    mean NES (and NES \> 0) is assigned as the factor’s `best_MP`.

This effectively asks: *“among all known recurrent TME/tumour
transcriptional programs for this cancer type, which one does this
latent factor’s Hallmark signature resemble most?”* — turning an
abstract NMF factor into a biologically named TME state (e.g. “hypoxia”,
“interferon response”, “stromal/CAF-like”).

Supported `cancer_type` values (bundled TCGA reference): `"skcm"`
(melanoma), `"blca"` (bladder cancer), `"luad"` (lung adenocarcinoma).

``` r

mp_mapping <- map_factors_to_metaprograms(
  gsea_study  = gsea_results,
  cancer_type = "skcm",        # match to your cancer type
  plot        = TRUE,
  file_name   = "Tutorial"
)
```

The result contains:

| Element | Description |
|----|----|
| `$factor_mapping` | Data frame with, per study factor: `factor`, `best_MP` (closest reference meta-program), `best_score` (its mean NES), `all_scores` (mean NES against every reference meta-program) |
| `$reference` | The reference Hallmarks x meta-programs NES matrix used for comparison |

A heatmap of factor-to-meta-program scores is saved to `Results/` when
`plot = TRUE`.

------------------------------------------------------------------------

## Annotate factors with TME immune subtypes

[`map_factors_to_TME()`](https://verapancaldilab.github.io/CellTFusion/reference/map_factors_to_TME.md)
uses TCGA sample-level TME subtype annotations to score each latent
factor for enrichment in specific immune environments (e.g.,
immune-desert, immune-excluded, inflamed).

``` r

tme_annotation <- map_factors_to_TME(
  cancer_name = "skcm",
  Z           = latent_spaces$Z,
  plot        = TRUE,
  file_name   = "Tutorial"
)
```

------------------------------------------------------------------------

## Derive meta-programs from GSEA (unsupervised)

If you have run CellTFusion across multiple cohorts and want to derive
consensus meta-programs from scratch rather than mapping to TCGA
references, use
[`derive_meta_programs()`](https://verapancaldilab.github.io/CellTFusion/reference/derive_meta_programs.md).
This clusters factors by their Hallmark NES profiles to identify
recurring biological programs.

``` r

meta_programs <- derive_meta_programs(
  gsea_results = gsea_results,
  k            = NULL,       # number of clusters; NULL for automatic selection
  file_name    = "Tutorial",
  plot         = TRUE
)
```

------------------------------------------------------------------------

## Putting it all together

In practice, this full annotation sequence is run automatically when
using the
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
wrapper. The outputs are accessible directly from the result object:

``` r

res <- CellTFusion(raw.counts = raw.counts, cancer_type = "skcm", ...)

# Access GSEA + meta-program mapping results
res$TME_states           # factor-to-meta-program mapping table
res$Metaprograms_reference  # reference NES matrix

# Access latent factors used as input to GSEA
head(res$Latent_spaces$Z)
```

## References

Gavish, Avishai, Michael Tyler, Alissa C. Greenwald, et al. 2023.
“Hallmarks of Transcriptional Intratumour Heterogeneity Across a
Thousand Tumours.” *Nature* 618 (7965): 598–606.
<https://doi.org/10.1038/s41586-023-06130-4>.

Korotkevich, Gennady, Vladimir Sukhov, Nikolay Budin, Boris Shpak, Maxim
N. Artyomov, and Alexey Sergushichev. 2021. “Fast Gene Set Enrichment
Analysis.” *bioRxiv*, ahead of print. <https://doi.org/10.1101/060012>.

Liberzon, Arthur, Chet Birger, Helga Thorvaldsdóttir, Mahmoud Ghandi,
Jill P. Mesirov, and Pablo Tamayo. 2015. “The Molecular Signatures
Database (MSigDB) Hallmark Gene Set Collection.” *Cell Systems* 1 (6):
417–25. <https://doi.org/10.1016/j.cels.2015.12.004>.
