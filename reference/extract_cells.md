# Extract cell types from deconvolution feature names

Reads the cell type of each deconvolution feature from its name, with
the cell type vocabulary of multideconv
([`multideconv::get_cell_type_nomenclature()`](https://verapancaldilab.github.io/multideconv/reference/get_cell_type_nomenclature.html),
plus `cells_extra`). The first vocabulary name found in the feature name
is kept, and a `_Subgroup.<n>` suffix (or `.Iteration.<n>` of
multideconv \< 0.2.0) is removed. Vocabulary names are matched literally
(a "." only matches a ".").

## Usage

``` r
extract_cells(features, cells_extra = NULL)
```

## Arguments

- features:

  A character vector of deconvolution feature names: raw features
  (`<method>_<signature>_<cell type>`, or `<method>_<cell type>`) or
  subgroups (`<cell type>_Subgroup.<n>`).

- cells_extra:

  Optional character vector of additional cell type names, added to the
  vocabulary.

## Value

A character vector of the distinct cell types found (`character(0)` if
none). Feature names without any vocabulary cell type are ignored.

## Examples

``` r
extract_cells(c("Epidish_BPRNACan_Monocytes", "Monocytes_Subgroup.2", "Quantiseq_B.cells"))
#> [1] "Monocytes" "B.cells"  
extract_cells("Method_Signature_Mesenchymal.cells", cells_extra = "Mesenchymal.cells")
#> [1] "Mesenchymal.cells"
```
