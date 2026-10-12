# Recursively retrieve the original deconvolution features of a subgroup

Replaces a deconvolution subgroup by its member features, going down
through nested subgroups, and returns the original deconvolution feature
names (e.g. `B.cells_Subgroup.1` gives `CBSX_BPRNACan_B.cells`,
`CBSX_BPRNACan3DProMet_B.cells`, ...). A name that is not a subgroup is
returned unchanged. The cell types are read from the returned names
afterwards with
[`extract_cells()`](https://verapancaldilab.github.io/CellTFusion/reference/extract_cells.md).

## Usage

``` r
get_all_cells(subgroup_name, cell_subgroups)
```

## Arguments

- subgroup_name:

  Character. A single deconvolution feature name: a subgroup or an
  original feature.

- cell_subgroups:

  A named list mapping subgroup names to vectors of member names (which
  can themselves be subgroup names), e.g. the flattened
  `Deconvolution subgroups composition` of
  [`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html).

## Value

A character vector of the unique original deconvolution feature names
contained in `subgroup_name` (the name itself if it is not a subgroup).
