# Recursively retrieve all base cell types for a subgroup

Walks a nested subgroup hierarchy and returns the leaf-level cell type
names that belong to the requested subgroup.

## Usage

``` r
get_all_cells(subgroup_name, cell_subgroups)
```

## Arguments

- subgroup_name:

  Character. Name of the subgroup or base cell type to resolve.

- cell_subgroups:

  A named list mapping subgroup names to vectors of member names (which
  can themselves be subgroup keys).

## Value

A character vector of unique base cell type names belonging to
`subgroup_name`.
