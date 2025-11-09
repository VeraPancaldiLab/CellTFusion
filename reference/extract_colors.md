# Extract colors

Extract TF module colors from cell type group names

## Usage

``` r
extract_colors(module_colors, cell_group_name)
```

## Arguments

- module_colors:

  Character vector of TF module colors, e.g., from compute.WTCNA()

- cell_group_name:

  Character vector or string of cell type group names to search within

## Value

Character vector of matched module colors found in cell_group_name,
ordered by their appearance. Returns NA if no matches are found.
