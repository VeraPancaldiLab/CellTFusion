# Classify samples by high or low deconvolution values in given cell groups

Classifies samples in `coldata` as "High" or "Low" based on whether
their deconvolution values in all specified cell groups exceed the
median.

## Usage

``` r
classify.deconvolution(coldata, deconvolution, group)
```

## Arguments

- coldata:

  A data frame with sample metadata.

- deconvolution:

  A numeric matrix/data frame with cell deconvolution features (samples
  x cell groups).

- group:

  A character vector specifying the cell groups (columns) to consider.

## Value

The input `coldata` with an additional factor column `Cells_level`
("High" or "Low").
