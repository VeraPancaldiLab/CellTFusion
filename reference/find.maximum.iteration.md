# Find maximum iteration from cell subgroups

Scans a nested list of cell subgroups and returns the highest iteration
number found across all elements.

## Usage

``` r
find.maximum.iteration(cells.groups)
```

## Arguments

- cells.groups:

  A nested list of cell subgroups. Names of inner elements must follow
  the pattern `*.Iteration.<n>` where `<n>` is an integer.

## Value

An integer giving the maximum iteration number across all subgroups.
