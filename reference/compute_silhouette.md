# Compute mean silhouette width for a clustering

Compute mean silhouette width for a clustering

## Usage

``` r
compute_silhouette(clusters, distance_matrix)
```

## Arguments

- clusters:

  An integer vector of cluster assignments (one per sample).

- distance_matrix:

  A numeric matrix used to compute Euclidean distances between samples.

## Value

A single numeric value: the mean silhouette width across all samples.
