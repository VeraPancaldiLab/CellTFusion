# Create TFs modules

This function re-create existing TF modules on a different TF activity
matrix.

## Usage

``` r
create_tfs_modules(TF.matrix, network_tfs)
```

## Arguments

- TF.matrix:

  TFs activity matrix with samples as rows and TFs as columns.

- network_tfs:

  A TF network object obtained from compute.WTCNA() from which TF
  modules need to be re-create.

## Value

A matrix of TFs modules scores across samples
