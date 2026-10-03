# Example Deconvolution Results

A toy dataset included with the package to illustrate usage of functions
in **CellTFusion**. This dataset contains the output of a deconvolution
run (cell-type proportions estimated with different methods and
signatures) on the tutorial samples.

## Usage

``` r
deconv.tuto
```

## Format

A data frame with 192 samples as rows and 399 deconvolution features as
columns (named `<method>_<signature>_<cell type>`, or
`<method>_<cell type>` for methods with a single signature such as
Quantiseq)

## Examples

``` r
data(deconv.tuto)
```
