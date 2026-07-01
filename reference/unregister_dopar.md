# Unregister a parallel backend registered with doParallel

Switches the foreach backend back to sequential execution and calls
[`gc()`](https://rdrr.io/r/base/gc.html) to release memory held by the
parallel workers.

## Usage

``` r
unregister_dopar()
```

## Value

Called for its side effect; returns `NULL` invisibly.
