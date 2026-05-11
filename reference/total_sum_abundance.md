# Estimates total sum abundance

This function estimates total sum abundance for numeric data in a
matrix, for (1) all features and (2) all features with complete data.

## Usage

``` r
total_sum_abundance(data, ztransform = TRUE)
```

## Arguments

- data:

  matrix, the 'omics data matrix. Samples in rows, features in columns

- ztransform:

  logical, should the feature data be z-transformed and absolute value
  minimum, mean shifted prior to summing the feature values. TRUE or
  FALSE.

## Value

a data frame of estimates for (1) total sum abundance and (2) total sum
abundance at complete features for each samples
