# Summary Statistics for Features

This function allows you to 'describe' 'omics features using the
describe() function from the psych package, as well as estimate
variance, a dispersion index, the coeficent of variation, and shapiro's
W-statistic. The output from
[`psych::describe()`](https://rdrr.io/pkg/psych/man/describe.html)
includes feature-level skewness and kurtosis estimates.

## Usage

``` r
feature_describe(data)
```

## Arguments

- data:

  matrix, the 'omics data matrix. Samples in row, features in columns

## Value

a data frame of summary statistics for features (columns) of a matrix
