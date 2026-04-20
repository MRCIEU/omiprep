# Read Nightingale Data (format 1)

Read Nightingale Data (format 1)

## Usage

``` r
read_nightingale(filepath, return_Omiprep = TRUE)
```

## Arguments

- filepath:

  character, commercial Nightingale excel sheet with extension .xls or
  .xlsx

- return_Omiprep:

  logical, if TRUE (default) return a Omiprep object, if FALSE return a
  list.

## Value

list or Omiprep object, list(data = matrix, samples = samples
data.frame, features = features data.frame)

## Examples

``` r
# version 1 data format
filepath1 <- system.file("extdata", "nightingale_v1_example.xlsx", package = "omiprep")
m <- read_nightingale(filepath1)

# version 2 data format
filepath2 <- system.file("extdata", "nightingale_v2_example.xlsx", package = "omiprep")
m <- read_nightingale(filepath2)
```
