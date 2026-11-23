# Feature Summary Statistics

This function estimates feature statistics for samples in a matrix of
'omics features.

## Usage

``` r
feature_summary(
  omiprep,
  source_layer = "input",
  outlier_udist = 5,
  tree_cut_height = 0.5,
  feature_selection = "max_var_exp",
  sample_ids = NULL,
  feature_ids = NULL,
  features_exclude = NULL,
  output = "data.frame",
  cores = NULL,
  fast = FALSE
)

## S7 method for class <omiprep::Omiprep>
feature_summary(
  omiprep,
  source_layer = "input",
  outlier_udist = 5,
  tree_cut_height = 0.5,
  feature_selection = "max_var_exp",
  sample_ids = NULL,
  feature_ids = NULL,
  features_exclude = NULL,
  output = "data.frame",
  cores = NULL,
  fast = FALSE
)
```

## Arguments

- omiprep:

  an object of class Omiprep

- source_layer:

  character, the data layer to summarise

- outlier_udist:

  the unit distance in SD or IQR from the mean or median estimate,
  respectively outliers are identified at. Default value is 5.

- tree_cut_height:

  numeric, the threshold for feature independence in hierarchical
  clustering. Default is 0.5.

- feature_selection:

  character, either 'max_var_exp' or 'least_missingness', how to select
  the independent feature within clusters

- sample_ids:

  character, vector of sample ids to work with

- feature_ids:

  character, vector of feature ids to work with

- features_exclude:

  character, vector of feature id indicating features to exclude from
  the sample and PCA summary analysis but keep in the data

- output:

  character, type of output, either 'object' to return the updated
  metaboprep object, or 'data.frame' to return the data.

- cores:

  number of cores available for parallelism; the default null will try
  find the maximum available cores - 1; set to 1 for linear, but
  potentially slow, computation of the correlation matrix.

- fast:

  If `TRUE`, accelerates correlation computation by imputing missing
  values to the column minimum, pre-ranking all columns, and computing
  Pearson correlation on ranked data (approximating Spearman).
  Substantially faster than exact Spearman at large feature dimensions
  (\\p \> 5000\\) but assumes values are missing because they fall below
  a lower limit of detection: every missing value is treated as the
  lowest value of its feature. If values are missing for other reasons
  (e.g. at random) the correlations are weakened and fewer features are
  clustered, so use the exact method (`fast = FALSE`). Features with
  high missingness will have many tied ranks at the minimum (ensure
  these are filtered out appropriately with the missingness option).
  Default `FALSE`.
