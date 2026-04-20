# Identify Independent Features in a Numeric Matrix

This function identifies independent features using Spearman's rho
correlation distances, and a dendrogram tree cut step.

## Usage

``` r
tree_and_independent_features(
  data,
  tree_cut_height = 0.5,
  features_exclude = NULL,
  feature_selection = "max_var_exp",
  cores = NULL,
  fast = FALSE
)
```

## Arguments

- data:

  matrix, the 'omics data matrix. samples in row, features in columns

- tree_cut_height:

  the tree cut height. A value of 0.2 (1-Spearman's rho) is equivalent
  to saying that features with a rho \>= 0.8 are NOT independent.

- features_exclude:

  character, vector of feature id indicating features to exclude from
  the sample and PCA summary analysis but keep in the data

- feature_selection:

  character. Method for selecting a representative feature from each
  correlated feature cluster.

- cores:

  number of cores available for parallelism; the default null will try
  find the maximum available cores - 1; set to 1 for linear, but
  potentially slow, computation of the correlation matrix.

- fast:

  If `TRUE`, accelerates correlation computation by imputing missing
  values to the column minimum, pre-ranking all columns, and computing
  Pearson correlation on ranked data (approximating Spearman).
  Substantially faster than exact Spearman at large feature dimensions
  (\\p \> 5000\\) but assumes missing data are missing at random.
  Features with high missingness will have inflated rank ties at the
  median (ensure these are filtered out appropriately with the
  missingness option). Default `FALSE`. One of:

  `"max_var_exp"`

  :   (Default) Selects the feature with the highest sum of absolute
      Spearman correlations to other features in the cluster;
      effectively the feature explaining the most shared variance.

  `"least_missingness"`

  :   Selects the feature with the fewest missing values within the
      cluster.

## Value

A list with the following components:

- data:

  A \`data.frame\` with:

  - \`feature_id\`: Feature (column) names from the input matrix.

  - \`k\`: The cluster index assigned to each feature after tree
    cutting.

  - \`independent_features\`: Logical indicator of whether the feature
    was selected as an independent (representative) feature.

- tree:

  A \`hclust\` object representing the hierarchical clustering of the
  features based on 1 - \|Spearman's rho\| distance.
