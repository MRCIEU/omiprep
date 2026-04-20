# Estimate Feature Skewness

Estimate per-feature skewness and optionally flag features beyond a
user-defined skewness threshold.

## Usage

``` r
feature_skewness(data, threshold = NULL, direction = "left")
```

## Arguments

- data:

  matrix, a numeric matrix with samples in rows and features in columns.

- threshold:

  numeric, optional skewness threshold. If \`NULL\`, only skewness is
  returned and no exclusion flag is calculated.

- direction:

  character, direction of skewness to flag. One of \`"left"\`,
  \`"right"\`, or \`"both"\`.

## Value

data.frame with columns \`feature_id\`, \`skew\`, and
\`exclude_by_skewness\` (logical; \`NA\` if \`threshold = NULL\`).
