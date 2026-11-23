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

  numeric, optional skewness threshold. If \`NULL\` (the default), only
  skewness is returned and no exclusion flag is calculated.

- direction:

  character, direction of skewness to flag. One of \`"left"\` (the
  default), \`"right"\`, or \`"both"\`.

## Value

data.frame with columns \`feature_id\`, \`skew\`, and
\`exclude_by_skewness\` (logical; \`NA\` if \`threshold = NULL\`).

## Details

Skewness is calculated for each feature across samples using
[`psych::describe()`](https://rdrr.io/pkg/psych/man/describe.html),
which by default returns the type 3 estimator \\b_1 = g_1 \cdot
((n-1)/n)^{3/2}\\, where \\g_1\\ is the standard moment coefficient of
skewness and \\n\\ is the number of non-missing values for that feature.

Missing values are dropped feature-wise, so \\n\\ varies between
features and no imputation is performed. No minimum sample-size is
imposed: a feature with only a handful of non-missing values will still
return an estimate, and such estimates are unstable. Features with zero
variance return \`NaN\` and are never flagged.

The threshold is a cutoff on the skewness statistic itself, not a
p-value, and is applied to the **signed** value unless \`direction =
"both"\`:

- \`direction = "left"\` flags features where \`skew \<= -threshold\`;

- \`direction = "right"\` flags features where \`skew \>= threshold\`;

- \`direction = "both"\` flags features where \`abs(skew) \>=
  threshold\`.

Skewness is computed on the values as supplied. It is not invariant to
transformation, so a feature's skewness on a raw intensity scale and on
a log scale can differ substantially.

A skewed distribution is a description of shape, not evidence of
technical failure. Genuine biological variation, subgroup structure and
low-frequency phenomena all produce skewed features. Treat the flag as a
prompt to inspect a feature rather than as grounds to remove it; see
\`vignette("skewness_qc")\` for guidance on interpretation.

## See also

[`quality_control()`](https://mrcieu.github.io/omiprep/reference/quality_control.md)
for applying this filter within the QC pipeline, and
[`feature_describe()`](https://mrcieu.github.io/omiprep/reference/feature_describe.md),
which also reports per-feature skewness.
