# Omiprep Object

A \`Omiprep\` object is a container for matrices of 'omics data, along
with associated metadata. It allows for efficient storage and
manipulation of data, supporting quality control, transformations, and
various analyses. This object facilitates easy access to data layers,
sample and feature summaries, outlier treatment, and more.

An \`Omiprep\` object can be built directly from any data source, not
only the commercial formats handled by the \`read\_\*()\` functions. It
requires a numeric matrix with samples in rows and features in columns,
plus two metadata data.frames containing the columns \`sample_id\` and
\`feature_id\` respectively. These two column names are matched exactly
and are case-sensitive; no alternative spellings are mapped. Ids are
coerced to character so that numeric identifiers read from flat files
match the matrix dimnames. See \`vignette("custom_data")\` for a worked
example starting from plain text files.

## Usage

``` r
Omiprep(
  data,
  samples,
  features,
  exclusions = list(samples = list(user_excluded = character(),
    extreme_sample_missingness = character(), user_defined_sample_missingness =
    character(), user_defined_sample_totalpeakarea = character(),
    user_defined_sample_pca_outlier = character()), features = list(user_excluded =
    character(), extreme_feature_missingness = character(),
    user_defined_feature_missingness = character(), user_defined_feature_skewness =
    character())),
  feature_summary = array(data = NA_real_, dim = c(0, 0, 0)),
  sample_summary = array(data = NA_real_, dim = c(0, 0, 0))
)
```

## Arguments

- data:

  numeric matrix, the data matrix containing 'omics values (not to be
  set directly).

- samples:

  data.frame, a data frame containing sample-related information (not to
  be set directly).

- features:

  data.frame, a data frame containing feature-related information (not
  to be set directly).

- exclusions:

  list, holds exclusion codes for data masking (not to be set directly).

- feature_summary:

  numeric matrix, summary statistics for features (not to be set
  directly).

- sample_summary:

  numeric matrix, summary statistics for samples (not to be set
  directly).

## Value

An object of class Omiprep, an S7 class.

## Slots

- `data`:

  numeric matrix, the 'omics data.

- `samples`:

  data.frame, the samples data frame

- `features`:

  data.frame, the features data frame

- `exclusions`:

  list, exclusion codes (mask for data).

- `feature_summary`:

  numeric matrix, feature summary statistics.

- `sample_summary`:

  numeric matrix, sample summary statistics.
