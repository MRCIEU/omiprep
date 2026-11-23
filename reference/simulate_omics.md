# Simulate an omics dataset with known quality problems

Generates a synthetic samples-by-features matrix on a log-abundance
scale, with correlated feature clusters, a chosen pattern of missing
values and, optionally, planted quality problems whose identity is
returned alongside the data. Because the truth is known, the output can
be used to check what
[`quality_control()`](https://mrcieu.github.io/omiprep/reference/quality_control.md)
removes and keeps under a given set of thresholds, to compare settings
(for example `fast = TRUE` against `fast = FALSE`), and as example data
for tests and tutorials.

## Usage

``` r
simulate_omics(
  n,
  p,
  n_factors = NULL,
  frac_clustered = 0.3,
  loading_range = c(0.6, 0.9),
  miss_rate = 0.1,
  miss_mech = c("lod", "mcar", "mixed", "none"),
  complete_frac = 0,
  plant = NULL,
  seed = NULL,
  equal_rates = FALSE,
  miss_ignores_shift = FALSE,
  as_omiprep = FALSE
)
```

## Arguments

- n:

  integer, the number of samples (rows).

- p:

  integer, the number of features (columns).

- n_factors:

  integer, the number of latent factors behind the correlated feature
  clusters. `NULL` (default) uses one factor per five clustered
  features.

- frac_clustered:

  numeric (0 to 1), the share of features that load on a latent factor.
  The remaining features are independent.

- loading_range:

  numeric vector of length 2, the range from which each clustered
  feature's loading on its factor is drawn. Values must lie between 0
  and 1.

- miss_rate:

  numeric (0 to \<1), the target share of missing values, as the mean of
  the per-feature rates.

- miss_mech:

  character, how values go missing: `"lod"` (default), the lowest values
  of each feature, as when they fall below a lower limit of detection;
  `"mcar"`, values missing completely at random; `"mixed"`, half of each
  feature's missing values by each mechanism; `"none"`, no missing
  values.

- complete_frac:

  numeric (0 to 1), the share of features that are left with no missing
  values from `miss_mech`.

- plant:

  named list of problems to plant, or `NULL` (default) for none. See the
  section "Planted problems".

- seed:

  a single integer, or `NULL` (default). If supplied,
  [`set.seed()`](https://rdrr.io/r/base/Random.html) is called with it
  so the same call returns the same data. This changes the random number
  state of the session.

- equal_rates:

  logical. `FALSE` (default): per-feature missingness rates vary around
  `miss_rate` (beta-distributed, at most 0.75), as in real data. Some
  features labelled `"clean"` will then exceed a missingness threshold
  by themselves (about 1 in 8 features is at least 20% missing when
  `miss_rate = 0.1`). `TRUE`: every feature has exactly `miss_rate`, so
  only planted features have high missingness.

- miss_ignores_shift:

  logical. `FALSE` (default): values below the detection limit are
  chosen after the dilution and outlier shifts, so a sample shifted down
  also gains missing values. `TRUE`: they are chosen on the unshifted
  values, so diluted and outlier samples have the same missingness as
  other samples and differ only in their values.

- as_omiprep:

  logical. `FALSE` (default) returns a list; `TRUE` returns an
  [`Omiprep()`](https://mrcieu.github.io/omiprep/reference/class_omiprep.md)
  object built from it.

## Value

If `as_omiprep = FALSE`, a list with three elements that can be passed
directly to
[`Omiprep()`](https://mrcieu.github.io/omiprep/reference/class_omiprep.md):

- `data`: numeric matrix, samples in rows and features in columns, with
  row and column names.

- `samples`: data frame with columns `sample_id`, `batch` (a random
  label, `"b1"` or `"b2"`, unrelated to the data) and `truth`, the label
  of what was planted in the sample (`"clean"` if nothing).

- `features`: data frame with columns `feature_id`, `platform`
  (constant), `truth`, the label of what was planted in the feature
  (`"clean"` if nothing), and `cluster_factor`, the latent factor a
  clustered feature loads on (`NA` for independent features).

If `as_omiprep = TRUE`, an `Omiprep` object holding the same data, with
the tables above in its `samples` and `features` properties.

## Details

The data are built in this order:

1.  **Base data.** Every value starts as an independent standard normal
    draw. The first `frac_clustered` of the features are then assigned
    to one of `n_factors` latent factors and mixed with it (loading
    drawn uniformly from `loading_range`), which gives clusters of
    correlated features. Each feature is then given its own mean
    (uniform on 8 to 14) and standard deviation (uniform on 0.4 to 1.5).

2.  **Planted problems and biological controls** (see the `plant`
    argument), applied to the complete data. Each sample or feature
    receives at most one of them.

3.  **Missing values**, according to `miss_mech` and `miss_rate`.

4.  **Extra missing values** for planted poor features and failed
    samples.

Nothing is planted unless requested: with `plant = NULL` every sample
and feature is labelled `"clean"`. The label of what was planted is
returned in the `truth` column of the sample and feature tables, so it
stays with the data when they are passed to
[`Omiprep()`](https://mrcieu.github.io/omiprep/reference/class_omiprep.md)
and through
[`quality_control()`](https://mrcieu.github.io/omiprep/reference/quality_control.md).
Note that the QC report treats `truth` like any other sample or feature
variable, so it will appear among the candidate batch variables.

## Planted problems (`plant`)

`plant` is a named list. Each `*_frac` entry is the share (0 to 1) of
samples or features given that problem; all default to 0. Entries not
supplied keep their default.

Technical problems, which quality control is expected to remove:

- `dilution_frac`, `dilution_shift` (default 1): diluted or concentrated
  samples. Every feature of the sample is shifted up or down (direction
  chosen at random per sample) by `dilution_shift` feature standard
  deviations. Truth label `"tech_dilution"`.

- `outlier_frac`, `outlier_shift` (default 1): outlier samples. Each is
  shifted along one random direction across the features, shared by the
  planted samples and not present in the rest of the data, scaled by
  `outlier_shift` feature standard deviations. Truth label
  `"tech_outlier"`.

- `failed_frac`, `failed_miss` (default `NULL`): failed samples, given
  extra values missing at random. `failed_miss` is the share of features
  set missing; if `NULL` it is drawn per sample from 0.5 to 0.9. Truth
  label `"tech_failed"`.

- `poor_feat_frac`, `poor_miss` (default `NULL`), `poor_mech` (default
  `"random"`): poorly detected features, given extra missing values.
  `poor_miss` is the share of samples set missing; if `NULL` it is drawn
  per feature from 0.5 to 0.95. With `poor_mech = "random"` the extra
  values are missing at random; with `"lod"` the lowest values still
  observed are set missing, as for a feature that is often below its
  detection limit. Truth label `"tech_poor_feature"`.

- `trunc_frac`, `trunc_quantile` (default 0.6): features capped at an
  upper limit. Values above the `trunc_quantile` quantile of the feature
  are set to that quantile, which gives a left-skewed feature. Truth
  label `"tech_truncated"`.

Biological controls, which quality control is expected to keep:

- `bio_subgroup_frac`, `bio_subgroup_shift` (default 2): a subgroup of
  samples shifted by `bio_subgroup_shift` feature standard deviations on
  the features of the first latent factor only. Requires
  `frac_clustered > 0`. Truth label `"bio_subgroup"`.

- `bio_tsa_frac`, `bio_tsa_shift` (default 0.5): a subgroup of samples
  with higher abundance, shifted by `bio_tsa_shift` feature standard
  deviations on a random 40% of the features. Truth label `"bio_tsa"`.

- `rare_frac`: features with a rare effect. A random 8% of samples are
  shifted up by 3 feature standard deviations, which gives a
  right-skewed feature. Truth label `"bio_rare_effect"`.

## Examples

``` r
# 200 samples by 200 features, 10% of each feature missing, with a few failed samples and
# poorly detected features planted
sim <- simulate_omics(
  n = 200, p = 200, miss_rate = 0.10, equal_rates = TRUE, seed = 1,
  plant = list(failed_frac = 0.02, poor_feat_frac = 0.05)
)
dim(sim$data)
#> [1] 200 200
table(sim$samples$truth)
#> 
#>       clean tech_failed 
#>         196           4 
table(sim$features$truth)
#> 
#>             clean tech_poor_feature 
#>               190                10 

# the same data as an Omiprep object; the truth labels are kept in the sample and feature tables
m <- simulate_omics(
  n = 200, p = 200, miss_rate = 0.10, equal_rates = TRUE, seed = 1, as_omiprep = TRUE,
  plant = list(failed_frac = 0.02, poor_feat_frac = 0.05)
)
table(m@samples$truth)
#> 
#>       clean tech_failed 
#>         196           4 

# run quality control and compare what was excluded with what was planted
if (FALSE) { # \dontrun{
m <- quality_control(m, source_layer = "input", cores = 1)
excluded_samples  <- unique(unlist(m@exclusions$samples))
excluded_features <- unique(unlist(m@exclusions$features))
table(truth = m@samples$truth,  excluded = m@samples$sample_id   %in% excluded_samples)
table(truth = m@features$truth, excluded = m@features$feature_id %in% excluded_features)
} # }
```
