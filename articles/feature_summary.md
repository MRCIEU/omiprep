# Feature summary

## Create Omiprep object

``` r

library(omiprep)

# import data
data     <- read.csv(system.file("extdata", "dummy_data.csv",     package = "omiprep"), header=T, row.names = 1) |> as.matrix()
samples  <- read.csv(system.file("extdata", "dummy_samples.csv",  package = "omiprep"), header=T, row.names = 1)
features <- read.csv(system.file("extdata", "dummy_features.csv", package = "omiprep"), header=T, row.names = 1)

# create object
mydata <- Omiprep(data = data, samples = samples, features = features)
```

## Summary of Omiprep object

``` r

summary(mydata)
#> Omiprep Object Summary
#> --------------------------
#> Samples      : 100
#> Features     : 20
#> Data Layers  : 1
#> Layer Names  : input
#> 
#> Sample Summary Layers : none
#> Feature Summary Layers: none
#> 
#> Sample Annotation (metadata):
#>   Columns: 5
#>   Names  : sample_id, age, sex, pos, neg
#> 
#> Feature Annotation (metadata):
#>   Columns: 5
#>   Names  : feature_id, platform, pathway, derived_feature, xenobiotic_feature
#> 
#> Exclusion Codes Summary:
#> 
#>   Sample Exclusions:
#> Exclusion | Count
#> -----------------
#> user_excluded                     | 0
#> extreme_sample_missingness        | 0
#> user_defined_sample_missingness   | 0
#> user_defined_sample_totalpeakarea | 0
#> user_defined_sample_pca_outlier   | 0
#> 
#>   Feature Exclusions:
#> Exclusion | Count
#> -----------------
#> user_excluded                    | 0
#> extreme_feature_missingness      | 0
#> user_defined_feature_missingness | 0
#> user_defined_feature_skewness    | 0
```

## Run standard quality control

``` r

mydata = quality_control(mydata)
#> 
#> ── Starting Omics QC Process ───────────────────────────────────────────────────
#> ℹ Validating input parameters
#> ✔ Validating input parameters [8ms]
#> 
#> ℹ Sample & Feature Summary Statistics for raw data
#> ℹ Number of informative PCs (Scree acceleration factor): 3
#> ℹ Sample & Feature Summary Statistics for raw data✔ Sample & Feature Summary Statistics for raw data [406ms]
#> 
#> ℹ Copying input data to new 'qc' data layer
#> ✔ Copying input data to new 'qc' data layer [24ms]
#> 
#> ℹ Assessing for extreme sample missingness >=80% - excluding 0 sample(s)
#> ✔ Assessing for extreme sample missingness >=80% - excluding 0 sample(s) [17ms]
#> 
#> ℹ Assessing for extreme feature missingness >=80% - excluding 0 feature(s)
#> ✔ Assessing for extreme feature missingness >=80% - excluding 0 feature(s) [17m…
#> 
#> ℹ Assessing for sample missingness at specified level of >=20% - excluding 0 sa…
#> ✔ Assessing for sample missingness at specified level of >=20% - excluding 0 sa…
#> 
#> ℹ Assessing for feature missingness at specified level of >=20% - excluding 0 f…
#> ✔ Assessing for feature missingness at specified level of >=20% - excluding 0 f…
#> 
#> ℹ Calculating total sum abundance outliers at +/- 5 Sdev - excluding 0 sample(s)
#> ✔ Calculating total sum abundance outliers at +/- 5 Sdev - excluding 0 sample(s…
#> 
#> ℹ Running sample data PCA outlier analysis at +/- 5 Sdev
#> ✔ Running sample data PCA outlier analysis at +/- 5 Sdev [23ms]
#> 
#> ℹ Sample PCA outlier analysis - re-identify feature independence and PC outlier…
#> ℹ Number of informative PCs (Scree acceleration factor): 3
#> ℹ Sample PCA outlier analysis - re-identify feature independence and PC outlier…! The stated max PCs [max_num_pcs=10] to use in PCA outlier assessment is greater than the number of available informative PCs [3]
#> ℹ Sample PCA outlier analysis - re-identify feature independence and PC outlier…✔ Sample PCA outlier analysis - re-identify feature independence and PC outlier…
#> 
#> ℹ Creating final QC dataset...
#> ℹ Number of informative PCs (Scree acceleration factor): 3
#> ℹ Creating final QC dataset...
#> ℹ Creating final QC dataset...── Step timings ──
#> ℹ Creating final QC dataset...
#> ℹ Creating final QC dataset...
#>                         step seconds   pct
#>                   validation    0.00   0.0
#>                summarise_raw    0.39  42.1
#>                   copy_layer    0.00   0.0
#>   extreme_sample_missingness    0.00   0.0
#>  extreme_feature_missingness    0.00   0.0
#>           sample_missingness    0.00   0.0
#>          total_sum_abundance    0.00   0.0
#>                summarise_pca    0.22  23.7
#>              summarise_final    0.11  11.9
#>                        total    0.93 100.3
#> ✔ Creating final QC dataset... [142ms]
#> 
#> ℹ 'Omics QC Process Completed
#> ✔ 'Omics QC Process Completed [17ms]
```

## Feature Summary

### View feature summary from the QC pipeline

``` r

# Note: the quality_control() ultimately returns the feature_summary attribute as a matrix.
df <- t( as.data.frame(mydata@feature_summary[, 1:5, "input"]) )
df <- as.data.frame( round( df , 3) )
df <- cbind(feature_id = rownames(df), df)

df |> knitr::kable( digits = 3, row.names = FALSE, align = "c") |>
  kableExtra::kable_styling(full_width = TRUE) 
```

| feature_id | missingness | outlier_count | n | mean | sd | median | min | max | range | skew | kurtosis | se | missing | var | disp_index | coef_variance | W | log10_W | k | independent_features |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| metab_id_1 | 0 | 0 | 100 | 0.511 | 0.293 | 0.530 | 0.000 | 0.993 | 0.992 | -0.123 | -1.231 | 0.029 | 0 | 0.086 | 0.168 | 0.574 | 0.949 | 0.744 | 1 | 1 |
| metab_id_2 | 0 | 0 | 100 | 0.521 | 0.310 | 0.547 | 0.018 | 0.993 | 0.975 | -0.150 | -1.404 | 0.031 | 0 | 0.096 | 0.184 | 0.594 | 0.924 | 0.834 | 2 | 1 |
| metab_id_3 | 0 | 0 | 100 | 0.488 | 0.283 | 0.504 | 0.001 | 0.995 | 0.994 | -0.036 | -1.109 | 0.028 | 0 | 0.080 | 0.165 | 0.580 | 0.963 | 0.749 | 3 | 1 |
| metab_id_4 | 0 | 0 | 100 | 0.464 | 0.286 | 0.466 | 0.004 | 0.992 | 0.988 | 0.092 | -1.199 | 0.029 | 0 | 0.082 | 0.177 | 0.617 | 0.954 | 0.833 | 4 | 1 |
| metab_id_5 | 0 | 0 | 100 | 0.521 | 0.293 | 0.547 | 0.004 | 0.976 | 0.972 | -0.219 | -1.161 | 0.029 | 0 | 0.086 | 0.164 | 0.561 | 0.945 | 0.782 | 5 | 1 |

### Manually run feature summary

While feature summary is run as a part of the quality_control() function
pipeline you can run the function yourself, on any layer you wish.

``` r

# NOTE:
# outlier_udist = number of IQRs from the median at which a value is flagged.
# 1.0 here is illustrative; in practice we favour 5.0, which is the default value
# for the quality_control() function.
feature_sum1 <- feature_summary(omiprep         = mydata, 
                                source_layer    = "input", 
                                outlier_udist   = 1.0,
                                tree_cut_height = 0.5,
                                output          = "data.frame", 
                                cores           = 1)
```

#### Table of feature summary

``` r

feature_sum1 |> 
  head(n = 10) |> 
  knitr::kable( digits = 3, row.names = FALSE, align = "c") |>
  kableExtra::kable_styling(full_width = FALSE) 
```

| feature_id | missingness | outlier_count | n | mean | sd | median | min | max | range | skew | kurtosis | se | missing | var | disp_index | coef_variance | W | log10_W | k | independent_features |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| metab_id_1 | 0 | 5 | 100 | 0.511 | 0.293 | 0.530 | 0.000 | 0.993 | 0.992 | -0.123 | -1.231 | 0.029 | 0 | 0.086 | 0.168 | 0.574 | 0.949 | 0.744 | 1 | TRUE |
| metab_id_2 | 0 | 0 | 100 | 0.521 | 0.310 | 0.547 | 0.018 | 0.993 | 0.975 | -0.150 | -1.404 | 0.031 | 0 | 0.096 | 0.184 | 0.594 | 0.924 | 0.834 | 2 | TRUE |
| metab_id_3 | 0 | 10 | 100 | 0.488 | 0.283 | 0.504 | 0.001 | 0.995 | 0.994 | -0.036 | -1.109 | 0.028 | 0 | 0.080 | 0.165 | 0.580 | 0.963 | 0.749 | 3 | TRUE |
| metab_id_4 | 0 | 5 | 100 | 0.464 | 0.286 | 0.466 | 0.004 | 0.992 | 0.988 | 0.092 | -1.199 | 0.029 | 0 | 0.082 | 0.177 | 0.617 | 0.954 | 0.833 | 4 | TRUE |
| metab_id_5 | 0 | 11 | 100 | 0.521 | 0.293 | 0.547 | 0.004 | 0.976 | 0.972 | -0.219 | -1.161 | 0.029 | 0 | 0.086 | 0.164 | 0.561 | 0.945 | 0.782 | 5 | TRUE |
| metab_id_6 | 0 | 7 | 100 | 0.490 | 0.259 | 0.473 | 0.007 | 0.993 | 0.986 | 0.007 | -1.006 | 0.026 | 0 | 0.067 | 0.137 | 0.528 | 0.973 | 0.803 | 6 | TRUE |
| metab_id_7 | 0 | 7 | 100 | 0.479 | 0.277 | 0.441 | 0.029 | 0.992 | 0.963 | 0.135 | -1.211 | 0.028 | 0 | 0.077 | 0.160 | 0.579 | 0.953 | 0.899 | 7 | TRUE |
| metab_id_8 | 0 | 0 | 100 | 0.476 | 0.312 | 0.491 | 0.001 | 0.999 | 0.998 | 0.059 | -1.350 | 0.031 | 0 | 0.097 | 0.205 | 0.656 | 0.936 | 0.796 | 8 | TRUE |
| metab_id_9 | 0 | 10 | 100 | 0.468 | 0.260 | 0.489 | 0.005 | 0.975 | 0.971 | 0.000 | -1.090 | 0.026 | 0 | 0.068 | 0.144 | 0.556 | 0.968 | 0.800 | 9 | TRUE |
| metab_id_10 | 0 | 0 | 100 | 0.524 | 0.290 | 0.532 | 0.019 | 0.993 | 0.974 | -0.158 | -1.252 | 0.029 | 0 | 0.084 | 0.161 | 0.554 | 0.945 | 0.841 | 10 | TRUE |

### Run feature summary on subset

Using the `sample_ids` and `feature_ids` arguments you can run the
summary for a subset of the data. Note: all rows will be return, however
summary data will only be returned for the specified ids.

``` r

## define a vector of sample IDs
sids <- mydata@samples[mydata@samples$sex == "female", "sample_id"] 

## define a vector of feature IDs
fids <- mydata@features[, "feature_id"] |> sample(10)

# NOTE:
# outlier_udist = number of IQRs from the median at which a value is flagged.
# 1.0 here is illustrative; in practice we favour 5.0, which is the default value
# for the quality_control() function.
feature_sum_subset <- feature_summary(omiprep         = mydata, 
                                      source_layer    = "input", 
                                      outlier_udist   = 1.0,
                                      tree_cut_height = 0.5,
                                      sample_ids      = sids,
                                      feature_ids     = fids,
                                      output          = "data.frame",
                                      cores           = 1)
```

#### Table of feature summary for subset

``` r

feature_sum_subset |> 
  na.omit() |>
  knitr::kable( digits = 3, row.names = FALSE, align = "c") |>
  kableExtra::kable_styling(full_width = FALSE) |>
  kableExtra::scroll_box(width = "100%", height = "500px")
```

| feature_id | missingness | outlier_count | n | mean | sd | median | min | max | range | skew | kurtosis | se | missing | var | disp_index | coef_variance | W | log10_W | k | independent_features |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| metab_id_4 | 0 | 5 | 55 | 0.498 | 0.273 | 0.478 | 0.009 | 0.992 | 0.983 | 0.094 | -1.159 | 0.037 | 0 | 0.075 | 0.150 | 0.549 | 0.960 | 0.836 | 7 | TRUE |
| metab_id_5 | 0 | 9 | 55 | 0.580 | 0.296 | 0.632 | 0.013 | 0.976 | 0.963 | -0.436 | -1.035 | 0.040 | 0 | 0.087 | 0.151 | 0.510 | 0.929 | 0.734 | 3 | TRUE |
| metab_id_6 | 0 | 6 | 55 | 0.520 | 0.284 | 0.560 | 0.007 | 0.993 | 0.986 | -0.273 | -1.103 | 0.038 | 0 | 0.080 | 0.155 | 0.546 | 0.952 | 0.757 | 4 | TRUE |
| metab_id_9 | 0 | 6 | 55 | 0.461 | 0.270 | 0.489 | 0.005 | 0.975 | 0.971 | 0.012 | -1.151 | 0.036 | 0 | 0.073 | 0.158 | 0.585 | 0.963 | 0.798 | 9 | TRUE |
| metab_id_12 | 0 | 0 | 55 | 0.518 | 0.296 | 0.498 | 0.006 | 0.989 | 0.983 | -0.125 | -1.233 | 0.040 | 0 | 0.088 | 0.170 | 0.573 | 0.946 | 0.801 | 2 | TRUE |
| metab_id_13 | 0 | 4 | 55 | 0.544 | 0.313 | 0.595 | 0.005 | 0.990 | 0.984 | -0.265 | -1.325 | 0.042 | 0 | 0.098 | 0.181 | 0.576 | 0.927 | 0.752 | 1 | TRUE |
| metab_id_15 | 0 | 7 | 55 | 0.495 | 0.280 | 0.522 | 0.032 | 0.991 | 0.959 | 0.054 | -1.139 | 0.038 | 0 | 0.078 | 0.158 | 0.566 | 0.961 | 0.874 | 5 | TRUE |
| metab_id_16 | 0 | 3 | 55 | 0.463 | 0.315 | 0.447 | 0.005 | 0.998 | 0.992 | 0.276 | -1.384 | 0.042 | 0 | 0.099 | 0.214 | 0.679 | 0.909 | 0.868 | 6 | TRUE |
| metab_id_17 | 0 | 3 | 55 | 0.449 | 0.287 | 0.437 | 0.025 | 0.972 | 0.947 | 0.185 | -1.280 | 0.039 | 0 | 0.083 | 0.184 | 0.639 | 0.938 | 0.902 | 8 | TRUE |
| metab_id_18 | 0 | 1 | 55 | 0.536 | 0.316 | 0.599 | 0.028 | 0.999 | 0.971 | -0.073 | -1.480 | 0.043 | 0 | 0.100 | 0.186 | 0.590 | 0.918 | 0.867 | 10 | TRUE |

## Additional feature_summary() attributes

``` r

## The attributes include column names, row names, and class for the feature summary table
## as well as a hierarchical cluster dendrogram or `input_tree` and the parameter values for 
## outlier_udist and input_tree_cut_height passed to the function. 
names( attributes(feature_sum1) )
#> [1] "names"                 "row.names"             "class"                
#> [4] "input_tree"            "input_outlier_udist"   "input_tree_cut_height"
```

### hierarchical cluster dendrogram

In addition to the summary data, the hierarchical cluster dendrogram is
appended to the returned `data.frame` as and `attribute`. This can be
accessed with the attribute name: `[source_layer]_tree`, in this case we
summarised the `input` data, therefore the attribute name is
`input_tree`.

``` r

suppressPackageStartupMessages(library(dendextend))

# extract tree from attributes
tree <- attr(feature_sum1, 'input_tree')
dend <- stats::as.dendrogram(tree)

# color the independent features blue
metab_color       <- feature_sum1[, c("feature_id", "independent_features")]
metab_color       <- metab_color[match(labels(dend), metab_color$feature_id), ]
metab_color$color <- ifelse(metab_color$independent_features==TRUE, "#477EB8", "grey")

# format dendrogram for ploting
dend <- dend |>
  dendextend::set("labels_cex", 0.75) |>
  dendextend::set("labels_col", metab_color$color) |>
  dendextend::set("branches_lwd", 1) |>
  dendextend::set("branches_k_color",  value = metab_color$color)

## plot the dendrogram
dend |> plot(main = "Feature clustering dendrogram")
abline(h = 0.5, col = "#E41A1C", lwd = 1.5)
```

![Decision tree showing feature importance in
dataset](feature_summary_files/figure-html/tree-1.png)

## Run sample & feature summaries together

``` r

# NOTE:
# outlier_udist = number of IQRs from the median at which a value is flagged.
# 1.0 here is illustrative; in practice we favour 5.0, which is the default value
# for the quality_control() function.
sf_sum <- summarise(omiprep         = mydata, 
                    source_layer    = "input", 
                    outlier_udist   = 1.0,
                    tree_cut_height = 0.5,
                    sample_ids      = sids, ## It is also possible to run on a subset of samples and/or features
                    feature_ids     = fids,
                    output          = "data.frame", 
                    cores           = 1)
#> ℹ Number of informative PCs (Scree acceleration factor): 2
```

## Table of feature summary from summarise() function

``` r

sf_sum$feature_summary |> 
  na.omit() |>
  knitr::kable( digits = 3, row.names = FALSE, align = "c") |>
  kableExtra::kable_styling(full_width = TRUE) |>
  kableExtra::scroll_box(width = "100%", height = "500px")
```

| feature_id | missingness | outlier_count | n | mean | sd | median | min | max | range | skew | kurtosis | se | missing | var | disp_index | coef_variance | W | log10_W | k | independent_features |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| metab_id_4 | 0 | 5 | 55 | 0.498 | 0.273 | 0.478 | 0.009 | 0.992 | 0.983 | 0.094 | -1.159 | 0.037 | 0 | 0.075 | 0.150 | 0.549 | 0.960 | 0.836 | 7 | TRUE |
| metab_id_5 | 0 | 9 | 55 | 0.580 | 0.296 | 0.632 | 0.013 | 0.976 | 0.963 | -0.436 | -1.035 | 0.040 | 0 | 0.087 | 0.151 | 0.510 | 0.929 | 0.734 | 3 | TRUE |
| metab_id_6 | 0 | 6 | 55 | 0.520 | 0.284 | 0.560 | 0.007 | 0.993 | 0.986 | -0.273 | -1.103 | 0.038 | 0 | 0.080 | 0.155 | 0.546 | 0.952 | 0.757 | 4 | TRUE |
| metab_id_9 | 0 | 6 | 55 | 0.461 | 0.270 | 0.489 | 0.005 | 0.975 | 0.971 | 0.012 | -1.151 | 0.036 | 0 | 0.073 | 0.158 | 0.585 | 0.963 | 0.798 | 9 | TRUE |
| metab_id_12 | 0 | 0 | 55 | 0.518 | 0.296 | 0.498 | 0.006 | 0.989 | 0.983 | -0.125 | -1.233 | 0.040 | 0 | 0.088 | 0.170 | 0.573 | 0.946 | 0.801 | 2 | TRUE |
| metab_id_13 | 0 | 4 | 55 | 0.544 | 0.313 | 0.595 | 0.005 | 0.990 | 0.984 | -0.265 | -1.325 | 0.042 | 0 | 0.098 | 0.181 | 0.576 | 0.927 | 0.752 | 1 | TRUE |
| metab_id_15 | 0 | 7 | 55 | 0.495 | 0.280 | 0.522 | 0.032 | 0.991 | 0.959 | 0.054 | -1.139 | 0.038 | 0 | 0.078 | 0.158 | 0.566 | 0.961 | 0.874 | 5 | TRUE |
| metab_id_16 | 0 | 3 | 55 | 0.463 | 0.315 | 0.447 | 0.005 | 0.998 | 0.992 | 0.276 | -1.384 | 0.042 | 0 | 0.099 | 0.214 | 0.679 | 0.909 | 0.868 | 6 | TRUE |
| metab_id_17 | 0 | 3 | 55 | 0.449 | 0.287 | 0.437 | 0.025 | 0.972 | 0.947 | 0.185 | -1.280 | 0.039 | 0 | 0.083 | 0.184 | 0.639 | 0.938 | 0.902 | 8 | TRUE |
| metab_id_18 | 0 | 1 | 55 | 0.536 | 0.316 | 0.599 | 0.028 | 0.999 | 0.971 | -0.073 | -1.480 | 0.043 | 0 | 0.100 | 0.186 | 0.590 | 0.918 | 0.867 | 10 | TRUE |
