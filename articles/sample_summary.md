# Sample summary

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
#> ℹ Sample & Feature Summary Statistics for raw data✔ Sample & Feature Summary Statistics for raw data [231ms]
#> 
#> ℹ Copying input data to new 'qc' data layer
#> ✔ Copying input data to new 'qc' data layer [28ms]
#> 
#> ℹ Assessing for extreme sample missingness >=80% - excluding 0 sample(s)
#> ✔ Assessing for extreme sample missingness >=80% - excluding 0 sample(s) [21ms]
#> 
#> ℹ Assessing for extreme feature missingness >=80% - excluding 0 feature(s)
#> ✔ Assessing for extreme feature missingness >=80% - excluding 0 feature(s) [20m…
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
#> ✔ Running sample data PCA outlier analysis at +/- 5 Sdev [19ms]
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
#>                summarise_raw    0.22  28.3
#>                   copy_layer    0.00   0.0
#>   extreme_sample_missingness    0.00   0.0
#>  extreme_feature_missingness    0.00   0.0
#>           sample_missingness    0.00   0.0
#>          total_sum_abundance    0.00   0.0
#>                summarise_pca    0.23  29.6
#>              summarise_final    0.11  14.1
#>                        total    0.78 100.3
#> ✔ Creating final QC dataset... [141ms]
#> 
#> ℹ 'Omics QC Process Completed
#> ✔ 'Omics QC Process Completed [17ms]
```

## Sample Summary

### View sample summary from the QC pipeline

``` r

# Note: the quality_control() ultimately returns the sample_summary attribute as a matrix.
df <- as.data.frame(mydata@sample_summary[1:3, 1:6, "input"])
df <- cbind(sample_id = rownames(df), df)

df |> knitr::kable( digits = 3, row.names = FALSE, align = "c") |>
  kableExtra::kable_styling(full_width = TRUE) 
```

| sample_id | missingness | non_na_feature_count | tsa_total | tsa_complete_features | complete_feature_count | outlier_count |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| id_100 | 0 | 20 | 38.699 | 38.699 | 20 | 0 |
| id_99 | 0 | 20 | 39.451 | 39.451 | 20 | 0 |
| id_98 | 0 | 20 | 45.135 | 45.135 | 20 | 0 |

### Manually run sample summary

While sample summary is run as a part of the quality_control() function
pipeline you can run the function yourself, on any layer you wish.

``` r

# NOTE:
# outlier_udist = number of IQRs from the median at which a value is flagged.
# 1.0 here is illustrative; in practice we favour 5.0, which is the default value
# for the quality_control() function.
sample_sum1 <- sample_summary(omiprep         = mydata, 
                              source_layer    = "input", 
                              outlier_udist   = 1.0,
                              output          = "data.frame")
```

#### Table of sample summary

``` r

sample_sum1 |> 
  head(n = 10) |> 
  knitr::kable( digits = 3, row.names = FALSE, align = "c") |>
  kableExtra::kable_styling(full_width = TRUE) 
```

| sample_id | missingness | non_na_feature_count | tsa_total | tsa_complete_features | complete_feature_count | outlier_count |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| id_100 | 0 | 20 | 38.699 | 38.699 | 20 | 0 |
| id_99 | 0 | 20 | 39.451 | 39.451 | 20 | 0 |
| id_98 | 0 | 20 | 45.135 | 45.135 | 20 | 1 |
| id_97 | 0 | 20 | 37.289 | 37.289 | 20 | 2 |
| id_96 | 0 | 20 | 31.545 | 31.545 | 20 | 3 |
| id_95 | 0 | 20 | 36.751 | 36.751 | 20 | 2 |
| id_94 | 0 | 20 | 37.384 | 37.384 | 20 | 0 |
| id_93 | 0 | 20 | 34.069 | 34.069 | 20 | 0 |
| id_92 | 0 | 20 | 38.942 | 38.942 | 20 | 0 |
| id_91 | 0 | 20 | 42.743 | 42.743 | 20 | 2 |

### Run sample summary on subset

Using the `sample_ids` and `feature_ids` arguments you can run the
summary for a subset of the data. Note: all rows will be return, however
summary data will only be returned for the specified ids.

``` r

## define a vector of sample IDs
sids <- mydata@samples[mydata@samples$sex == "female", "sample_id"] 

## define a vector of feature IDs
## extract only those features run on `pos`itive ion mode. 
fids <- mydata@features[mydata@features$platform == "pos", "feature_id"] 

# run sample summary on subset
sample_sum_subset <- sample_summary(omiprep       = mydata, 
                                    source_layer  = "input", 
                                    outlier_udist = 1.0,
                                    sample_ids    = sids,
                                    feature_ids   = fids,
                                    output        = "data.frame")
```

#### Table of sample summary on subset

``` r

sample_sum_subset |> 
  na.omit() |>
  head(n = 10) |> 
  knitr::kable( digits = 3, row.names = FALSE, align = "c") |>
  kableExtra::kable_styling(full_width = TRUE) 
```

| sample_id | missingness | non_na_feature_count | tsa_total | tsa_complete_features | complete_feature_count | outlier_count |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| id_98 | 0 | 4 | 6.590 | 6.590 | 4 | 0 |
| id_97 | 0 | 4 | 4.484 | 4.484 | 4 | 0 |
| id_96 | 0 | 4 | 7.318 | 7.318 | 4 | 0 |
| id_94 | 0 | 4 | 7.641 | 7.641 | 4 | 0 |
| id_93 | 0 | 4 | 7.053 | 7.053 | 4 | 0 |
| id_92 | 0 | 4 | 8.123 | 8.123 | 4 | 0 |
| id_90 | 0 | 4 | 4.924 | 4.924 | 4 | 0 |
| id_89 | 0 | 4 | 9.921 | 9.921 | 4 | 0 |
| id_88 | 0 | 4 | 6.454 | 6.454 | 4 | 1 |
| id_87 | 0 | 4 | 7.254 | 7.254 | 4 | 0 |

## Principal Componet Analysis

### View PCs from the QC pipeline

PCs and outliers are available as a part of the `quality_countrol()`
function pipeline.

``` r

# Note: the quality_control() ultimately returns the sample_summary attribute as a matrix.
df <- as.data.frame(mydata@sample_summary[1:3, -c(1:6), "input"])
df <- cbind(sample_id = rownames(df), df)

df |> knitr::kable( digits = 3, row.names = FALSE, align = "c") |>
  kableExtra::kable_styling(full_width = TRUE) 
```

| sample_id | pc1 | pc2 | pc3 | pc4 | pc5 | pc6 | pc7 | pc8 | pc9 | pc10 | pc1_3_sd_outlier | pc2_3_sd_outlier | pc3_3_sd_outlier | pc1_4_sd_outlier | pc2_4_sd_outlier | pc3_4_sd_outlier | pc1_5_sd_outlier | pc2_5_sd_outlier | pc3_5_sd_outlier |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| id_100 | 1.087 | 0.321 | 0.003 | 0.477 | 0.837 | 0.256 | -0.162 | -0.276 | 0.577 | -1.466 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| id_99 | 0.418 | 0.573 | 1.581 | 0.110 | -1.154 | -0.256 | -0.103 | 0.115 | 1.961 | -1.201 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| id_98 | -2.249 | 0.131 | 0.040 | -0.605 | -0.720 | 0.763 | -0.614 | 0.076 | -0.521 | -1.457 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

### Manually run PCA analysis

You can derive PCs and identify outlier independent of the
[`quality_control()`](https://mrcieu.github.io/omiprep/reference/quality_control.md)
function, however.

[`pc_and_outliers()`](https://mrcieu.github.io/omiprep/reference/pc_and_outliers.md)
performs principal component analysis. Missing data is imputed to the
median and used to identify the number of informative or ‘significant’
PCs by (1) an acceleration analysis, and (2) a parallel analysis.
Finally the number of sample outliers are determined at 3, 4, and 5
standard deviations from the mean on the top PCs as determined by the
acceleration factor analysis.

``` r

pc_analysis <- pc_and_outliers(omiprep      = mydata, 
                               source_layer = "input",
                               sample_ids   = sids, ## It is also possible to run on a subset of samples and/or features
                               feature_ids  = NULL
                               )
#> ℹ Number of informative PCs (Scree acceleration factor): 3
```

#### Table of PCA analysis results

Returned are the PC eigenvectors for the top 10 PCs, and outlier counts
at 3, 4, and 5 standard deviations from the mean for the top two PCs

``` r

pc_analysis |> 
  head(n = 10) |> 
  knitr::kable( digits = 3, row.names = FALSE, align = "c") |>
  kableExtra::kable_styling(full_width = TRUE) 
```

| sample_id | pc1 | pc2 | pc3 | pc4 | pc5 | pc6 | pc7 | pc8 | pc9 | pc10 | pc1_3_sd_outlier | pc2_3_sd_outlier | pc3_3_sd_outlier | pc1_4_sd_outlier | pc2_4_sd_outlier | pc3_4_sd_outlier | pc1_5_sd_outlier | pc2_5_sd_outlier | pc3_5_sd_outlier |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| id_98 | 1.504 | -0.899 | -0.072 | -0.987 | 1.221 | 0.429 | -1.603 | 0.586 | -1.086 | -0.193 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| id_97 | 0.826 | 0.458 | 1.623 | 0.188 | -1.195 | 1.666 | -1.318 | -2.481 | -1.081 | 1.007 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| id_96 | -2.067 | 0.785 | -1.024 | 1.720 | 1.435 | 0.663 | 1.057 | -0.073 | -1.554 | -0.479 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| id_94 | -1.635 | 1.103 | 0.791 | 1.214 | -0.920 | -0.986 | 1.635 | 0.401 | -0.521 | 1.214 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| id_93 | -1.399 | 1.202 | -0.168 | -2.277 | 1.068 | -1.482 | -0.129 | -1.497 | 0.834 | 1.243 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| id_92 | 2.299 | 0.640 | 1.312 | 0.308 | -0.094 | 0.492 | 0.595 | -0.751 | 0.994 | -1.484 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| id_90 | -0.723 | 1.214 | 1.251 | 1.834 | -1.257 | 0.548 | -1.184 | -1.282 | 0.829 | 0.076 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| id_89 | -0.292 | 0.097 | 2.914 | 1.079 | -0.734 | -1.654 | 0.601 | 0.015 | -0.721 | 0.914 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| id_88 | 2.664 | 1.542 | -2.028 | 0.080 | 0.419 | -0.416 | 1.305 | -1.710 | 0.352 | -0.213 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| id_87 | -2.053 | 0.898 | -0.883 | 0.960 | 0.578 | 1.764 | -0.703 | 0.158 | 0.812 | 0.708 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

#### PCA plot

``` r

library(ggplot2)
## extract variance explained attribute
varexp <- attr(pc_analysis, 'input_varexp')
## extract PC data
pcs = pc_analysis[, c("pc1","pc2")]
## generate plot
pcs |> ggplot(aes(x = pc1, y = pc2)) +
  geom_point(size = 2, col = "#377EB8") +
  labs(x = paste0("PC1; VarExp =  ", varexp[1]*100, "%"), 
       y = paste0("PC2; VarExp =  ", varexp[2]*100, "%")) +
  theme_bw()
```

![](sample_summary_files/figure-html/pca_plot-1.png)

### Additional pc_and_outliers() attributes

In addition, the variance explained vector is appended to the returned
`data.frame` as and `attribute`. This can be accessed with the attribute
name: `[source_layer]_varexp`, in this case we used the `input` data,
therefore the attribute name is `input_varexp`. In a similar way, the
results of the acceleration analysis (`input_num_pcs_scree`) and a
parallel analysis (`input_num_pcs_parallel`) can also be extracted.

``` r

library(ggplot2)

# extract varexp from attributes
varexp <- attr(pc_analysis, 'input_varexp')

# subset to top 100 for nicer plotting
if (length(varexp) > 100) varexp <- varexp[1:100]

# get acceleration and parallel analysis results
af <- attr(pc_analysis, 'input_num_pcs_scree')

# as data.frame
x_labs <- sub("(?i)pc","", names(varexp))
ve     <- data.frame("pc"      = factor(x_labs, levels=x_labs),
                     "var_exp" = varexp)
lines  <- data.frame("Analysis" = c("Acceleration"), 
                     "pc"       = c(af))   

# plot
ggplot(ve, aes(x = pc, y = var_exp)) +
  geom_line(color = "grey") +
  geom_point(shape = 21, fill = "#377EB8", size = 2) +
  geom_vline(data = lines, aes(xintercept = pc, color = Analysis), inherit.aes = FALSE) +
  scale_color_manual(values = c("Acceleration"="#E41A1C")) +
  scale_x_discrete(labels = function(x) ifelse(seq_along(x) %% 10 == 0 | x==1, x, "")) +
  labs(x = "PC", y = "Variance explained") +
  theme_classic() +
  theme(legend.position = "top")
```

![Variance explained](sample_summary_files/figure-html/ScreePlot-1.png)

## Run sample & feature summaries together

``` r

sf_sum <- summarise(omiprep         = mydata, 
                    source_layer    = "input", 
                    outlier_udist   = 1.0,
                    tree_cut_height = 0.5,
                    sample_ids      = sids, ## It is also possible to run on a subset of samples and/or features
                    feature_ids     = NULL,
                    output          = "data.frame", 
                    cores           = 1)
#> ℹ Number of informative PCs (Scree acceleration factor): 3

## two data frames are returned as a list object
names(sf_sum)
#> [1] "sample_summary"  "feature_summary"
```

### Table of sample summary on subset

Note that when the summarise() function is used the sample summary now
includes PCA derived summary data. This is not the case when the
sample_summary() function is run alone, as seen above. The reason for
the difference is because the PCA data is dependent upon the
feature_summary() analysis.

Also, please note that when running on a subset, you are returned the
full summary for all samples and features, but only the summary data for
the specified subset will be populated, the rest will be `NA`.

| sample_id | missingness | non_na_feature_count | tsa_total | tsa_complete_features | complete_feature_count | outlier_count | pc1 | pc2 | pc3 | pc4 | pc5 | pc6 | pc7 | pc8 | pc9 | pc10 | pc1_3_sd_outlier | pc2_3_sd_outlier | pc3_3_sd_outlier | pc1_4_sd_outlier | pc2_4_sd_outlier | pc3_4_sd_outlier | pc1_5_sd_outlier | pc2_5_sd_outlier | pc3_5_sd_outlier |
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| id_100 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| id_99 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| id_98 | 0 | 20 | 45.392 | 45.392 | 20 | 1 | 1.504 | -0.899 | -0.072 | -0.987 | 1.221 | 0.429 | -1.603 | 0.586 | -1.086 | -0.193 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| id_97 | 0 | 20 | 37.473 | 37.473 | 20 | 3 | 0.826 | 0.458 | 1.623 | 0.188 | -1.195 | 1.666 | -1.318 | -2.481 | -1.081 | 1.007 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| id_96 | 0 | 20 | 32.042 | 32.042 | 20 | 3 | -2.067 | 0.785 | -1.024 | 1.720 | 1.435 | 0.663 | 1.057 | -0.073 | -1.554 | -0.479 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| id_95 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| id_94 | 0 | 20 | 37.615 | 37.615 | 20 | 2 | -1.635 | 1.103 | 0.791 | 1.214 | -0.920 | -0.986 | 1.635 | 0.401 | -0.521 | 1.214 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| id_93 | 0 | 20 | 34.238 | 34.238 | 20 | 1 | -1.399 | 1.202 | -0.168 | -2.277 | 1.068 | -1.482 | -0.129 | -1.497 | 0.834 | 1.243 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| id_92 | 0 | 20 | 39.200 | 39.200 | 20 | 1 | 2.299 | 0.640 | 1.312 | 0.308 | -0.094 | 0.492 | 0.595 | -0.751 | 0.994 | -1.484 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| id_91 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |
