# Import Olink Proteomic Data

## Import Olink data

Read in the Olink data using the `read_olink` function. Here we will
read in the example data provided with the package, as a list object.

``` r
library(omiprep)

# example file
olink_file <- system.file("extdata", "olink_v1_example.txt", package="omiprep")

# import
dat <- read_olink(olink_file, 
                  return_Omiprep = FALSE )
```

## Quick look to identify the types of data imported

Note that we would advocate that you read in olink data as a list object
first, as it will often have additional information that you may wish to
inspect prior to creating the `Omiprep` object.

``` r
names(dat)
#> [1] "data"             "samples"          "features"         "controls"        
#> [5] "control_metadata"
```

## Create Omiprep object

Once imported, we pass the data to the Omiprep() function to build the
`Omiprep` class object.

``` r
mydata <- Omiprep(data      = dat$data, 
                  features  = dat$features, 
                  samples   = dat$samples)
```

## Quick summary of the Omiprep object

``` r
summary(mydata)
#> Omiprep Object Summary
#> --------------------------
#> Samples      : 90
#> Features     : 100
#> Data Layers  : 1
#> Layer Names  : input
#> 
#> Sample Summary Layers : none
#> Feature Summary Layers: none
#> 
#> Sample Annotation (metadata):
#>   Columns: 3
#>   Names  : sample_id, PlateID, QC_Warning
#> 
#> Feature Annotation (metadata):
#>   Columns: 8
#>   Names  : feature_id, UniProt, Assay, MissingFreq, Panel, Panel_Version, LOD, Normalization
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

## QC Olink data

Perform the QC steps using the `quality_control` function.

``` r
mydata <- mydata |>
  quality_control(source_layer        = "input", 
                  sample_missingness  = 0.2, 
                  feature_missingness = 0.2, 
                  total_peak_area_sd  = 5, 
                  outlier_udist       = 5, 
                  outlier_treatment   = "leave_be", 
                  winsorize_quantile  = 1.0, 
                  tree_cut_height     = 0.5, 
                  pc_outlier_sd       = 5,
                  feature_selection   = "max_var_exp",
                  features_exclude_but_keep = NULL,
                  cores               = 1
                  )
#> 
#> ── Starting Omics QC Process ───────────────────────────────────────────────────
#> ℹ Validating input parameters
#> 
#> ℹ Validating input parameters── Starting 'Omics QC Process ──────────────────────────────────────────────────
#> ℹ Validating input parameters✔ Validating input parameters [18ms]
#> 
#> ℹ Validating input parameters
#> ✔ Validating input parameters [14ms]
#> 
#> ℹ Sample & Feature Summary Statistics for raw data
#> AF =  17
#> ✔ Sample & Feature Summary Statistics for raw data [580ms]
#> 
#> ℹ Copying input data to new 'qc' data layer
#> ✔ Copying input data to new 'qc' data layer [23ms]
#> 
#> ℹ Assessing for extreme sample missingness >=80% - excluding 0 sample(s)
#> ✔ Assessing for extreme sample missingness >=80% - excluding 0 sample(s) [19ms]
#> 
#> ℹ Assessing for extreme feature missingness >=80% - excluding 0 feature(s)
#> ✔ Assessing for extreme feature missingness >=80% - excluding 0 feature(s) [23m…
#> 
#> ℹ Assessing for sample missingness at specified level of >=20% - excluding 0 sa…
#> ✔ Assessing for sample missingness at specified level of >=20% - excluding 0 sa…
#> 
#> ℹ Assessing for feature missingness at specified level of >=20% - excluding 0 f…
#> ✔ Assessing for feature missingness at specified level of >=20% - excluding 0 f…
#> 
#> ℹ Calculating total peak abundance outliers at +/- 5 Sdev - excluding 0 sample(…
#> ✔ Calculating total peak abundance outliers at +/- 5 Sdev - excluding 0 sample(…
#> 
#> ℹ Running sample data PCA outlier analysis at +/- 5 Sdev
#> ✔ Running sample data PCA outlier analysis at +/- 5 Sdev [17ms]
#> 
#> ℹ Sample PCA outlier analysis - re-identify feature independence and PC outlier…
#> AF =  17
#> ✔ Sample PCA outlier analysis - re-identify feature independence and PC outlier…
#> 
#> ℹ Creating final QC dataset...
#> AF =  17
#> 
#> ℹ Creating final QC dataset...── Step timings ──
#> ℹ Creating final QC dataset...
#> ℹ Creating final QC dataset...
#>                         step seconds  pct
#>                   validation    0.02  1.0
#>                summarise_raw    0.55 27.8
#>                   copy_layer    0.00  0.0
#>   extreme_sample_missingness    0.00  0.0
#>  extreme_feature_missingness    0.00  0.0
#>           sample_missingness    0.00  0.0
#>              total_peak_area    0.00  0.0
#>                summarise_pca    0.63 31.8
#>              summarise_final    0.56 28.3
#>                        total    1.98 99.9
#> ✔ Creating final QC dataset... [598ms]
#> 
#> ℹ 'Omics QC Process Completed
#> ✔ 'Omics QC Process Completed [14ms]
```

## Quick summary of the Omiprep object following QC

``` r
summary(mydata)
#> Omiprep Object Summary
#> --------------------------
#> Samples      : 90
#> Features     : 100
#> Data Layers  : 2
#> Layer Names  : input, qc
#> 
#> Sample Summary Layers : input, qc
#> Feature Summary Layers: input, qc
#> 
#> Sample Annotation (metadata):
#>   Columns: 5
#>   Names  : sample_id, PlateID, QC_Warning, reason_excluded, excluded
#> 
#> Feature Annotation (metadata):
#>   Columns: 10
#>   Names  : feature_id, UniProt, Assay, MissingFreq, Panel, Panel_Version, LOD, Normalization, reason_excluded, excluded
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
