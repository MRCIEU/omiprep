# Batch Normalise

## Create Omiprep object

``` r

library(omiprep)

# example file
filepath <- system.file("extdata", "metabolon_v1.1_example.xlsx", package = "omiprep")

## Read in Metabolon v1.1 example and directly create Omiprep object
mydata <- read_metabolon(filepath, sheet = "OrigScale" , return_Omiprep = TRUE )


## summary of Omiprep object
summary(mydata)
#> Omiprep Object Summary
#> --------------------------
#> Samples      : 125
#> Features     : 253
#> Data Layers  : 1
#> Layer Names  : input
#> 
#> Sample Summary Layers : none
#> Feature Summary Layers: none
#> 
#> Sample Annotation (metadata):
#>   Columns: 14
#>   Names  : sample_id, neg, pos, run_day, box_id, lot, lcms_polar, lcms_neg, lcms_pos_early, lcms_pos_late, client_identifier, sex, age, bmi
#> 
#> Feature Annotation (metadata):
#>   Columns: 14
#>   Names  : feature_id, metabolite_id, comp_id, platform, super_pathway, sub_pathway, kegg, chem_id, pathway_sortorder, ri, mass, cas, pubchem, group_hmdb
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

## Run batch normalisation

Early Metabolon releases only provided the original scale data and a
batch normalised and minimally imputed data set. We would not advocate
using minimally imputed data, particularly not for quality control /
pre-processing. As such, you may have to (1) batch normalize the
original scale data, as we are doing here. Or (2) convert the minimally
imputed cells in the ScaledImp tab of the release back into NAs.

Again, here we are illustrating how to run a mass-spectroscopy style
median batch normalization.

``` r

mydata <- mydata |>
  batch_normalise(run_mode_col          = "platform",  ## Column names in feature sheet that defines the run mode for each feature
                  run_mode_colmap       = 
                    c("LC/MS Neg"       = 'lcms_neg' , ## string from feature sheet column = 
                                                       ## column name in sample sheet that holds the batch information
                      "LC/MS Polar"     = 'lcms_polar',
                      "LC/MS Pos Early" = 'lcms_pos_early',
                      "LC/MS Pos Late"  = 'lcms_pos_late' ) 
                  )
```

## Accessing data

### Raw input data

``` r

mydata@data[1:5, 1:5, "input"]
#>          48719    43532     46639       606     62279
#> sam_1 39145556       NA 1667316.8 129141879  317923.4
#> sam_2 46799204       NA  490729.6 136439235 2111209.0
#> sam_3       NA 320490.0        NA  96088900 3407171.5
#> sam_4  4094175 494174.2  214858.1  41081520        NA
#> sam_5  3908151 127978.8  238216.9  42057264  623911.2
```

### Batch normalised data

``` r

mydata@data[1:5, 1:5, "batch_normalised"]
#>           48719    43532      46639       606     62279
#> sam_1 1.5376548       NA 0.26359902 0.9248251 0.5003213
#> sam_2 1.8382935       NA 0.07758324 0.9770837 1.0718382
#> sam_3        NA 1.147501         NA 0.6881224 5.3619212
#> sam_4 0.9071567 7.172319 0.41313486 1.0725661        NA
#> sam_5 0.8659388 1.857452 0.45804982 1.0980411 0.9818592
```

## Now any quality control step should be run on the “batch_normalised” layer.

``` r

## identify the xenos
xenos <- mydata@features[!is.na(mydata@features$super_pathway) & mydata@features$super_pathway == "Xenobiotics", "feature_id"]

## run the QC on the batch_normalised layer
mydata <- mydata |> quality_control(source_layer               = "batch_normalised", 
                                    sample_missingness         = 0.2, 
                                    feature_missingness        = 0.2, 
                                    total_sum_abundance_sd     = 5, 
                                    outlier_udist              = 5, 
                                    outlier_treatment          = "leave_be", 
                                    winsorize_quantile         = 1.0, 
                                    tree_cut_height            = 0.5, 
                                    pc_outlier_sd              = 5,
                                    feature_selection          = "least_missingness", 
                                    features_exclude_but_keep  = xenos, 
                                    cores                      = 1
                                    )
#> 
#> ── Starting Omics QC Process ───────────────────────────────────────────────────
#> ℹ Validating input parameters
#> ✔ Validating input parameters [8ms]
#> 
#> ℹ Excluding 0 features from sample summary analysis but keeping in output data
#> ✔ Excluding 39 features from sample summary analysis but keeping in output data…
#> 
#> ℹ Sample & Feature Summary Statistics for raw data
#> ℹ Number of informative PCs (Scree acceleration factor): 2
#> ℹ Sample & Feature Summary Statistics for raw data✔ Sample & Feature Summary Statistics for raw data [1.7s]
#> 
#> ℹ Copying batch_normalised data to new 'qc' data layer
#> ✔ Copying batch_normalised data to new 'qc' data layer [20ms]
#> 
#> ℹ Assessing for extreme sample missingness >=80% - excluding 0 sample(s)
#> ✔ Assessing for extreme sample missingness >=80% - excluding 0 sample(s) [17ms]
#> 
#> ℹ Assessing for extreme feature missingness >=80% - excluding 0 feature(s)
#> ✔ Assessing for extreme feature missingness >=80% - excluding 4 feature(s) [18m…
#> 
#> ℹ Assessing for sample missingness at specified level of >=20% - excluding 0 sa…
#> ✔ Assessing for sample missingness at specified level of >=20% - excluding 1 sa…
#> 
#> ℹ Assessing for feature missingness at specified level of >=20% - excluding 0 f…
#> ✔ Assessing for feature missingness at specified level of >=20% - excluding 37 …
#> 
#> ℹ Calculating total sum abundance outliers at +/- 5 Sdev - excluding 0 sample(s)
#> ✔ Calculating total sum abundance outliers at +/- 5 Sdev - excluding 0 sample(s…
#> 
#> ℹ Running sample data PCA outlier analysis at +/- 5 Sdev
#> ✔ Running sample data PCA outlier analysis at +/- 5 Sdev [17ms]
#> 
#> ℹ Sample PCA outlier analysis - re-identify feature independence and PC outlier…
#> ℹ Number of informative PCs (Scree acceleration factor): 2
#> ℹ Sample PCA outlier analysis - re-identify feature independence and PC outlier…! The stated max PCs [max_num_pcs=10] to use in PCA outlier assessment is greater than the number of available informative PCs [2]
#> ℹ Sample PCA outlier analysis - re-identify feature independence and PC outlier…✔ Sample PCA outlier analysis - re-identify feature independence and PC outlier…
#> 
#> ℹ Creating final QC dataset...
#> ℹ Number of informative PCs (Scree acceleration factor): 2
#> ℹ Creating final QC dataset...
#> ℹ Creating final QC dataset...── Step timings ──
#> ℹ Creating final QC dataset...
#> ℹ Creating final QC dataset...
#>                         step seconds   pct
#>                   validation    0.00   0.0
#>                summarise_raw    1.70  36.6
#>                   copy_layer    0.00   0.0
#>   extreme_sample_missingness    0.00   0.0
#>  extreme_feature_missingness    0.00   0.0
#>           sample_missingness    0.00   0.0
#>          total_sum_abundance    0.01   0.2
#>                summarise_pca    1.42  30.6
#>              summarise_final    1.29  27.8
#>                        total    4.64 100.0
#> ✔ Creating final QC dataset... [1.3s]
#> 
#> ℹ 'Omics QC Process Completed
#> ✔ 'Omics QC Process Completed [21ms]
```

### QC data

``` r

mydata@data[1:5, 1:5, "qc"]
#>       48719    43532      46639       606     62279
#> sam_1    NA       NA 0.26359902 0.9248251 0.5003213
#> sam_2    NA       NA 0.07758324 0.9770837 1.0718382
#> sam_3    NA 1.147501         NA 0.6881224 5.3619212
#> sam_4    NA 7.172319 0.41313486 1.0725661        NA
#> sam_5    NA 1.857452 0.45804982 1.0980411 0.9818592
```

## Quick view of omiprep object

### data array layers

``` r

dimnames(mydata@data)[[3]]
#> [1] "input"            "batch_normalised" "qc"
```

### omiprep summary

``` r

summary(mydata)
#> Omiprep Object Summary
#> --------------------------
#> Samples      : 125
#> Features     : 253
#> Data Layers  : 3
#> Layer Names  : input, batch_normalised, qc
#> 
#> Sample Summary Layers : batch_normalised, qc
#> Feature Summary Layers: batch_normalised, qc
#> 
#> Sample Annotation (metadata):
#>   Columns: 16
#>   Names  : sample_id, neg, pos, run_day, box_id, lot, lcms_polar, lcms_neg, lcms_pos_early, lcms_pos_late, client_identifier, sex, age, bmi, reason_excluded, excluded
#> 
#> Feature Annotation (metadata):
#>   Columns: 16
#>   Names  : feature_id, metabolite_id, comp_id, platform, super_pathway, sub_pathway, kegg, chem_id, pathway_sortorder, ri, mass, cas, pubchem, group_hmdb, reason_excluded, excluded
#> 
#> Exclusion Codes Summary:
#> 
#>   Sample Exclusions:
#> Exclusion | Count
#> -----------------
#> user_excluded                     | 0
#> extreme_sample_missingness        | 0
#> user_defined_sample_missingness   | 1
#> user_defined_sample_totalpeakarea | 0
#> user_defined_sample_pca_outlier   | 2
#> 
#>   Feature Exclusions:
#> Exclusion | Count
#> -----------------
#> user_excluded                    |  0
#> extreme_feature_missingness      |  4
#> user_defined_feature_missingness | 37
#> user_defined_feature_skewness    |  0
```
