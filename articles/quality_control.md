# quality_control

``` r

library(omiprep)

# example file
filepath <- system.file("extdata", "metabolon_v2_example.xlsx", package = "omiprep")

# import data as a list object rather than directly as a Omiprep object
mydata <- read_metabolon(filepath,  
                         sheet = 'Batch-normalized Data',        ## read in the batch normalized but NOT imputed tab
                         feature_sheet = "Chemical Annotation",  ## tab name of the feature data
                         feature_id_col = "COMP_ID",             ## in this instance features in batch normalized.. 
                                                                 ##  ..tab are annotated to column name COMP_ID. Verify in your data.
                         sample_sheet = "Sample Meta Data",      ## tab name of the sample data
                         sample_id_col = "PARENT_SAMPLE_NAME",   ## column name of sample IDs 
                         return_Omiprep = TRUE                   ## return omiprep object
                         )
```

## Identify the Xenobiotics to exclude from the QC steps

``` r

xenos <- mydata@features[!is.na(mydata@features$SUPER_PATHWAY) & 
                           mydata@features$SUPER_PATHWAY == "Xenobiotics", "feature_id"]

## how many xenobiotics identified
length(xenos)
#> [1] 39
```

### Run the quality control pipeline

``` r

# run QC
mydata <- mydata |>
  quality_control(source_layer              = "input", 
                  sample_missingness        = 0.2, 
                  feature_missingness       = 0.2, 
                  total_sum_abundance_sd    = 5, 
                  outlier_udist             = 5, 
                  outlier_treatment         = "leave_be", 
                  winsorize_quantile        = 1.0, 
                  tree_cut_height           = 0.5, 
                  pc_outlier_sd             = 5,
                  feature_selection         = "least_missingness", ## We suggest using `least_missingness` 
                                                                   ## when working with data, like Metabolon, 
                                                                   ## with high missingness. 
                                                                   ## Default is "max_var_exp".
                  features_exclude_but_keep = xenos, ## Exclude xenobiotics from QC, but retain them in the final dataset.
                                                     ## We suggest this as xenos can have extreme missingness and are 
                                                     ## commonly qc'd from a data set. However, they may be appropriate 
                                                     ## to model as a binary present|absent trait. A choice for the researcher.
                  cores                     = 1
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
#> ℹ Sample & Feature Summary Statistics for raw data✔ Sample & Feature Summary Statistics for raw data [1.6s]
#> 
#> ℹ Copying input data to new 'qc' data layer
#> ✔ Copying input data to new 'qc' data layer [20ms]
#> 
#> ℹ Assessing for extreme sample missingness >=80% - excluding 0 sample(s)
#> ✔ Assessing for extreme sample missingness >=80% - excluding 0 sample(s) [17ms]
#> 
#> ℹ Assessing for extreme feature missingness >=80% - excluding 0 feature(s)
#> ✔ Assessing for extreme feature missingness >=80% - excluding 4 feature(s) [19m…
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
#> ✔ Running sample data PCA outlier analysis at +/- 5 Sdev [22ms]
#> 
#> ℹ Sample PCA outlier analysis - re-identify feature independence and PC outlier…
#> ℹ Number of informative PCs (Scree acceleration factor): 2
#> ℹ Sample PCA outlier analysis - re-identify feature independence and PC outlier…ℹ Sample PCA outlier analysis - re-identify feature independence and PC outlier…
#> ! The stated max PCs [max_num_pcs=10] to use in PCA outlier assessment is greater than the number of available informative PCs [2]
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
#>                summarise_raw    1.56  31.3
#>                   copy_layer    0.00   0.0
#>   extreme_sample_missingness    0.00   0.0
#>  extreme_feature_missingness    0.00   0.0
#>           sample_missingness    0.00   0.0
#>          total_sum_abundance    0.01   0.2
#>                summarise_pca    1.66  33.3
#>              summarise_final    1.52  30.5
#>                        total    4.98 100.0
#> ✔ Creating final QC dataset... [1.6s]
#> 
#> ℹ 'Omics QC Process Completed
#> ✔ 'Omics QC Process Completed [14ms]
```

### View a summary of the Omiprep object

``` r

# view summary
summary(mydata)
#> Omiprep Object Summary
#> --------------------------
#> Samples      : 125
#> Features     : 253
#> Data Layers  : 2
#> Layer Names  : input, qc
#> 
#> Sample Summary Layers : input, qc
#> Feature Summary Layers: input, qc
#> 
#> Sample Annotation (metadata):
#>   Columns: 12
#>   Names  : sample_id, PARENT_SAMPLE_NAME, CLIENT_IDENTIFIER, LC/MS Polar, LC/MS Neg, LC/MS Pos Early, LC/MS Pos Late, sex, age, bmi, reason_excluded, excluded
#> 
#> Feature Annotation (metadata):
#>   Columns: 16
#>   Names  : feature_id, COMP_ID, CHEM_ID, CHEMICAL_NAME, SUPER_PATHWAY, SUB_PATHWAY, PATHWAY_SORTORDER, PLATFORM, RI, MASS, CAS, PUBCHEM, KEGG, HMDB, reason_excluded, excluded
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
