# quality_control

``` r

library(omiprep)

# import data 
m <- read_metabolon(system.file("extdata", "metabolon_v1.1_example.xlsx", package = "omiprep"), 
                    sheet = "OrigScale",      ## The name of the sheet in the excel file to read in
                    return_Omiprep = TRUE  ## Whether to return a Omiprep object (TRUE) or a list (FALSE)
                    )
```

### Run the quality control pipeline

``` r

# run QC
m <- quality_control(m, 
                     source_layer = "input", 
                     sample_missingness  = 0.2, 
                     feature_missingness = 0.2, 
                     feature_skewness_threshold = NULL,
                     feature_skewness_direction = "left",
                     total_sum_abundance_sd  = 5, 
                     outlier_udist       = 5, 
                     outlier_treatment   = "leave_be", 
                     winsorize_quantile  = 1.0, 
                     tree_cut_height     = 0.5, 
                     pc_outlier_sd       = 5, 
                     sample_ids          = NULL, 
                     feature_ids         = NULL, 
                     cores               = 1)
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
#> AF =  2
#> ✔ Sample & Feature Summary Statistics for raw data [561ms]
#> 
#> ℹ Copying input data to new 'qc' data layer
#> ✔ Copying input data to new 'qc' data layer [24ms]
#> 
#> ℹ Assessing for extreme sample missingness >=80% - excluding 0 sample(s)
#> ✔ Assessing for extreme sample missingness >=80% - excluding 0 sample(s) [25ms]
#> 
#> ℹ Assessing for extreme feature missingness >=80% - excluding 0 feature(s)
#> ✔ Assessing for extreme feature missingness >=80% - excluding 0 feature(s) [17m…
#> 
#> ℹ Assessing for sample missingness at specified level of >=20% - excluding 0 sa…
#> ✔ Assessing for sample missingness at specified level of >=20% - excluding 2 sa…
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
#> AF =  2
#> ! The stated max PCs [max_num_pcs=10] to use in PCA outlier assessment is greater than the number of available informative PCs [2]
#> ℹ Sample PCA outlier analysis - re-identify feature independence and PC outlier…✔ Sample PCA outlier analysis - re-identify feature independence and PC outlier…
#> 
#> ℹ Creating final QC dataset...
#> AF =  2
#> 
#> ℹ Creating final QC dataset...── Step timings ──
#> ℹ Creating final QC dataset...
#> ℹ Creating final QC dataset...
#>                         step seconds   pct
#>                   validation    0.02   1.0
#>                summarise_raw    0.54  28.1
#>                   copy_layer    0.00   0.0
#>   extreme_sample_missingness    0.00   0.0
#>  extreme_feature_missingness    0.00   0.0
#>           sample_missingness    0.00   0.0
#>          total_sum_abundance    0.00   0.0
#>                summarise_pca    0.62  32.2
#>              summarise_final    0.51  26.5
#>                        total    1.93 100.3
#> ✔ Creating final QC dataset... [557ms]
#> 
#> ℹ 'Omics QC Process Completed
#> ✔ 'Omics QC Process Completed [30ms]
```

### View a summary of the Omiprep object

``` r

# view summary
summary(m)
#> Omiprep Object Summary
#> --------------------------
#> Samples      : 100
#> Features     : 100
#> Data Layers  : 2
#> Layer Names  : input, qc
#> 
#> Sample Summary Layers : input, qc
#> Feature Summary Layers: input, qc
#> 
#> Sample Annotation (metadata):
#>   Columns: 8
#>   Names  : sample_id, neg, pos, run_day, box_id, lot, reason_excluded, excluded
#> 
#> Feature Annotation (metadata):
#>   Columns: 9
#>   Names  : feature_id, metabolite_id, comp_id, platform, pathway, kegg, group_hmdb, reason_excluded, excluded
#> 
#> Exclusion Codes Summary:
#> 
#>   Sample Exclusions:
#> Exclusion | Count
#> -----------------
#> user_excluded                     | 0
#> extreme_sample_missingness        | 0
#> user_defined_sample_missingness   | 2
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
