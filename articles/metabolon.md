# Import Metabolon Metabolomic Data

## Import Metabolon data

Read in the Metabolon data using the `read_metabolon` function. Here we
will read in the example data provided with the package, as a list
object.

``` r

library(omiprep)

# example file
filepath <- system.file("extdata", "metabolon_v2_example.xlsx", package = "omiprep")

# import data as a list object rather than directly as a Omiprep object
datain <- read_metabolon(filepath,  
                         sheet = 'Batch-normalized Data',        ## read in the batch normalized but NOT imputed tab
                         feature_sheet = "Chemical Annotation",  ## tab name of the feature data
                         feature_id_col = "COMP_ID",             ## in this instance features in batch normalized.. 
                                                                 ##  ..tab are annotated to column name COMP_ID. Verify in your data.
                         sample_sheet = "Sample Meta Data",      ## tab name of the sample data
                         sample_id_col = "PARENT_SAMPLE_NAME",   ## column name of sample IDs 
                         return_Omiprep = FALSE                  ## return data as a list of data, feature, and samples, ..
                                                                 ##  .. NOT an omiprep object  
                         )
```

## Quick look at data structure of the imported data

``` r

str(datain)
#> List of 3
#>  $ data    : num [1:125, 1:253] 1.538 1.838 NA 0.907 0.866 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:125] "sam_1" "sam_2" "sam_3" "sam_4" ...
#>   .. ..$ : chr [1:253] "48719" "43532" "46639" "606" ...
#>  $ samples :'data.frame':    125 obs. of  10 variables:
#>   ..$ PARENT_SAMPLE_NAME: chr [1:125] "sam_1" "sam_2" "sam_3" "sam_4" ...
#>   ..$ CLIENT_IDENTIFIER : chr [1:125] "1_serum_1848" "2_serum_1844" "3_serum_1089" "4_serum_1467" ...
#>   ..$ LC/MS Polar       : chr [1:125] "batch1" "batch1" "batch2" "batch1" ...
#>   ..$ LC/MS Neg         : chr [1:125] "batch2" "batch2" "batch2" "batch1" ...
#>   ..$ LC/MS Pos Early   : chr [1:125] "batch1" "batch2" "batch1" "batch2" ...
#>   ..$ LC/MS Pos Late    : chr [1:125] "batch1" "batch2" "batch1" "batch1" ...
#>   ..$ sex               : chr [1:125] "male" "male" "male" "female" ...
#>   ..$ age               : num [1:125] 50 42 43 54 51 35 59 36 40 57 ...
#>   ..$ bmi               : num [1:125] 17.1 27.3 30.6 29.2 24.4 ...
#>   ..$ sample_id         : chr [1:125] "sam_1" "sam_2" "sam_3" "sam_4" ...
#>  $ features:'data.frame':    253 obs. of  14 variables:
#>   ..$ COMP_ID          : num [1:253] 48719 43532 46639 606 62279 ...
#>   ..$ CHEM_ID          : num [1:253] NA 1.00e+08 NA 5.35e+02 1.00e+08 ...
#>   ..$ CHEMICAL_NAME    : chr [1:253] "X - 23145" "veratric acid" "X - 12851" "uridine" ...
#>   ..$ SUPER_PATHWAY    : chr [1:253] NA "Xenobiotics" NA "Nucleotide" ...
#>   ..$ SUB_PATHWAY      : chr [1:253] NA "Drug - Neurological" NA "Pyrimidine Metabolism, Uracil containing" ...
#>   ..$ PATHWAY_SORTORDER: num [1:253] NA 4533 NA 3428 4751 ...
#>   ..$ PLATFORM         : chr [1:253] "LC/MS Neg" "LC/MS Neg" "LC/MS Neg" "LC/MS Neg" ...
#>   ..$ RI               : chr [1:253] "5642" "1890" "4844" "1457.6" ...
#>   ..$ MASS             : chr [1:253] "283.26420000000002" "181.0506" "311.66579999999999" "243.06229999999999" ...
#>   ..$ CAS              : chr [1:253] NA "93-07-2" NA "58-96-8" ...
#>   ..$ PUBCHEM          : num [1:253] NA 7121 NA 6029 NA ...
#>   ..$ KEGG             : chr [1:253] NA NA NA "C00299" ...
#>   ..$ HMDB             : chr [1:253] NA NA NA "HMDB00296" ...
#>   ..$ feature_id       : chr [1:253] "48719" "43532" "46639" "606" ...
```

## Create Omiprep object

Once imported, we pass the data to the Omiprep() function to build the
`Omiprep` class object.

``` r

## This step could be avoided by defining return_Omiprep = TRUE in read_metabolon() function above.
mydata <- Omiprep(data     = datain$data, 
                  features = datain$features, 
                  samples  = datain$samples)
```

## Quick summary of the Omiprep object

``` r

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
#>   Columns: 10
#>   Names  : PARENT_SAMPLE_NAME, CLIENT_IDENTIFIER, LC/MS Polar, LC/MS Neg, LC/MS Pos Early, LC/MS Pos Late, sex, age, bmi, sample_id
#> 
#> Feature Annotation (metadata):
#>   Columns: 14
#>   Names  : COMP_ID, CHEM_ID, CHEMICAL_NAME, SUPER_PATHWAY, SUB_PATHWAY, PATHWAY_SORTORDER, PLATFORM, RI, MASS, CAS, PUBCHEM, KEGG, HMDB, feature_id
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

## Identify the Xenobiotics to exclude from the QC steps

Use the feature data just imported to identify xenobiotic metabolites.
It may be best to excluded these features from the quality-control (QC)
process. Xenobiotics typically exhibit much higher levels of missingness
than endogenous metabolites, and including them in QC can result in
excessive exclusion of both features and samples. This step will allow
you to retain these features in the final dataset, by excluding them
from QC filtering steps.

``` r

xenos <- mydata@features[!is.na(mydata@features$SUPER_PATHWAY) & 
                           mydata@features$SUPER_PATHWAY == "Xenobiotics", "feature_id"]

## how many xenobiotics identified
length(xenos)
#> [1] 39
```

## QC Metabolon

Perform the QC steps using the `quality_control` function, specifying
the xenobiotics to exclude from the QC steps.

``` r

## Given the high missingness in metabolon data, 
##   we suggest using the `least_missingness` feature selection method
##   for the identification of principle variable that will then be 
##   used in the construction of PCs. 

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
#> ✔ Validating input parameters [9ms]
#> 
#> ℹ Excluding 0 features from sample summary analysis but keeping in output data
#> ✔ Excluding 39 features from sample summary analysis but keeping in output data…
#> 
#> ℹ Sample & Feature Summary Statistics for raw data
#> ℹ Number of informative PCs (Scree acceleration factor): 2
#> ℹ Sample & Feature Summary Statistics for raw data✔ Sample & Feature Summary Statistics for raw data [1.5s]
#> 
#> ℹ Copying input data to new 'qc' data layer
#> ✔ Copying input data to new 'qc' data layer [19ms]
#> 
#> ℹ Assessing for extreme sample missingness >=80% - excluding 0 sample(s)
#> ✔ Assessing for extreme sample missingness >=80% - excluding 0 sample(s) [18ms]
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
#> ✔ Running sample data PCA outlier analysis at +/- 5 Sdev [18ms]
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
#>                summarise_raw    1.50  32.0
#>                   copy_layer    0.00   0.0
#>   extreme_sample_missingness    0.00   0.0
#>  extreme_feature_missingness    0.00   0.0
#>           sample_missingness    0.00   0.0
#>          total_sum_abundance    0.01   0.2
#>                summarise_pca    1.69  36.0
#>              summarise_final    1.25  26.6
#>                        total    4.69 100.0
#> ✔ Creating final QC dataset... [1.3s]
#> 
#> ℹ 'Omics QC Process Completed
#> ✔ 'Omics QC Process Completed [14ms]
```

## Quick summary of the Omiprep object following QC

``` r

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
