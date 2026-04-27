# Import Metabolon Metabolomic Data

## Import Metabolon data

Read in the Metabolon data using the `read_metabolon` function. Here we
will read in the example data provided with the package, as a list
object.

``` r
library(omiprep)

# example file
filepath <- system.file("extdata", "metabolon_v1.2_example.xlsx", package = "omiprep")

# import data as a list object rather than directly as a Omiprep object
dat <- read_metabolon(filepath, 
                      sheet = 'OrigScale', 
                      return_Omiprep = FALSE)
```

## Quick look at data structure of the imported data

``` r
str(dat)
#> List of 3
#>  $ data    : num [1:100, 1:104] 98551 43695 44899 37811 36825 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:100] "ind1" "ind2" "ind3" "ind4" ...
#>   .. ..$ : chr [1:104] "123" "124" "125" "126" ...
#>  $ samples :'data.frame':    100 obs. of  5 variables:
#>   ..$ sample_id        : chr [1:100] "ind1" "ind2" "ind3" "ind4" ...
#>   ..$ parent_sample_id : chr [1:100] "ps_id" "ps_id" "ps_id" "ps_id" ...
#>   ..$ client_identifier: chr [1:100] "FR01234" "FR01235" "FR01236" "FR01237" ...
#>   ..$ pair             : chr [1:100] "99999" "99999" "99999" "99999" ...
#>   ..$ volume_extracted : chr [1:100] "100" "100" "100" "100" ...
#>  $ features:'data.frame':    104 obs. of  14 variables:
#>   ..$ feature_id       : chr [1:104] "123" "124" "125" "126" ...
#>   ..$ pathway_sortorder: chr [1:104] "1" "2" "3" "4" ...
#>   ..$ biochemical      : chr [1:104] "(N(1) + N(8))-acetylspermidine" "1,2,3-benzenetriol sulfate (2)" "1,2-dilinoleoyl-GPC (18:2/18:2)" "1,2-dilinoleoyl-GPE (18:2/18:2)*" ...
#>   ..$ super_pathway    : chr [1:104] "Amino Acid" "Xenobiotics" "Lipid" "Lipid" ...
#>   ..$ sub_pathway      : chr [1:104] "Polyamine Metabolism" "Chemical" "Phosphatidylcholine (PC)" "Phosphatidylethanolamine (PE)" ...
#>   ..$ comp_id          : chr [1:104] "123" "124" "125" "126" ...
#>   ..$ platform         : chr [1:104] "LC/MS Pos Early" "LC/MS Neg" "LC/MS Pos Late" "LC/MS Pos Late" ...
#>   ..$ chemical_id      : chr [1:104] "1111" "1112" "1113" "1114" ...
#>   ..$ ri               : chr [1:104] "2221" "2222" "2223" "2224" ...
#>   ..$ mass             : chr [1:104] "111.111" "111.11199999999999" "111.113" "111.114" ...
#>   ..$ cas              : chr [1:104] NA NA "111-11-1" NA ...
#>   ..$ pubchem          : chr [1:104] NA NA "11111" "11112" ...
#>   ..$ kegg             : chr [1:104] NA NA NA NA ...
#>   ..$ group_hmdb       : chr [1:104] NA NA "HMDB123" "HMDB124" ...
```

## Create Omiprep object

Once imported, we pass the data to the Omiprep() function to build the
`Omiprep` class object.

``` r
## This step could be avoided by defining return_Omiprep = TRUE in read_metabolon() function above.
mydata <- Omiprep(data     = dat$data, 
                  features = dat$features, 
                  samples  = dat$samples)
```

## Quick summary of the Omiprep object

``` r
summary(mydata)
#> Omiprep Object Summary
#> --------------------------
#> Samples      : 100
#> Features     : 104
#> Data Layers  : 1
#> Layer Names  : input
#> 
#> Sample Summary Layers : none
#> Feature Summary Layers: none
#> 
#> Sample Annotation (metadata):
#>   Columns: 5
#>   Names  : sample_id, parent_sample_id, client_identifier, pair, volume_extracted
#> 
#> Feature Annotation (metadata):
#>   Columns: 14
#>   Names  : feature_id, pathway_sortorder, biochemical, super_pathway, sub_pathway, comp_id, platform, chemical_id, ri, mass, cas, pubchem, kegg, group_hmdb
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
xenos <- mydata@features[!is.na(mydata@features$super_pathway) & 
                           mydata@features$super_pathway == "Xenobiotics", "feature_id"]

## how many xenobiotics identified
length(xenos)
#> [1] 7
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
  quality_control(source_layer        = "input", 
                  sample_missingness  = 0.2, 
                  feature_missingness = 0.2, 
                  total_peak_area_sd  = 5, 
                  outlier_udist       = 5, 
                  outlier_treatment   = "leave_be", 
                  winsorize_quantile  = 1.0, 
                  tree_cut_height     = 0.5, 
                  pc_outlier_sd       = 5,
                  feature_selection   = "least_missingness", ## We suggest using `least_missingness` when working with data, like Metabolon, with high missingness. Default is "max_var_exp".
                  features_exclude_but_keep = xenos, ## exclude xenobiotics from QC, but retain them in the final dataset
                  cores               = 1
                  )
#> 
#> ── Starting Omics QC Process ───────────────────────────────────────────────────
#> ℹ Validating input parameters
#> 
#> ℹ Validating input parameters── Starting 'Omics QC Process ──────────────────────────────────────────────────
#> ℹ Validating input parameters✔ Validating input parameters [17ms]
#> 
#> ℹ Validating input parameters
#> ✔ Validating input parameters [14ms]
#> 
#> ℹ Excluding 0 features from sample summary analysis but keeping in output data
#> ✔ Excluding 7 features from sample summary analysis but keeping in output data …
#> 
#> ℹ Sample & Feature Summary Statistics for raw data
#> AF =  2
#> ✔ Sample & Feature Summary Statistics for raw data [546ms]
#> 
#> ℹ Copying input data to new 'qc' data layer
#> ✔ Copying input data to new 'qc' data layer [27ms]
#> 
#> ℹ Assessing for extreme sample missingness >=80% - excluding 0 sample(s)
#> ✔ Assessing for extreme sample missingness >=80% - excluding 1 sample(s) [22ms]
#> 
#> ℹ Assessing for extreme feature missingness >=80% - excluding 0 feature(s)
#> ✔ Assessing for extreme feature missingness >=80% - excluding 0 feature(s) [20m…
#> 
#> ℹ Assessing for sample missingness at specified level of >=20% - excluding 0 sa…
#> ✔ Assessing for sample missingness at specified level of >=20% - excluding 0 sa…
#> 
#> ℹ Assessing for feature missingness at specified level of >=20% - excluding 0 f…
#> ✔ Assessing for feature missingness at specified level of >=20% - excluding 1 f…
#> 
#> ℹ Calculating total peak abundance outliers at +/- 5 Sdev - excluding 0 sample(…
#> ✔ Calculating total peak abundance outliers at +/- 5 Sdev - excluding 0 sample(…
#> 
#> ℹ Running sample data PCA outlier analysis at +/- 5 Sdev
#> ✔ Running sample data PCA outlier analysis at +/- 5 Sdev [17ms]
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
#>                summarise_raw    0.53  27.1
#>                   copy_layer    0.00   0.0
#>   extreme_sample_missingness    0.00   0.0
#>  extreme_feature_missingness    0.00   0.0
#>           sample_missingness    0.00   0.0
#>              total_peak_area    0.00   0.0
#>                summarise_pca    0.60  30.7
#>              summarise_final    0.53  27.1
#>                        total    1.96 100.2
#> ✔ Creating final QC dataset... [592ms]
#> 
#> ℹ 'Omics QC Process Completed
#> ✔ 'Omics QC Process Completed [22ms]
```

## Quick summary of the Omiprep object following QC

``` r
summary(mydata)
#> Omiprep Object Summary
#> --------------------------
#> Samples      : 100
#> Features     : 104
#> Data Layers  : 2
#> Layer Names  : input, qc
#> 
#> Sample Summary Layers : input, qc
#> Feature Summary Layers: input, qc
#> 
#> Sample Annotation (metadata):
#>   Columns: 7
#>   Names  : sample_id, parent_sample_id, client_identifier, pair, volume_extracted, reason_excluded, excluded
#> 
#> Feature Annotation (metadata):
#>   Columns: 16
#>   Names  : feature_id, pathway_sortorder, biochemical, super_pathway, sub_pathway, comp_id, platform, chemical_id, ri, mass, cas, pubchem, kegg, group_hmdb, reason_excluded, excluded
#> 
#> Exclusion Codes Summary:
#> 
#>   Sample Exclusions:
#> Exclusion | Count
#> -----------------
#> user_excluded                     | 0
#> extreme_sample_missingness        | 1
#> user_defined_sample_missingness   | 0
#> user_defined_sample_totalpeakarea | 0
#> user_defined_sample_pca_outlier   | 0
#> 
#>   Feature Exclusions:
#> Exclusion | Count
#> -----------------
#> user_excluded                    | 0
#> extreme_feature_missingness      | 0
#> user_defined_feature_missingness | 1
#> user_defined_feature_skewness    | 0
```
