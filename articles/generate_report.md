# Generate QC Report

## Import example metabolomics data

``` r

library(omiprep)

## example data file path
filepath <- system.file("extdata", "metabolon_v2_example.xlsx", package = "omiprep")

# import data directly as a Omiprep object
mydata <- read_metabolon(filepath,  
                         sheet = 'Batch-normalized Data',        ## read in the batch normalized but NOT imputed tab
                         feature_sheet = "Chemical Annotation",  ## tab name of the feature data
                         feature_id_col = "COMP_ID",             ## column name of feature IDs
                         sample_sheet = "Sample Meta Data",      ## tab name of the sample data
                         sample_id_col = "PARENT_SAMPLE_NAME",   ## column name of sample IDs 
                         return_Omiprep = TRUE
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

## QC the example Metabolon data

``` r

mydata <- mydata |>
  quality_control(source_layer        = "input", 
                  sample_missingness  = 0.2, 
                  feature_missingness = 0.2, 
                  total_sum_abundance_sd  = 5, 
                  outlier_udist       = 5, 
                  outlier_treatment   = "leave_be", 
                  winsorize_quantile  = 1.0, 
                  tree_cut_height     = 0.5, 
                  pc_outlier_sd       = 5,
                  feature_selection   = "least_missingness", 
                  features_exclude_but_keep = xenos,
                  cores               = 1
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
#>                summarise_raw    1.45  28.8
#>                   copy_layer    0.00   0.0
#>   extreme_sample_missingness    0.00   0.0
#>  extreme_feature_missingness    0.00   0.0
#>           sample_missingness    0.00   0.0
#>          total_sum_abundance    0.01   0.2
#>                summarise_pca    2.03  40.4
#>              summarise_final    1.31  26.1
#>                        total    5.03 100.1
#> ✔ Creating final QC dataset... [1.3s]
#> 
#> ℹ 'Omics QC Process Completed
#> ✔ 'Omics QC Process Completed [14ms]
```

## Generate the Omiprep report

``` r

# render report
project     <- "metabolon_example"
report_html <- paste0(project, "_omiprep_qc_report.html")   # name generate_report() writes

generate_report(mydata,
                project         = project,
                output_dir      = getwd(),
                output_filename = NULL,
                format          = "html",
                template        = "qc_report")
#> processing file: skeleton.Rmd
#> Warning in call_block(x): The chunk 'unnamed-chunk-1' has the 'child' option,
#> and this code chunk must be empty. Its code will be ignored.
#> Warning in call_block(x): The chunk 'unnamed-chunk-2' has the 'child' option,
#> and this code chunk must be empty. Its code will be ignored.
#> Warning in call_block(x): The chunk 'unnamed-chunk-3' has the 'child' option,
#> and this code chunk must be empty. Its code will be ignored.
#> output file: /home/runner/work/omiprep/omiprep/vignettes/skeleton.knit.md
#> /opt/hostedtoolcache/pandoc/3.8.3/x64/pandoc +RTS -K512m -RTS /home/runner/work/omiprep/omiprep/vignettes/skeleton.knit.md --to html4 --from markdown+autolink_bare_uris+tex_math_single_backslash --output /home/runner/work/omiprep/omiprep/vignettes/metabolon_example_omiprep_qc_report.html --lua-filter /home/runner/work/_temp/Library/rmarkdown/rmarkdown/lua/pagebreak.lua --lua-filter /home/runner/work/_temp/Library/rmarkdown/rmarkdown/lua/latex-div.lua --lua-filter /home/runner/work/_temp/Library/rmarkdown/rmarkdown/lua/table-classes.lua --embed-resources --standalone --variable bs3=TRUE --section-divs --table-of-contents --toc-depth 2 --template /home/runner/work/_temp/Library/rmarkdown/rmd/h/default.html --syntax-highlighting none --variable highlightjs=1 --number-sections --variable theme=bootstrap --css styles.css --mathjax --variable 'mathjax-url=https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --include-in-header /tmp/RtmpR45bty/rmarkdown-str1d7b5bb837ec.html
#> 
#> Output created: metabolon_example_omiprep_qc_report.html
```

    #> [1] TRUE
