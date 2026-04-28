# Generate QC Report

## Import example metabolomics data

``` r
library(omiprep)

# import data directly as a Omiprep object
mydata <- read_metabolon(system.file("extdata", "metabolon_v1.1_example.xlsx", package = "omiprep"), 
                         sheet="OrigScale",
                         return_Omiprep = TRUE)
```

## Identify the Xenobiotics to exclude from the QC steps

``` r
xenos <- mydata@features[mydata@features$pathway == "xenobiotic", "feature_id"]

## how many xenobiotics identified
length(xenos)
#> [1] 3
```

## QC the example Metabolon data

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
                  feature_selection   = "least_missingness", 
                  features_exclude_but_keep = xenos,
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
#> ✔ Validating input parameters [13ms]
#> 
#> ℹ Excluding 0 features from sample summary analysis but keeping in output data
#> ✔ Excluding 3 features from sample summary analysis but keeping in output data …
#> 
#> ℹ Sample & Feature Summary Statistics for raw data
#> AF =  2
#> ✔ Sample & Feature Summary Statistics for raw data [523ms]
#> 
#> ℹ Copying input data to new 'qc' data layer
#> ✔ Copying input data to new 'qc' data layer [33ms]
#> 
#> ℹ Assessing for extreme sample missingness >=80% - excluding 0 sample(s)
#> ✔ Assessing for extreme sample missingness >=80% - excluding 0 sample(s) [20ms]
#> 
#> ℹ Assessing for extreme feature missingness >=80% - excluding 0 feature(s)
#> ✔ Assessing for extreme feature missingness >=80% - excluding 0 feature(s) [16m…
#> 
#> ℹ Assessing for sample missingness at specified level of >=20% - excluding 0 sa…
#> ✔ Assessing for sample missingness at specified level of >=20% - excluding 2 sa…
#> 
#> ℹ Assessing for feature missingness at specified level of >=20% - excluding 0 f…
#> ✔ Assessing for feature missingness at specified level of >=20% - excluding 0 f…
#> 
#> ℹ Calculating total peak abundance outliers at +/- 5 Sdev - excluding 0 sample(…
#> ✔ Calculating total peak abundance outliers at +/- 5 Sdev - excluding 0 sample(…
#> 
#> ℹ Running sample data PCA outlier analysis at +/- 5 Sdev
#> ✔ Running sample data PCA outlier analysis at +/- 5 Sdev [16ms]
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
#>                         step seconds  pct
#>                   validation    0.02  1.1
#>                summarise_raw    0.51 27.5
#>                   copy_layer    0.00  0.0
#>   extreme_sample_missingness    0.00  0.0
#>  extreme_feature_missingness    0.00  0.0
#>           sample_missingness    0.00  0.0
#>              total_peak_area    0.00  0.0
#>                summarise_pca    0.58 31.3
#>              summarise_final    0.49 26.4
#>                        total    1.85 99.8
#> ✔ Creating final QC dataset... [542ms]
#> 
#> ℹ 'Omics QC Process Completed
#> ✔ 'Omics QC Process Completed [26ms]
```

## Generate the Omiprep report

``` r
# render report
generate_report(mydata,
                project         = "myproject",
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
#> /opt/hostedtoolcache/pandoc/3.1.11/x64/pandoc +RTS -K512m -RTS /home/runner/work/omiprep/omiprep/vignettes/skeleton.knit.md --to html4 --from markdown+autolink_bare_uris+tex_math_single_backslash --output /home/runner/work/omiprep/omiprep/vignettes/myproject_omiprep_qc_report.html --lua-filter /home/runner/work/_temp/Library/rmarkdown/rmarkdown/lua/pagebreak.lua --lua-filter /home/runner/work/_temp/Library/rmarkdown/rmarkdown/lua/latex-div.lua --embed-resources --standalone --variable bs3=TRUE --section-divs --table-of-contents --toc-depth 2 --template /home/runner/work/_temp/Library/rmarkdown/rmd/h/default.html --no-highlight --variable highlightjs=1 --number-sections --variable theme=bootstrap --css styles.css --mathjax --variable 'mathjax-url=https://mathjax.rstudio.com/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' --include-in-header /tmp/RtmpvdcWZI/rmarkdown-str1dad47e1362e.html
#> 
#> Output created: myproject_omiprep_qc_report.html
```

    #> [1] TRUE
