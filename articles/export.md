# Export Data

`Omiprep` can export data to various formats.

## Setup

### load the omiprep library

``` r
library(omiprep)
```

### Read in the data and make a Omiprep object

Create a `Omiprep` object as described in the [Getting
Started](https://mrcieu.github.io/omiprep/articles/index.md) vignette.

``` r
# read in the metabolon data as a list object
datain     <- read_metabolon(system.file("extdata", "metabolon_v1.1_example.xlsx", package = "omiprep"), 
                            sheet="OrigScale", 
                            return_Omiprep = FALSE)

# build the Omiprep class object
mydata      <- Omiprep(data = datain$data, samples = datain$samples, features = datain$features)
```

## Run the quality control

``` r
## Adding suppressWarnings() to avoid deparse() error when rendering vignette with S7 method warnings
mydata         <- suppressWarnings( quality_control(mydata, cores = 1) )
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
#> ℹ Sample & Feature Summary Statistics for raw data
#> AF =  2
#> ✔ Sample & Feature Summary Statistics for raw data [497ms]
#> 
#> ℹ Copying input data to new 'qc' data layer
#> ✔ Copying input data to new 'qc' data layer [25ms]
#> 
#> ℹ Assessing for extreme sample missingness >=80% - excluding 0 sample(s)
#> ✔ Assessing for extreme sample missingness >=80% - excluding 0 sample(s) [23ms]
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
#> ℹ Calculating total peak abundance outliers at +/- 5 Sdev - excluding 0 sample(…
#> ✔ Calculating total peak abundance outliers at +/- 5 Sdev - excluding 0 sample(…
#> 
#> ℹ Running sample data PCA outlier analysis at +/- 5 Sdev
#> ✔ Running sample data PCA outlier analysis at +/- 5 Sdev [24ms]
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
#>                   validation    0.02   1.2
#>                summarise_raw    0.48  27.8
#>                   copy_layer    0.00   0.0
#>   extreme_sample_missingness    0.00   0.0
#>  extreme_feature_missingness    0.00   0.0
#>           sample_missingness    0.00   0.0
#>              total_peak_area    0.00   0.0
#>                summarise_pca    0.54  31.3
#>              summarise_final    0.46  26.7
#>                        total    1.73 100.2
#> ✔ Creating final QC dataset... [505ms]
#> 
#> ℹ 'Omics QC Process Completed
#> ✔ 'Omics QC Process Completed [26ms]
```

## Export Omiprep

``` r
# where to put the files
output_dir <- file.path(getwd(), "output")

# run export
export(mydata, directory = output_dir, format = "omiprep")
#> Exporting in omiprep format to: 
#>      /home/runner/work/omiprep/omiprep/vignettes/output

# view output directory files
files <- list.files(output_dir, full.names = TRUE, recursive = TRUE)
unname(sapply(files, function(path) {
  parts <- strsplit(path, .Platform$file.sep)[[1]]
  paste(tail(parts, 4), collapse = .Platform$file.sep)
}))
#>  [1] "output/omiprep_export_2026_04_20/input/config.yml"         
#>  [2] "output/omiprep_export_2026_04_20/input/data.tsv"           
#>  [3] "output/omiprep_export_2026_04_20/input/feature_summary.tsv"
#>  [4] "output/omiprep_export_2026_04_20/input/features.tsv"       
#>  [5] "output/omiprep_export_2026_04_20/input/sample_summary.tsv" 
#>  [6] "output/omiprep_export_2026_04_20/input/samples.tsv"        
#>  [7] "output/omiprep_export_2026_04_20/qc/config.yml"            
#>  [8] "output/omiprep_export_2026_04_20/qc/data.tsv"              
#>  [9] "output/omiprep_export_2026_04_20/qc/feature_summary.tsv"   
#> [10] "output/omiprep_export_2026_04_20/qc/feature_tree.RDS"      
#> [11] "output/omiprep_export_2026_04_20/qc/features.tsv"          
#> [12] "output/omiprep_export_2026_04_20/qc/sample_summary.tsv"    
#> [13] "output/omiprep_export_2026_04_20/qc/samples.tsv"           
#> [14] "output/omiprep_export_2026_04_20/qc/var_exp.tsv"           
#> [15] "output/omiprep_export_2026_04_28/input/config.yml"         
#> [16] "output/omiprep_export_2026_04_28/input/data.tsv"           
#> [17] "output/omiprep_export_2026_04_28/input/feature_summary.tsv"
#> [18] "output/omiprep_export_2026_04_28/input/features.tsv"       
#> [19] "output/omiprep_export_2026_04_28/input/sample_summary.tsv" 
#> [20] "output/omiprep_export_2026_04_28/input/samples.tsv"        
#> [21] "output/omiprep_export_2026_04_28/qc/config.yml"            
#> [22] "output/omiprep_export_2026_04_28/qc/data.tsv"              
#> [23] "output/omiprep_export_2026_04_28/qc/feature_summary.tsv"   
#> [24] "output/omiprep_export_2026_04_28/qc/feature_tree.RDS"      
#> [25] "output/omiprep_export_2026_04_28/qc/features.tsv"          
#> [26] "output/omiprep_export_2026_04_28/qc/sample_summary.tsv"    
#> [27] "output/omiprep_export_2026_04_28/qc/samples.tsv"           
#> [28] "output/omiprep_export_2026_04_28/qc/var_exp.tsv"
```
