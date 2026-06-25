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
#> ✔ Validating input parameters [9ms]
#> 
#> ℹ Sample & Feature Summary Statistics for raw data
#> ℹ Number of informative PCs (Scree acceleration factor): 2
#> ℹ Sample & Feature Summary Statistics for raw data✔ Sample & Feature Summary Statistics for raw data [1.9s]
#> 
#> ℹ Copying input data to new 'qc' data layer
#> ✔ Copying input data to new 'qc' data layer [50ms]
#> 
#> ℹ Assessing for extreme sample missingness >=80% - excluding 0 sample(s)
#> ✔ Assessing for extreme sample missingness >=80% - excluding 0 sample(s) [24ms]
#> 
#> ℹ Assessing for extreme feature missingness >=80% - excluding 0 feature(s)
#> ✔ Assessing for extreme feature missingness >=80% - excluding 5 feature(s) [19m…
#> 
#> ℹ Assessing for sample missingness at specified level of >=20% - excluding 0 sa…
#> ✔ Assessing for sample missingness at specified level of >=20% - excluding 1 sa…
#> 
#> ℹ Assessing for feature missingness at specified level of >=20% - excluding 0 f…
#> ✔ Assessing for feature missingness at specified level of >=20% - excluding 46 …
#> 
#> ℹ Calculating total sum abundance outliers at +/- 5 Sdev - excluding 0 sample(s)
#> ✔ Calculating total sum abundance outliers at +/- 5 Sdev - excluding 0 sample(s…
#> 
#> ℹ Running sample data PCA outlier analysis at +/- 5 Sdev
#> ✔ Running sample data PCA outlier analysis at +/- 5 Sdev [20ms]
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
#>                summarise_raw    1.85  32.4
#>                   copy_layer    0.01   0.2
#>   extreme_sample_missingness    0.00   0.0
#>  extreme_feature_missingness    0.00   0.0
#>           sample_missingness    0.00   0.0
#>          total_sum_abundance    0.01   0.2
#>                summarise_pca    1.98  34.7
#>              summarise_final    1.61  28.2
#>                        total    5.71 100.0
#> ✔ Creating final QC dataset... [1.7s]
#> 
#> ℹ 'Omics QC Process Completed
#> ✔ 'Omics QC Process Completed [30ms]
```

## Export Omiprep

This will generate a folder containing an input and qc directory, each
containing flat text versions of data, features, samples, feature
summaries, and sample summaries.

``` r

# where to put the files
output_dir <- file.path(tempdir(), "output")

# run export
export(mydata, directory = output_dir, format = "omiprep")
#> Exporting in omiprep format to: 
#>      /tmp/RtmpN0SHrC/output

# view output directory files
files <- list.files(output_dir, full.names = TRUE, recursive = TRUE)
unname(sapply(files, function(path) {
  parts <- strsplit(path, .Platform$file.sep)[[1]]
  paste(tail(parts, 4), collapse = .Platform$file.sep)
}))
#>  [1] "output/omiprep_export_2026_06_25/input/config.yml"         
#>  [2] "output/omiprep_export_2026_06_25/input/data.tsv"           
#>  [3] "output/omiprep_export_2026_06_25/input/feature_summary.tsv"
#>  [4] "output/omiprep_export_2026_06_25/input/features.tsv"       
#>  [5] "output/omiprep_export_2026_06_25/input/sample_summary.tsv" 
#>  [6] "output/omiprep_export_2026_06_25/input/samples.tsv"        
#>  [7] "output/omiprep_export_2026_06_25/qc/config.yml"            
#>  [8] "output/omiprep_export_2026_06_25/qc/data.tsv"              
#>  [9] "output/omiprep_export_2026_06_25/qc/feature_summary.tsv"   
#> [10] "output/omiprep_export_2026_06_25/qc/feature_tree.RDS"      
#> [11] "output/omiprep_export_2026_06_25/qc/features.tsv"          
#> [12] "output/omiprep_export_2026_06_25/qc/sample_summary.tsv"    
#> [13] "output/omiprep_export_2026_06_25/qc/samples.tsv"           
#> [14] "output/omiprep_export_2026_06_25/qc/var_exp.tsv"
```

## Export to Comets format

This will generate a Comets formated excel sheet

``` r

# where to put the files
output_dir <- file.path(tempdir(), "output")

# run export
export(mydata, directory = output_dir, format = "comets")
#> Exporting data layer `qc` in comets format to /tmp/RtmpN0SHrC/output/omiprep_comets_export_2026_06_25.xlsx
#> Export complete.

# view output directory files
files <- list.files(output_dir, full.names = TRUE, recursive = TRUE)
w = grep("comets", files)
unname(sapply(files[w], function(path) {
  parts <- strsplit(path, .Platform$file.sep)[[1]]
  paste(tail(parts, 2), collapse = .Platform$file.sep)
}))
#> [1] "output/omiprep_comets_export_2026_06_25.xlsx"
```

## Export to Metaboanalyst format

This will generate a csv file with your feature
abundance/concentration/intensity data ready for upload to
(MetaboAnalyst)\[<https://www.metaboanalyst.ca>\].

``` r

# where to put the files
output_dir <- file.path(tempdir(), "output")

# run export
export(mydata, directory = output_dir, format = "metaboanalyst")
#> Exporting in export_metaboanalyst format to /tmp/RtmpN0SHrC/output
#> Exporting data layer `qc` in metaboanalyst format to /tmp/RtmpN0SHrC/output/omiprep_metaboanalyst_export_2026_06_25.csv

# view output directory files
files <- list.files(output_dir, full.names = TRUE, recursive = TRUE)
w = grep("metaboanalyst", files)
unname(sapply(files[w], function(path) {
  parts <- strsplit(path, .Platform$file.sep)[[1]]
  paste(tail(parts, 2), collapse = .Platform$file.sep)
}))
#> [1] "output/omiprep_metaboanalyst_export_2026_06_25.csv"
```
