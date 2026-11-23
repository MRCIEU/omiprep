# Importing In-House & Custom Data

``` r

library(omiprep)
```

## Who this vignette is for

`omiprep` ships dedicated readers for four commercial platforms —
[`read_metabolon()`](https://mrcieu.github.io/omiprep/reference/read_metabolon.md),
[`read_nightingale()`](https://mrcieu.github.io/omiprep/reference/read_nightingale.md),
[`read_olink()`](https://mrcieu.github.io/omiprep/reference/read_olink.md)
and
[`read_somalogic()`](https://mrcieu.github.io/omiprep/reference/read_somalogic.md).
Those readers exist only to translate a vendor’s particular spreadsheet
layout into a common structure. They are a convenience, **not** a
requirement.

If your data were generated in-house, come from a platform without a
dedicated reader, or have already been through an upstream pipeline, you
can build an `Omiprep` object directly from flat text files. Everything
downstream —
[`quality_control()`](https://mrcieu.github.io/omiprep/reference/quality_control.md),
[`summarise()`](https://mrcieu.github.io/omiprep/reference/summarise.md),
[`generate_report()`](https://mrcieu.github.io/omiprep/reference/generate_report.md),
[`export()`](https://mrcieu.github.io/omiprep/reference/export.md) —
then behaves identically, because all the readers do is produce the same
three tables you will assemble here by hand.

## The three tables

An `Omiprep` object is built from exactly three pieces:

| Piece | Shape | Requirement |
|----|----|----|
| `data` | numeric matrix, **samples in rows, features in columns** | `rownames` are sample ids, `colnames` are feature ids |
| `samples` | data.frame, one row per sample | must contain a column named `sample_id` |
| `features` | data.frame, one row per feature | must contain a column named `feature_id` |

Beyond those two id columns, both metadata tables are free-form: add as
many annotation columns as you like (batch, run order, sex, age,
pathway, assay panel, LOD), or none at all. `omiprep` carries them
through QC and reporting untouched.

### Are the column names strict?

**Yes.** `sample_id` and `feature_id` are matched exactly and are
**case-sensitive**. There is no automatic mapping or fuzzy matching:
`Sample_ID`, `sampleID` and `SampleId` are all rejected, with an error
naming the column it expected.

This is a deliberate design choice — silently guessing which of a user’s
columns holds the identifier is a good way to mis-assign an entire
dataset. Conforming is a one-line rename, shown below.

## A worked example from flat text files

The package bundles a small example dataset as three plain CSV files, in
exactly the form an in-house pipeline might produce them.

``` r

pkg <- "omiprep"
data_file     <- system.file("extdata", "dummy_data.csv",     package = pkg)
samples_file  <- system.file("extdata", "dummy_samples.csv",  package = pkg)
features_file <- system.file("extdata", "dummy_features.csv", package = pkg)
```

### 1) Read the measurement matrix

The measurement table has sample ids in the first column and feature ids
as headers. Read it with `row.names = 1` and coerce to a matrix.
`check.names = FALSE` keeps feature ids exactly as written — without it,
R prefixes purely numeric column names with `X`.

``` r

data <- as.matrix(read.csv(data_file, row.names = 1, check.names = FALSE))

dim(data)
#> [1] 125 253
data[1:4, 1:3]
#>           48719    43532      46639
#> sam_1 1.5376548       NA 0.26359902
#> sam_2 1.8382935       NA 0.07758324
#> sam_3        NA 1.147501         NA
#> sam_4 0.9071567 7.172319 0.41313486
```

For a tab-delimited file, use
[`read.delim()`](https://rdrr.io/r/utils/read.table.html) instead;
nothing else changes.

### 2) Read the sample and feature metadata

``` r

samples  <- read.csv(samples_file,  row.names = 1, check.names = FALSE)
features <- read.csv(features_file, row.names = 1, check.names = FALSE)

names(samples)
#>  [1] "sample_id"         "parent_sample_id"  "client_identifier"
#>  [4] "sex"               "age"               "bmi"              
#>  [7] "LC/MS Polar"       "LC/MS Neg"         "LC/MS Pos Early"  
#> [10] "LC/MS Pos Late"
names(features)
#>  [1] "feature_id"        "pathway_sortorder" "biochemical"      
#>  [4] "super_pathway"     "sub_pathway"       "comp_id"          
#>  [7] "platform"          "chemical_id"       "ri"               
#> [10] "mass"              "cas"               "pubchem"          
#> [13] "kegg"              "group_hmdb"
```

Both tables already carry the required id columns alongside their
annotation.

``` r

head(samples[, c("sample_id", "sex", "age", "bmi")], 3)
#>   sample_id  sex age   bmi
#> 1     sam_1 male  50 17.09
#> 2     sam_2 male  42 27.33
#> 3     sam_3 male  43 30.56
head(features[, c("feature_id", "biochemical", "super_pathway")], 3)
#>       feature_id   biochemical super_pathway
#> 48719      48719     X - 23145          <NA>
#> 43532      43532 veratric acid   Xenobiotics
#> 46639      46639     X - 12851          <NA>
```

### 3) Build the object

``` r

mydata <- Omiprep(data = data, samples = samples, features = features)

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
#>   Names  : sample_id, parent_sample_id, client_identifier, sex, age, bmi, LC/MS Polar, LC/MS Neg, LC/MS Pos Early, LC/MS Pos Late
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

That is the whole import. From here the object is indistinguishable from
one produced by
[`read_metabolon()`](https://mrcieu.github.io/omiprep/reference/read_metabolon.md).

## Renaming your own columns

If your files use different identifier names, rename them before
constructing the object:

``` r

# a metadata table as an in-house pipeline might label it
inhouse <- data.frame(
  SampleName = samples$sample_id,
  RunBatch   = samples$`LC/MS Polar`
)

names(inhouse)[names(inhouse) == "SampleName"] <- "sample_id"
names(inhouse)
#> [1] "sample_id" "RunBatch"
```

The helper
[`clean_names()`](https://mrcieu.github.io/omiprep/reference/clean_names.md)
is useful when your headers contain spaces, dots, dashes or percent
signs. It lowercases and standardises them, which turns a header like
`Sample ID` into `sample_id` directly:

``` r

clean_names(c("Sample ID", "Feature-Name.1", "Concentration %"))
#> [1] "sample_id"         "featurename_1"     "concentration_pct"
```

Note that
[`clean_names()`](https://mrcieu.github.io/omiprep/reference/clean_names.md)
lowercases everything, so apply it to the whole header vector rather
than to one column, and check the result before relying on it.

## Identifier types: a common pitfall

Matrix `dimnames` in R are *always* character. If your ids are numeric
and you read the metadata without specifying types,
[`read.csv()`](https://rdrr.io/r/utils/read.table.html) will type them
as integer, and the ids will not compare equal to the matrix names.

[`Omiprep()`](https://mrcieu.github.io/omiprep/reference/class_omiprep.md)
guards against this by coercing `sample_id` and `feature_id` to
character for you, so the naive read used above works. If you assemble
metadata by some other route, it is still worth being explicit:

``` r

features <- read.csv(features_file, colClasses = c(feature_id = "character"))
```

The ids must match the matrix names as a **set**; order does not matter.
The constructor reorders both metadata tables to follow the matrix,
carrying all annotation columns with them. Any id present in the matrix
but absent from the metadata is an error, and the message names the
offending ids.

## A minimal object

Annotation is optional. If you have nothing but a matrix, generate the
two id columns from its dimnames:

``` r

minimal <- Omiprep(
  data     = data,
  samples  = data.frame(sample_id  = rownames(data)),
  features = data.frame(feature_id = colnames(data))
)

dim(minimal@data)
#> [1] 125 253   1
```

## Running the pipeline

Everything downstream now works as documented elsewhere. Quality
control:

``` r

mydata <- quality_control(
  mydata,
  source_layer        = "input",
  sample_missingness  = 0.2,
  feature_missingness = 0.2,
  cores               = 2
)
#>                         step seconds   pct
#>                   validation    0.00   0.0
#>                summarise_raw    1.80  33.8
#>                   copy_layer    0.00   0.0
#>   extreme_sample_missingness    0.00   0.0
#>  extreme_feature_missingness    0.00   0.0
#>           sample_missingness    0.00   0.0
#>          total_sum_abundance    0.01   0.2
#>                summarise_pca    1.68  31.5
#>              summarise_final    1.60  30.0
#>                        total    5.33 100.1
```

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
#>   Names  : sample_id, parent_sample_id, client_identifier, sex, age, bmi, LC/MS Polar, LC/MS Neg, LC/MS Pos Early, LC/MS Pos Late, reason_excluded, excluded
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
#> extreme_sample_missingness        | 0
#> user_defined_sample_missingness   | 1
#> user_defined_sample_totalpeakarea | 0
#> user_defined_sample_pca_outlier   | 0
#> 
#>   Feature Exclusions:
#> Exclusion | Count
#> -----------------
#> user_excluded                    |  0
#> extreme_feature_missingness      |  5
#> user_defined_feature_missingness | 46
#> user_defined_feature_skewness    |  0
```

Note that `total_sum_abundance_sd` assumes features share a comparable
abundance scale. For in-house data that has already been normalised,
standardised, or combined across assays, set it to `NA` to skip that
step, as total abundance across such features is not interpretable.

## Where to go next

- [`vignette("quality_control")`](https://mrcieu.github.io/omiprep/articles/quality_control.md)
  — what each filtering step does and how to set thresholds
- [`vignette("skewness_qc")`](https://mrcieu.github.io/omiprep/articles/skewness_qc.md)
  — the optional skewness filter and when it is appropriate
- [`vignette("sample_summary")`](https://mrcieu.github.io/omiprep/articles/sample_summary.md)
  and
  [`vignette("feature_summary")`](https://mrcieu.github.io/omiprep/articles/feature_summary.md)
  — summaries without running QC
- [`vignette("generate_report")`](https://mrcieu.github.io/omiprep/articles/generate_report.md)
  — rendering an HTML or PDF QC report
- [`vignette("export")`](https://mrcieu.github.io/omiprep/articles/export.md)
  — writing results out for downstream tools
