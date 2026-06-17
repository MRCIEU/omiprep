# Read and Process SomaLogic adat file

This function reads and processes a commercial SomaLogic \`.adat\` file.
It extracts RFU (Relative Florecent Units) data for samples and
controls, along with their respective metadata and feature (protein)
metadata. The function returns a structured list suitable for further
analysis.

## Usage

``` r
read_somalogic(filepath, return_Omiprep = FALSE)
```

## Arguments

- filepath:

  A string specifying the path to the SomaLogic \`.adat\` file.

- return_Omiprep:

  logical, if TRUE return a Omiprep object, if FALSE (default) return a
  list.

## Value

A omiprep object or a named list with the following elements:

- data:

  A matrix of RFU values for experimental samples, with \`SampleId\` as
  row names and \`SeqId\` columns as column names.

- samples:

  A data frame containing metadata for experimental samples, with
  \`sample_id\` renamed from \`SampleId\`.

- features:

  A data frame containing feature-level metadata, including a newly
  created \`feature_id\` column derived from \`SeqId\`.

- controls:

  A matrix of RFU values for control samples, specifically "Calibrator"
  samples, with \`SampleId\` as row names and \`SeqId\` columns as
  column names. If duplicate control \`SampleId\` values are present,
  control \`SampleId\` is replaced with \`paste0(SampleId, "\_",
  SlideId)\`.

- control_metadata:

  A data frame containing metadata for control samples, specifically
  "Calibrator" samples, with \`sample_id\` renamed from \`SampleId\`.
