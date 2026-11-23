# Metaboprep 1 pipeline

This function runs the original metaboprep1 pipeline using the old
parameter file input format. The function requires access to the
internet as the old package (default github tag \`v1.0.1\`) will be
dynamically downloaded and used to process the data.

## Usage

``` r
run_metaboprep1(parameter_file, gitcommit = "v1.0.1", attempt_report = FALSE)
```

## Arguments

- parameter_file:

  character, full file path to the metaboprep 1 parameter file

- gitcommit:

  character, Github commit - default pinned to the metaboprep 1 release
  \`v1.0.1\`

- attempt_report:

  logical, whether to attempt metaboprep1 report generation.
  Default=FALSE as this can lead to errors on some operating systems.
