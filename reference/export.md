# Export Data from a Omiprep Object

Exports all data from a \`Omiprep\` object to a structured directory
format. For each data layer, the function creates a subdirectory
containing: - the primary data matrix (\`data.tsv\`), - associated
feature and sample metadata (\`features.tsv\`, \`samples.tsv\`), -
feature and sample summaries (if present, \`feature_summary.tsv\`,
\`sample_summary.tsv\`), - a serialized feature tree (if present), - and
a \`config.yml\` file with additional metadata and processing
parameters.

## Usage

``` r
export(omiprep, directory, format = "omiprep", ...)
```

## Arguments

- omiprep:

  A \`Omiprep\` object containing the data to be exported.

- directory:

  character, string specifying the path to the directory where the data
  should be written.

- format:

  character, string specifying the format of the exported data - one of
  "omiprep", "comets", or "metaboanalyst".

- ...:

  Arguments passed on to
  [`export_comets`](https://mrcieu.github.io/omiprep/reference/export_comets.md),
  [`export_metaboanalyst`](https://mrcieu.github.io/omiprep/reference/export_metaboanalyst.md)

  `layer`

  :   character, the name of the \`omiprep@data\` layer (3rd array
      dimension) to write out

  `group_col`

  :   character, the column name in the \`omiprep@samples\` data
      identifying the group for one-factor analysis

## Value

the \`Omiprep\` object, invisibly, for use in pipes
