# Export Data to \`MetaboAnalyst\` format

Export Data to \`MetaboAnalyst\` format

## Usage

``` r
export_metaboanalyst(omiprep, directory, layer = NULL, group_col = NULL)
```

## Arguments

- omiprep:

  A \`Omiprep\` object containing the data to be exported.

- directory:

  character, string specifying the path to the directory where the data
  should be written.

- layer:

  character, the name of the \`omiprep@data\` layer (3rd array
  dimension) to write out

- group_col:

  character, the column name in the \`omiprep@samples\` data identifying
  the group for one-factor analysis

## Value

the \`Omiprep\` object, invisibly, for use in pipes
