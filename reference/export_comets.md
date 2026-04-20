# Export Data to \`COMETS\` format

Export Data to \`COMETS\` format

## Usage

``` r
export_comets(omiprep, directory, layer = NULL)
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

## Value

the \`Omiprep\` object, invisibly, for use in pipes
