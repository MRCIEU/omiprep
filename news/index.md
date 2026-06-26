# Changelog

## omiprep 0.9.0

First public, peer-review release. `omiprep` is the successor to
**metaboprep**, generalised from Nightingale/Metabolon metabolomics to
mixed metabolomics and proteomics platforms.

### Data import

- Unified readers into a single `Omiprep` object:
  [`read_metabolon()`](https://mrcieu.github.io/omiprep/reference/read_metabolon.md)
  (v1.1, v1.2, and v2 multi-tab formats),
  [`read_nightingale()`](https://mrcieu.github.io/omiprep/reference/read_nightingale.md)
  (single- and multi-sheet formats),
  [`read_olink()`](https://mrcieu.github.io/omiprep/reference/read_olink.md),
  and
  [`read_somalogic()`](https://mrcieu.github.io/omiprep/reference/read_somalogic.md).
- Synthetic example datasets bundled for the Metabolon and Nightingale
  formats.

### Quality control

- [`quality_control()`](https://mrcieu.github.io/omiprep/reference/quality_control.md)
  orchestrates sample/feature missingness filtering, total-sum-abundance
  outliers, PC-based sample outliers, feature-tree reduction, and
  winsorisation, recording all exclusions on the object.
- Building blocks exported individually:
  [`missingness()`](https://mrcieu.github.io/omiprep/reference/missingness.md),
  [`outlier_detection()`](https://mrcieu.github.io/omiprep/reference/outlier_detection.md),
  [`pc_and_outliers()`](https://mrcieu.github.io/omiprep/reference/pc_and_outliers.md),
  [`total_sum_abundance()`](https://mrcieu.github.io/omiprep/reference/total_sum_abundance.md),
  [`tree_and_independent_features()`](https://mrcieu.github.io/omiprep/reference/tree_and_independent_features.md),
  [`feature_summary()`](https://mrcieu.github.io/omiprep/reference/feature_summary.md),
  [`sample_summary()`](https://mrcieu.github.io/omiprep/reference/sample_summary.md),
  [`summarise()`](https://mrcieu.github.io/omiprep/reference/summarise.md).

### Normalisation & layers

- [`batch_normalise()`](https://mrcieu.github.io/omiprep/reference/batch_normalise.md)
  for platform / run-mode batch normalisation.
- Multi-layer `Omiprep@data` array with
  [`add_layer()`](https://mrcieu.github.io/omiprep/reference/add_layer.md).

### Reporting & export

- [`generate_report()`](https://mrcieu.github.io/omiprep/reference/generate_report.md)
  renders HTML/PDF QC reports
  ([`available_report_templates()`](https://mrcieu.github.io/omiprep/reference/available_report_templates.md)).
- [`export_omiprep()`](https://mrcieu.github.io/omiprep/reference/export_omiprep.md),
  [`export_metaboanalyst()`](https://mrcieu.github.io/omiprep/reference/export_metaboanalyst.md),
  and
  [`export_comets()`](https://mrcieu.github.io/omiprep/reference/export_comets.md)
  for downstream tools.
- [`shiny_app()`](https://mrcieu.github.io/omiprep/reference/shiny_app.md)
  interactive explorer.

### Compatibility

- [`run_metaboprep1()`](https://mrcieu.github.io/omiprep/reference/run_metaboprep1.md)
  preserves the original metaboprep workflow.
