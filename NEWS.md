# omiprep 0.9.0

First public, peer-review release. `omiprep` is the successor to **metaboprep**,
generalised from Nightingale/Metabolon metabolomics to mixed metabolomics and
proteomics platforms.

## Data import

* Unified readers into a single `Omiprep` object: `read_metabolon()` (v1.1, v1.2,
  and v2 multi-tab formats), `read_nightingale()` (single- and multi-sheet
  formats), `read_olink()`, and `read_somalogic()`.
* Synthetic example datasets bundled for the Metabolon and Nightingale formats.

## Quality control

* `quality_control()` orchestrates sample/feature missingness filtering,
  total-sum-abundance outliers, PC-based sample outliers, feature-tree reduction,
  and winsorisation, recording all exclusions on the object.
* Building blocks exported individually: `missingness()`, `outlier_detection()`,
  `pc_and_outliers()`, `total_sum_abundance()`, `tree_and_independent_features()`,
  `feature_summary()`, `sample_summary()`, `summarise()`.

## Normalisation & layers

* `batch_normalise()` for platform / run-mode batch normalisation.
* Multi-layer `Omiprep@data` array with `add_layer()`.

## Reporting & export

* `generate_report()` renders HTML/PDF QC reports
  (`available_report_templates()`).
* `export_omiprep()`, `export_metaboanalyst()`, and `export_comets()` for
  downstream tools.
* `shiny_app()` interactive explorer.

## Compatibility

* `run_metaboprep1()` preserves the original metaboprep workflow.
