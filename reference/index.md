# Package index

## Omiprep Object

The main S7 class behind the omiprep R package along with class helper
functions. A container for measurement, sample, and feature data.

- [`Omiprep()`](https://mrcieu.github.io/omiprep/reference/class_omiprep.md)
  : Omiprep Object
- [`summary.Omiprep`](https://mrcieu.github.io/omiprep/reference/summary.Omiprep.md)
  : Summary Method for Omiprep Object
- [`add_layer()`](https://mrcieu.github.io/omiprep/reference/add_layer.md)
  : Add a Layer of Data (internal use)

## Importing & Exporting Data

Functions importing data from a variety of data sources and formats.

- [`read_metabolon()`](https://mrcieu.github.io/omiprep/reference/read_metabolon.md)
  : Read Metabolon Data
- [`read_nightingale()`](https://mrcieu.github.io/omiprep/reference/read_nightingale.md)
  : Read Nightingale Data (format 1)
- [`read_olink()`](https://mrcieu.github.io/omiprep/reference/read_olink.md)
  : Read and Process Olink NPX Data File
- [`read_somalogic()`](https://mrcieu.github.io/omiprep/reference/read_somalogic.md)
  : Read and Process SomaLogic adat file
- [`available_data_formats()`](https://mrcieu.github.io/omiprep/reference/available_data_formats.md)
  : List Available Data Formats
- [`export()`](https://mrcieu.github.io/omiprep/reference/export.md) :
  Export Data from a Omiprep Object
- [`export_comets()`](https://mrcieu.github.io/omiprep/reference/export_comets.md)
  : Export Data to \`COMETS\` format
- [`export_metaboanalyst()`](https://mrcieu.github.io/omiprep/reference/export_metaboanalyst.md)
  : Export Data to \`MetaboAnalyst\` format
- [`export_omiprep()`](https://mrcieu.github.io/omiprep/reference/export_omiprep.md)
  : Export Data to \`Omiprep\` format

## Summary functions

Functions to summarise data.

- [`summarise()`](https://mrcieu.github.io/omiprep/reference/summarise.md)
  : Summary Statistics
- [`feature_summary()`](https://mrcieu.github.io/omiprep/reference/feature_summary.md)
  : Feature Summary Statistics
- [`sample_summary()`](https://mrcieu.github.io/omiprep/reference/sample_summary.md)
  : Sample Summary Statistics
- [`pc_and_outliers()`](https://mrcieu.github.io/omiprep/reference/pc_and_outliers.md)
  : Principal Component Analysis
- [`tree_and_independent_features()`](https://mrcieu.github.io/omiprep/reference/tree_and_independent_features.md)
  : Identify Independent Features in a Numeric Matrix
- [`feature_skewness()`](https://mrcieu.github.io/omiprep/reference/feature_skewness.md)
  : Estimate Feature Skewness

## Quality control & Reporting

Functions to run the quality control pipeline and generate a report.

- [`quality_control()`](https://mrcieu.github.io/omiprep/reference/quality_control.md)
  : Omics Quality Control
- [`available_report_templates()`](https://mrcieu.github.io/omiprep/reference/available_report_templates.md)
  : List Available Report Templates
- [`generate_report()`](https://mrcieu.github.io/omiprep/reference/generate_report.md)
  : Generate Output Report
- [`run_metaboprep1()`](https://mrcieu.github.io/omiprep/reference/run_metaboprep1.md)
  : Metaboprep 1 pipeline
- [`shiny_app()`](https://mrcieu.github.io/omiprep/reference/shiny_app.md)
  : Omiprep Shiny App

## Other processing tools

Stand alone data processing or filtering tools.

- [`batch_normalise()`](https://mrcieu.github.io/omiprep/reference/batch_normalise.md)
  : Batch Normalisation

## Helper functions

Helper functions, mainly used internally.

- [`feature_describe()`](https://mrcieu.github.io/omiprep/reference/feature_describe.md)
  : Summary Statistics for Features
- [`missingness()`](https://mrcieu.github.io/omiprep/reference/missingness.md)
  : Estimate Missingness
- [`outlier_detection()`](https://mrcieu.github.io/omiprep/reference/outlier_detection.md)
  : Identify indexes of outliers in data
- [`outliers()`](https://mrcieu.github.io/omiprep/reference/outliers.md)
  : Identify Outliers
- [`total_peak_area()`](https://mrcieu.github.io/omiprep/reference/total_peak_area.md)
  : Estimates total peak abundance
- [`continuous_power_plot()`](https://mrcieu.github.io/omiprep/reference/continuous_power_plot.md)
  : continuous trait power analysis plot
- [`multivariate_anova()`](https://mrcieu.github.io/omiprep/reference/multivariate_anova.md)
  : multivariate analysis
- [`cramerV()`](https://mrcieu.github.io/omiprep/reference/cramerV.md) :
  Cramer's V (phi)
- [`eval.power.binary.imbalanced()`](https://mrcieu.github.io/omiprep/reference/eval.power.binary.imbalanced.md)
  : Estimate power for a binary variable in an imbalanced design
- [`eval.power.cont()`](https://mrcieu.github.io/omiprep/reference/eval.power.cont.md)
  : estimate power for continuous variable
- [`find.PA.effect.sizes.2.sim()`](https://mrcieu.github.io/omiprep/reference/find.PA.effect.sizes.2.sim.md)
  : identify effect sizes
- [`find.cont.effect.sizes.2.sim()`](https://mrcieu.github.io/omiprep/reference/find.cont.effect.sizes.2.sim.md)
  : identify continuos trait effect sizes
- [`imbalanced_power_plot()`](https://mrcieu.github.io/omiprep/reference/imbalanced_power_plot.md)
  : binary trait imbalanced design power analysis plot
- [`variable_by_factor()`](https://mrcieu.github.io/omiprep/reference/variable_by_factor.md)
  : ggplot2 violin plot
- [`clean_names()`](https://mrcieu.github.io/omiprep/reference/clean_names.md)
  : Standardize Column or Feature Names
