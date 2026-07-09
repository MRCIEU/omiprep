## Import

Two ways to load 'omics data into omiprep.

### 1 · Base R data structures

\# build matrices yourself samples \<- data.frame('sample_id' = ...)
features \<- data.frame('feature_id' = ...) data \<- matrix(input_dat,
nrow = nrow(features), ncol = nrow(samples), dimnames =
list(samples\$sample_id, features\$feature_id))

### 2 · Commercial format read functions

read_metabolon(filepath, sheet, feature_sheet, feature_id_col,
sample_sheet, sample_id_col, return_Omiprep) read_nightingale(filepath,
return_Omiprep) read_olink(filepath, return_Omiprep)
read_somalogic(filepath, return_Omiprep)

## Quality control pipeline

Wrapper that runs the full omiprep QC pipeline on a layer.

quality_control(omiprep, source_layer = "input", sample_missingness =
0.2, feature_missingness = 0.2, feature_skewness_threshold = NULL,
feature_skewness_direction = "left", total_sum_abundance_sd = 5,
outlier_udist = 5, outlier_treatment = "leave_be", winsorize_quantile =
1.0, tree_cut_height = 0.5, feature_selection = "max_var_exp",
pc_outlier_sd = 5, max_num_pcs = 10, sample_ids = NULL, feature_ids =
NULL, features_exclude_but_keep = NULL, cores = NULL, fast = FALSE)

## Omiprep object

Construct the S7 container for omics data, samples, and features.

Omiprep(data = data, samples = samples, features = features) \# inspect
& layer-stack helpers summary(omiprep) add_layer(current, layer,
layer_name, force = FALSE) batch_normalise(omiprep, run_mode_col,
run_mode_colmap, source_layer = "input", dest_layer =
"batch_normalised")

## Reporting

Render an HTML or PDF QC report for an Omiprep object.

generate_report(omiprep, output_dir, output_filename = NULL, project =
"Project", format = "pdf", template = "qc_report")

## Summarise

Per-sample, per-feature, or combined summary statistics.

feature_summary(omiprep, source_layer, outlier_udist, tree_cut_height,
feature_selection, sample_ids, feature_ids, features_exclude, output,
cores, fast) sample_summary(omiprep, source_layer, outlier_udist,
sample_ids, feature_ids, output) summarise(omiprep, source_layer,
outlier_udist, tree_cut_height, feature_selection, sample_ids,
feature_ids, features_exclude, output, cores, fast)

pc_and_outliers(omiprep, source_layer, sample_ids, feature_ids)
tree_and_independent_features(data, tree_cut_height = 0.5,
features_exclude = NULL, feature_selection = "max_var_exp", cores =
NULL, fast = FALSE) feature_skewness(data, threshold = NULL, direction =
"left")

## Export

Write Omiprep data to disk as flat files, or in COMETS / MetaboAnalyst
format.

export(omiprep, directory, format = "omiprep", ...) \# format: "omiprep"
\| "comets" \| "metaboanalyst" export_omiprep(omiprep, directory, ...)
export_comets(omiprep, directory, layer = NULL)
export_metaboanalyst(omiprep, directory, layer = NULL, group_col = NULL)

\# typical pipe mydata \|\> quality_control(source_layer = "input") \|\>
generate_report(output_dir = "out/", format = "html") \|\>
export("out/", format = "omiprep")

## Pipeline

High-level flow from raw vendor file to QC'd, reportable Omiprep object.

**1 · Read data**

- read_metabolon / nightingale
- read_olink / somalogic

**2 · Pre-processing**

- Normalisation
- batch_normalise()

**3 · Summary statistics  
(raw data)**

**4 · Data filtering**

- Sample / feature missingness
- Feature skewness
- Total sum abundance outliers
- PCA outliers (max_num_pcs)

**5 · Summary statistics  
(filtered data)**

**6 · Generate report & export**

- generate_report()
- export()

`function()` Outlier treatment: `leave_be` · `turn_NA` · `winsorize`
