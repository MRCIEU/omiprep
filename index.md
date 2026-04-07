---
pagetitle: "metaboprep"
---

```{=html}
<div class="text-center py-4">
  <h1 class="display-5 fw-bold">metaboprep</h1>
  <p class="lead text-muted">Metabolomics &amp; proteomics data preparation and quality control pipeline for R</p>
  <div class="d-flex justify-content-center gap-2 flex-wrap mt-3">
    <a href="https://lifecycle.r-lib.org/articles/stages.html#experimental" target="_blank">
      <img src="https://img.shields.io/badge/lifecycle-experimental-orange.svg" alt="Lifecycle: experimental"/>
    </a>
  </div>
</div>
```

## Overview

`metaboprep` supports the full data-preparation workflow for untargeted and targeted omics data:

1. **Import** raw data from Metabolon, Nightingale Health, Olink, and SomaLogic platforms (Excel / flat-text)  
2. **Summarise** sample- and feature-level statistics  
3. **Filter** using a standard QC pipeline with user-defined thresholds  
4. **Report** results as an interactive HTML or PDF document  
5. **Export** cleaned data for downstream analysis

## Installation

```r
# install.packages("pak")
pak::pak("MRCIEU/metaboprep")
```

## Quick start

```r
library(metaboprep)

# 1. Read data
mydata <- read_metabolon(
  system.file("extdata", "metabolon_v1.1_example.xlsx", package = "metaboprep"),
  sheet             = "OrigScale",
  return_Metaboprep = TRUE
)

# 2. Run QC pipeline
mydata <- mydata |> quality_control(
  source_layer        = "input",
  sample_missingness  = 0.2,
  feature_missingness = 0.2,
  total_peak_area_sd  = 5,
  outlier_udist       = 5,
  outlier_treatment   = "leave_be"
)

# 3. Summarise
summary(mydata)

# 4. Generate HTML report
generate_report(mydata, output_dir = ".")
```

---

## Articles

```{=html}
<!-- ── Importing Data ───────────────────────────────────────────────── -->
<h3 class="mt-4 mb-3">Importing Data</h3>
<div class="row row-cols-1 row-cols-md-2 row-cols-lg-4 g-3 mb-4">

  <div class="col">
    <div class="card h-100 shadow-sm">
      <div class="card-body">
        <h5 class="card-title">
          <a href="articles/metabolon.html" class="stretched-link text-decoration-none">Metabolon</a>
        </h5>
        <p class="card-text text-muted small">Import untargeted metabolomics data from Metabolon Excel sheets.</p>
      </div>
      <div class="card-footer text-end">
        <a href="articles/metabolon.html" class="btn btn-sm btn-outline-primary">Read &rarr;</a>
      </div>
    </div>
  </div>

  <div class="col">
    <div class="card h-100 shadow-sm">
      <div class="card-body">
        <h5 class="card-title">
          <a href="articles/nightingale.html" class="stretched-link text-decoration-none">Nightingale Health</a>
        </h5>
        <p class="card-text text-muted small">Import NMR-based metabolomic data from Nightingale Health.</p>
      </div>
      <div class="card-footer text-end">
        <a href="articles/nightingale.html" class="btn btn-sm btn-outline-primary">Read &rarr;</a>
      </div>
    </div>
  </div>

  <div class="col">
    <div class="card h-100 shadow-sm">
      <div class="card-body">
        <h5 class="card-title">
          <a href="articles/olink.html" class="stretched-link text-decoration-none">Olink</a>
        </h5>
        <p class="card-text text-muted small">Import proximity extension assay proteomic data from Olink.</p>
      </div>
      <div class="card-footer text-end">
        <a href="articles/olink.html" class="btn btn-sm btn-outline-primary">Read &rarr;</a>
      </div>
    </div>
  </div>

  <div class="col">
    <div class="card h-100 shadow-sm">
      <div class="card-body">
        <h5 class="card-title">
          <a href="articles/somalogic.html" class="stretched-link text-decoration-none">SomaLogic</a>
        </h5>
        <p class="card-text text-muted small">Import aptamer-based proteomic data from SomaLogic SomaScan.</p>
      </div>
      <div class="card-footer text-end">
        <a href="articles/somalogic.html" class="btn btn-sm btn-outline-primary">Read &rarr;</a>
      </div>
    </div>
  </div>

</div>

<!-- ── Summaries ──────────────────────────────────────────────────── -->
<h3 class="mt-2 mb-3">Summaries &amp; QC</h3>
<div class="row row-cols-1 row-cols-md-2 row-cols-lg-3 g-3 mb-4">

  <div class="col">
    <div class="card h-100 shadow-sm">
      <div class="card-body">
        <h5 class="card-title">
          <a href="articles/sample_summary.html" class="stretched-link text-decoration-none">Sample Summary</a>
        </h5>
        <p class="card-text text-muted small">Compute per-sample statistics: missingness, total peak area, and PCA-based outlier detection.</p>
      </div>
      <div class="card-footer text-end">
        <a href="articles/sample_summary.html" class="btn btn-sm btn-outline-primary">Read &rarr;</a>
      </div>
    </div>
  </div>

  <div class="col">
    <div class="card h-100 shadow-sm">
      <div class="card-body">
        <h5 class="card-title">
          <a href="articles/feature_summary.html" class="stretched-link text-decoration-none">Feature Summary</a>
        </h5>
        <p class="card-text text-muted small">Compute per-feature statistics: missingness, variance, and independent feature trees.</p>
      </div>
      <div class="card-footer text-end">
        <a href="articles/feature_summary.html" class="btn btn-sm btn-outline-primary">Read &rarr;</a>
      </div>
    </div>
  </div>

  <div class="col">
    <div class="card h-100 shadow-sm">
      <div class="card-body">
        <h5 class="card-title">
          <a href="articles/quality_control.html" class="stretched-link text-decoration-none">QC Pipeline</a>
        </h5>
        <p class="card-text text-muted small">Run the full quality control pipeline with configurable thresholds for missingness, outliers, and more.</p>
      </div>
      <div class="card-footer text-end">
        <a href="articles/quality_control.html" class="btn btn-sm btn-outline-primary">Read &rarr;</a>
      </div>
    </div>
  </div>

</div>

<!-- ── Reports & Export ───────────────────────────────────────────── -->
<h3 class="mt-2 mb-3">Reports &amp; Export</h3>
<div class="row row-cols-1 row-cols-md-2 row-cols-lg-3 g-3 mb-4">

  <div class="col">
    <div class="card h-100 shadow-sm">
      <div class="card-body">
        <h5 class="card-title">
          <a href="articles/generate_report.html" class="stretched-link text-decoration-none">Generate HTML / PDF Report</a>
        </h5>
        <p class="card-text text-muted small">Produce a fully annotated, interactive QC report in HTML or PDF format.</p>
      </div>
      <div class="card-footer text-end">
        <a href="articles/generate_report.html" class="btn btn-sm btn-outline-primary">Read &rarr;</a>
      </div>
    </div>
  </div>

  <div class="col">
    <div class="card h-100 shadow-sm">
      <div class="card-body">
        <h5 class="card-title">
          <a href="articles/export.html" class="stretched-link text-decoration-none">Export Data</a>
        </h5>
        <p class="card-text text-muted small">Export processed data and summary tables to Excel or tab-delimited flat files.</p>
      </div>
      <div class="card-footer text-end">
        <a href="articles/export.html" class="btn btn-sm btn-outline-primary">Read &rarr;</a>
      </div>
    </div>
  </div>

  <div class="col">
    <div class="card h-100 shadow-sm">
      <div class="card-body">
        <h5 class="card-title">
          <a href="articles/batch_normalise.html" class="stretched-link text-decoration-none">Batch Normalisation</a>
        </h5>
        <p class="card-text text-muted small">Correct for run-order and batch effects using quantile or rank-based normalisation.</p>
      </div>
      <div class="card-footer text-end">
        <a href="articles/batch_normalise.html" class="btn btn-sm btn-outline-primary">Read &rarr;</a>
      </div>
    </div>
  </div>

</div>
```
