library(testthat)

# Integration tests for the platform readers (R/read_*.R) against the example
# datasets bundled in inst/extdata. These pin:
#   - the universal reader contract (numeric named data matrix; sample_id /
#     feature_id metadata columns; data dimnames covered by the metadata),
#   - exact parsed dimensions per example file (regression guard against
#     parser drift), and
#   - that the returned (data, samples, features) triple builds a valid
#     Omiprep object.
#
# Readers prompt via readline() when sheet = NULL, so every call passes an
# explicit sheet/argument set (taken from each reader's @examples).

ex_path <- function(file) system.file("extdata", file, package = "omiprep")

# Universal contract every reader's list return must satisfy.
expect_reader_contract <- function(r) {
  expect_type(r, "list")
  expect_true(all(c("data", "samples", "features") %in% names(r)))

  expect_true(is.matrix(r$data))
  expect_true(is.numeric(r$data))
  expect_false(is.null(rownames(r$data)))
  expect_false(is.null(colnames(r$data)))
  expect_false(anyNA(rownames(r$data)))
  expect_false(anyNA(colnames(r$data)))

  expect_true("sample_id"  %in% names(r$samples))
  expect_true("feature_id" %in% names(r$features))

  # every data axis label is described in the metadata tables
  expect_true(all(colnames(r$data) %in% r$features$feature_id))
  expect_true(all(rownames(r$data) %in% r$samples$sample_id))
}

# The triple must build a valid Omiprep whose input layer is the data matrix.
expect_builds_omiprep <- function(r) {
  o <- suppressMessages(Omiprep(r$data, r$samples, r$features))
  expect_s3_class(o, "omiprep::Omiprep")
  expect_identical(dim(o@data)[1:2], dim(r$data))
  expect_identical(dimnames(o@data)[[3]], "input")
  expect_equal(o@data[, , "input"], r$data)
}

# ---- Metabolon (three on-disk format variants) ------------------------------

test_that("read_metabolon parses the v1.1 single-sheet format", {
  fp <- ex_path("metabolon_v1.1_example.xlsx")
  skip_if(fp == "", "example file not bundled")
  r <- suppressWarnings(read_metabolon(fp, sheet = 2))

  expect_reader_contract(r)
  expect_identical(dim(r$data), c(100L, 100L))
  expect_identical(r$samples$sample_id, rownames(r$data))   # reader asserts this
  expect_identical(r$features$feature_id, colnames(r$data))
  expect_builds_omiprep(r)
})

test_that("read_metabolon parses the v1.2 named-sheet format", {
  fp <- ex_path("metabolon_v1.2_example.xlsx")
  skip_if(fp == "", "example file not bundled")
  r <- suppressWarnings(read_metabolon(fp, sheet = "OrigScale"))

  expect_reader_contract(r)
  expect_identical(dim(r$data), c(100L, 104L))
  expect_builds_omiprep(r)
})

test_that("read_metabolon parses the v2 multi-tab format", {
  fp <- ex_path("metabolon_v2_example.xlsx")
  skip_if(fp == "", "example file not bundled")
  r <- suppressWarnings(read_metabolon(
    fp,
    sheet          = "Batch-normalized Data",
    feature_sheet  = "Chemical Annotation",
    feature_id_col = "CHEM_ID",
    sample_sheet   = "Sample Meta Data",
    sample_id_col  = "PARENT_SAMPLE_NAME"
  ))

  expect_reader_contract(r)
  expect_identical(dim(r$data), c(10L, 8L))
  expect_builds_omiprep(r)
})

test_that("read_metabolon requires the metadata args for the multi-tab format", {
  fp <- ex_path("metabolon_v2_example.xlsx")
  skip_if(fp == "", "example file not bundled")
  # sheet present but feature/sample sheet+id args omitted
  expect_error(
    suppressWarnings(read_metabolon(fp, sheet = "Batch-normalized Data")),
    "feature_sheet|feature_id_col|sample_sheet|sample_id_col"
  )
})

test_that("read_metabolon errors on a feature_sheet that is not a tab", {
  fp <- ex_path("metabolon_v2_example.xlsx")
  skip_if(fp == "", "example file not bundled")
  expect_error(
    suppressWarnings(read_metabolon(
      fp,
      sheet          = "Batch-normalized Data",
      feature_sheet  = "Nonexistent Sheet",
      feature_id_col = "CHEM_ID",
      sample_sheet   = "Sample Meta Data",
      sample_id_col  = "PARENT_SAMPLE_NAME"
    )),
    "not found"
  )
})

test_that("read_metabolon rejects a non-excel extension", {
  expect_error(read_metabolon("not_a_spreadsheet.txt"),
               "\\.xls|excel")
})

test_that("read_metabolon carries ids as rownames on samples and features", {
  fp <- ex_path("metabolon_v1.1_example.xlsx")
  skip_if(fp == "", "example file not bundled")

  r <- suppressWarnings(read_metabolon(fp, sheet = 2))
  expect_identical(rownames(r$samples), r$samples$sample_id)
  expect_identical(rownames(r$features), r$features$feature_id)

  # and the ids survive as rownames through Omiprep construction
  o <- suppressWarnings(read_metabolon(fp, sheet = 2, return_Omiprep = TRUE))
  expect_identical(rownames(o@samples), o@samples$sample_id)
  expect_identical(rownames(o@features), o@features$feature_id)
})

# ---- Nightingale (single-sheet and multi-sheet formats) ---------------------

test_that("read_nightingale parses the v1 example", {
  fp <- ex_path("nightingale_v1_example.xlsx")
  skip_if(fp == "", "example file not bundled")
  r <- suppressWarnings(read_nightingale(fp))

  expect_reader_contract(r)
  expect_identical(dim(r$data), c(50L, 12L))
  expect_identical(r$features$feature_id, colnames(r$data))
  expect_builds_omiprep(r)
})

test_that("read_nightingale parses the v2 example", {
  fp <- ex_path("nightingale_v2_example.xlsx")
  skip_if(fp == "", "example file not bundled")
  r <- suppressWarnings(read_nightingale(fp))

  expect_reader_contract(r)
  expect_identical(dim(r$data), c(50L, 15L))
  expect_builds_omiprep(r)
})

test_that("read_nightingale rejects a non-excel extension", {
  expect_error(read_nightingale("not_a_spreadsheet.txt"),
               "\\.xls|excel")
})

# ---- Olink ------------------------------------------------------------------

test_that("read_olink parses the NPX example and splits samples vs controls", {
  skip_if_not_installed("OlinkAnalyze")
  fp <- ex_path("olink_v1_example.txt")
  skip_if(fp == "", "example file not bundled")
  r <- suppressWarnings(suppressMessages(read_olink(fp)))

  expect_reader_contract(r)
  expect_identical(dim(r$data), c(90L, 100L))
  expect_true(all(c("controls", "control_metadata") %in% names(r)))
  # controls mirror data's numeric-matrix shape (samples x assays)
  expect_true(is.matrix(r$controls))
  expect_true(is.numeric(r$controls))
  # samples and controls are disjoint sets of ids
  expect_length(intersect(rownames(r$data), rownames(r$controls)), 0)
  expect_builds_omiprep(r)

  # regression: the per-assay QC_Warning column must NOT duplicate samples.
  # One row per sample, no duplicate ids, and QC summarised to sample level.
  expect_identical(nrow(r$samples), nrow(r$data))
  expect_false(any(duplicated(r$samples$sample_id)))
  expect_setequal(r$samples$sample_id, rownames(r$data))
  expect_true(all(c("qc_n_warnings", "qc_any_warning") %in% names(r$samples)))
  expect_type(r$samples$qc_n_warnings, "integer")
  expect_type(r$samples$qc_any_warning, "logical")
  # Sample_00001 had both Pass and Warning assay rows -> flagged at sample level
  expect_true(r$samples$qc_any_warning[r$samples$sample_id == "Sample_00001"])
})

test_that("read_olink rejects an unsupported extension", {
  expect_error(read_olink("data.pdf"), "Olink|extension")
})

# ---- SomaLogic --------------------------------------------------------------

test_that("read_somalogic parses the adat example and splits samples vs controls", {
  skip_if_not_installed("SomaDataIO")
  fp <- ex_path("somalogic_v1_example.adat")
  skip_if(fp == "", "example file not bundled")
  r <- suppressWarnings(suppressMessages(read_somalogic(fp)))

  expect_reader_contract(r)
  expect_identical(dim(r$data), c(9L, 5284L))
  expect_true(all(grepl("^seq\\.", colnames(r$data))))
  expect_true(all(c("controls", "control_metadata") %in% names(r)))
  expect_builds_omiprep(r)
})

test_that("read_somalogic rejects a non-adat extension", {
  expect_error(read_somalogic("data.txt"), "adat")
})
