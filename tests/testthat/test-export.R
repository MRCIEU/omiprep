library(testthat)

# Tests for the export methods (R/export.R): the dispatching `export()` wrapper
# and the three back-ends export_omiprep / export_comets / export_metaboanalyst.
# Each test writes to a throwaway temp directory and reads the output back.

# ---- shared fixture ---------------------------------------------------------

make_export_omiprep <- function(n_samples = 6, n_features = 4) {
  data <- matrix(as.numeric(seq_len(n_samples * n_features)),
                 nrow = n_samples, ncol = n_features)
  rownames(data) <- paste0("Sample_", seq_len(n_samples))
  colnames(data) <- paste0("Met_", seq_len(n_features))
  samples <- data.frame(
    sample_id = rownames(data),
    sex       = factor(rep(c("M", "F"), length.out = n_samples)),
    age       = seq(40, by = 2, length.out = n_samples)
  )
  features <- data.frame(feature_id = colnames(data))
  list(
    data = data,
    m    = suppressMessages(Omiprep(data, samples, features))
  )
}

# ---- export() dispatch ------------------------------------------------------

test_that("export() defaults to omiprep format and returns the object invisibly", {
  fx  <- make_export_omiprep()
  dir <- withr::local_tempdir()
  res <- expect_invisible(suppressMessages(export(fx$m, dir)))
  expect_s3_class(res, "omiprep::Omiprep")
  # the omiprep-format export directory was created
  expect_true(any(grepl("^omiprep_export_", list.files(dir))))
})

test_that("export() rejects an unknown format", {
  fx  <- make_export_omiprep()
  dir <- withr::local_tempdir()
  expect_error(suppressMessages(export(fx$m, dir, format = "not_a_format")),
               "should be one of|arg")
})

# ---- export_omiprep ---------------------------------------------------------

test_that("export_omiprep writes one subdir per layer with the expected files", {
  fx  <- make_export_omiprep()
  dir <- withr::local_tempdir()
  suppressMessages(export_omiprep(fx$m, dir))

  export_dir <- list.files(dir, pattern = "^omiprep_export_", full.names = TRUE)
  expect_length(export_dir, 1)
  # input-only object -> a single "input" layer subdir
  expect_identical(list.files(export_dir), "input")

  written <- list.files(file.path(export_dir, "input"))
  expect_true(all(c("data.tsv", "samples.tsv", "features.tsv", "config.yml")
                  %in% written))
})

test_that("export_omiprep data.tsv round-trips the input matrix", {
  fx  <- make_export_omiprep()
  dir <- withr::local_tempdir()
  suppressMessages(export_omiprep(fx$m, dir))
  data_path <- list.files(dir, pattern = "data.tsv$",
                          recursive = TRUE, full.names = TRUE)

  back <- read.table(data_path, sep = "\t", header = TRUE,
                     check.names = FALSE, stringsAsFactors = FALSE)
  expect_identical(colnames(back), c("sample_id", colnames(fx$data)))
  expect_identical(back[["sample_id"]], rownames(fx$data))

  values <- as.matrix(back[, -1])
  rownames(values) <- back[["sample_id"]]
  expect_equal(values, fx$data)
})

test_that("export_omiprep config records correct sample/feature counts", {
  fx  <- make_export_omiprep()
  dir <- withr::local_tempdir()
  suppressMessages(export_omiprep(fx$m, dir))
  cfg_path <- list.files(dir, pattern = "config.yml$",
                         recursive = TRUE, full.names = TRUE)
  cfg <- yaml::read_yaml(cfg_path)

  expect_identical(cfg$layer_name, "input")
  expect_identical(cfg$total_n_features, ncol(fx$data))
  expect_identical(cfg$total_n_samples, nrow(fx$data))
  expect_identical(cfg$included_n_features, ncol(fx$data))
  expect_identical(cfg$included_n_samples, nrow(fx$data))
})

# ---- export_comets ----------------------------------------------------------

test_that("export_comets writes an xlsx with the five COMETS sheets", {
  fx  <- make_export_omiprep()
  dir <- withr::local_tempdir()
  fp  <- suppressMessages(export_comets(fx$m, dir))

  expect_true(file.exists(fp))
  expect_match(fp, "omiprep_comets_export_.*\\.xlsx$")
  expect_setequal(
    openxlsx::getSheetNames(fp),
    c("Metabolites", "SubjectMetabolites", "SubjectData", "VarMap", "Models")
  )
})

test_that("export_comets SubjectMetabolites sheet round-trips the data", {
  fx  <- make_export_omiprep()
  dir <- withr::local_tempdir()
  fp  <- suppressMessages(export_comets(fx$m, dir))

  sheet <- openxlsx::read.xlsx(fp, sheet = "SubjectMetabolites")
  expect_identical(sheet[["SAMPLE_ID"]], rownames(fx$data))
  values <- as.matrix(sheet[, colnames(fx$data)])
  rownames(values) <- sheet[["SAMPLE_ID"]]
  expect_equal(values, fx$data)
})

test_that("export_comets errors on a layer that does not exist", {
  fx  <- make_export_omiprep()
  dir <- withr::local_tempdir()
  expect_error(suppressMessages(export_comets(fx$m, dir, layer = "qc")),
               "not found")
})

# ---- export_metaboanalyst ---------------------------------------------------

test_that("export_metaboanalyst writes a csv with a default group column", {
  fx  <- make_export_omiprep()
  dir <- withr::local_tempdir()
  suppressMessages(export_metaboanalyst(fx$m, dir))

  fp <- list.files(dir, pattern = "^omiprep_metaboanalyst_export_.*\\.csv$",
                   full.names = TRUE)
  expect_length(fp, 1)

  back <- read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE)
  expect_true(all(c("PatientID", "group") %in% names(back)))
  expect_setequal(back[["PatientID"]], rownames(fx$data))
  expect_true(all(colnames(fx$data) %in% names(back)))
  expect_true(all(back[["group"]] == 0))
})

test_that("export_metaboanalyst uses a supplied group_col", {
  fx  <- make_export_omiprep()
  dir <- withr::local_tempdir()
  suppressMessages(export_metaboanalyst(fx$m, dir, group_col = "sex"))

  fp <- list.files(dir, pattern = "^omiprep_metaboanalyst_export_.*\\.csv$",
                   full.names = TRUE)
  back <- read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE)
  expect_true("sex" %in% names(back))
  # row order is by PatientID after merge(); compare as a keyed lookup
  expected_sex <- as.character(fx$m@samples[["sex"]])
  names(expected_sex) <- fx$m@samples[["sample_id"]]
  expect_equal(back[["sex"]], unname(expected_sex[back[["PatientID"]]]))
})
