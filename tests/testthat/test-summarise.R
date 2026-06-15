library(testthat)

# Tests for summarise() (R/summarise.R), the QC backbone that orchestrates
# feature_summary(), sample_summary(), pc_and_outliers() and
# tree_and_independent_features(), then either returns a list of data.frames
# or writes the feature/sample summary array layers back into the object.
#
# The fixture is sized (60 x 24) so pc_and_outliers() takes the deterministic
# prcomp() branch rather than the stochastic irlba() one. Correlated feature
# blocks are injected so the tree-cut independence step has clusters to
# collapse, plus one outlier cell and one NA for missingness/outlier checks.

make_summarise_fixture <- function() {
  withr::with_seed(7, {
    ns <- 60L; nf <- 24L
    X <- matrix(rnorm(ns * nf), ns, nf)
    # two correlated feature blocks -> non-trivial independence structure
    X[, 2] <- X[, 1] + rnorm(ns, sd = 0.05)
    X[, 3] <- X[, 1] + rnorm(ns, sd = 0.05)
    X[, 5] <- X[, 4] + rnorm(ns, sd = 0.05)
    X[6, 7] <- 50      # an extreme outlier cell (sample S6, feature F7)
    X[2, 9] <- NA       # a single missing cell (feature F9)
    rownames(X) <- paste0("S", seq_len(ns))
    colnames(X) <- paste0("F", seq_len(nf))
  })
  m <- suppressMessages(
    Omiprep(X, data.frame(sample_id = rownames(X)),
            data.frame(feature_id = colnames(X)))
  )
  list(X = X, m = m)
}

run_summarise <- function(m, ...) {
  suppressMessages(summarise(m, cores = 1, ...))
}

# ---- data.frame mode: top-level shape --------------------------------------

test_that("summarise(data.frame) returns sample and feature summary data.frames", {
  fx  <- make_summarise_fixture()
  res <- run_summarise(fx$m, output = "data.frame")

  expect_type(res, "list")
  expect_named(res, c("sample_summary", "feature_summary"))
  expect_s3_class(res$sample_summary, "data.frame")
  expect_s3_class(res$feature_summary, "data.frame")
})

# ---- feature summary --------------------------------------------------------

test_that("feature_summary has one row per feature, in data order", {
  fx  <- make_summarise_fixture()
  res <- run_summarise(fx$m, output = "data.frame")
  fs  <- res$feature_summary

  expect_identical(nrow(fs), ncol(fx$X))
  expect_identical(fs$feature_id, colnames(fx$X))
  expect_true(all(c("missingness", "mean", "sd", "outlier_count",
                    "k", "independent_features") %in% names(fs)))
})

test_that("feature_summary captures missingness and per-feature outliers", {
  fx  <- make_summarise_fixture()
  res <- run_summarise(fx$m, output = "data.frame")
  fs  <- res$feature_summary

  # F9 has exactly one missing value of 60; F1 is fully present
  expect_equal(fs["F9", "missingness"], 1 / 60)
  expect_equal(fs["F1", "missingness"], 0)
  # the injected extreme cell makes F7 carry at least one outlier
  expect_gte(fs["F7", "outlier_count"], 1)
})

test_that("feature_summary independence collapses correlated blocks", {
  fx  <- make_summarise_fixture()
  res <- run_summarise(fx$m, output = "data.frame")
  fs  <- res$feature_summary

  expect_type(fs$independent_features, "logical")
  # not every feature is independent once correlated blocks collapse
  expect_lt(sum(fs$independent_features, na.rm = TRUE), ncol(fx$X))
  # each injected correlated block yields exactly one representative
  expect_equal(sum(fs[c("F1", "F2", "F3"), "independent_features"]), 1L)
  expect_equal(sum(fs[c("F4", "F5"), "independent_features"]), 1L)
})

# ---- sample summary ---------------------------------------------------------

test_that("sample_summary has one row per sample with PCA outlier columns", {
  fx  <- make_summarise_fixture()
  res <- run_summarise(fx$m, output = "data.frame")
  ss  <- res$sample_summary

  expect_identical(nrow(ss), nrow(fx$X))
  expect_identical(ss$sample_id, rownames(fx$X))
  expect_true(all(c("missingness", "tsa_total", "outlier_count", "pc1")
                  %in% names(ss)))
  # PCA outlier flag columns are present for the 3/4/5 SD thresholds
  expect_true(any(grepl("_3_sd_outlier$", names(ss))))
  expect_true(any(grepl("_4_sd_outlier$", names(ss))))
  expect_true(any(grepl("_5_sd_outlier$", names(ss))))
  # the sample carrying the extreme cell shows a feature-level outlier
  expect_gte(ss["S6", "outlier_count"], 1)
})

# ---- object mode ------------------------------------------------------------

test_that("summarise(object) writes feature and sample summary array layers", {
  fx <- make_summarise_fixture()
  o  <- run_summarise(fx$m, output = "object")

  expect_s3_class(o, "omiprep::Omiprep")

  # feature_summary: stats x features x layers
  expect_length(dim(o@feature_summary), 3)
  expect_identical(dim(o@feature_summary)[2], ncol(fx$X))
  expect_identical(dimnames(o@feature_summary)[[2]], colnames(fx$X))
  expect_true("input" %in% dimnames(o@feature_summary)[[3]])

  # sample_summary: samples x stats x layers
  expect_length(dim(o@sample_summary), 3)
  expect_identical(dim(o@sample_summary)[1], nrow(fx$X))
  expect_identical(dimnames(o@sample_summary)[[1]], rownames(fx$X))
  expect_true("input" %in% dimnames(o@sample_summary)[[3]])
})

test_that("object and data.frame modes produce consistent summary values", {
  fx  <- make_summarise_fixture()
  res <- run_summarise(fx$m, output = "data.frame")
  o   <- run_summarise(fx$m, output = "object")

  # object feature layer == transpose of the data.frame (minus the id column)
  fs_mat <- t(as.matrix(res$feature_summary[, setdiff(names(res$feature_summary),
                                                      "feature_id")]))
  expect_equal(o@feature_summary[, , "input"], fs_mat, ignore_attr = TRUE)

  ss_mat <- as.matrix(res$sample_summary[, setdiff(names(res$sample_summary),
                                                   "sample_id")])
  expect_equal(o@sample_summary[, , "input"], ss_mat, ignore_attr = TRUE)
})

# ---- processing attributes --------------------------------------------------

test_that("summarise stores tree, outlier and PCA processing attributes", {
  fx  <- make_summarise_fixture()
  res <- run_summarise(fx$m, output = "data.frame")

  expect_s3_class(attr(res$feature_summary, "input_tree"), "hclust")
  expect_equal(attr(res$feature_summary, "input_outlier_udist"), 5)
  expect_equal(attr(res$feature_summary, "input_tree_cut_height"), 0.5)

  expect_type(attr(res$sample_summary, "input_varexp"), "double")
  expect_gt(length(attr(res$sample_summary, "input_varexp")), 0)
  expect_gte(attr(res$sample_summary, "input_num_pcs_scree"), 2)
  expect_equal(attr(res$sample_summary, "input_outlier_udist"), 5)
})

test_that("outlier_udist is threaded through to the summaries", {
  fx  <- make_summarise_fixture()
  res <- run_summarise(fx$m, output = "data.frame", outlier_udist = 3)
  expect_equal(attr(res$feature_summary, "input_outlier_udist"), 3)
  expect_equal(attr(res$sample_summary,  "input_outlier_udist"), 3)
})

# ---- input validation -------------------------------------------------------

test_that("summarise validates its arguments", {
  fx <- make_summarise_fixture()
  expect_error(run_summarise(fx$m, output = "nope"), "should be one of|arg")
  expect_error(run_summarise(fx$m, source_layer = "bogus"), "should be one of|arg")
  expect_error(run_summarise(fx$m, sample_ids = "ghost_sample"),
               "must all be found")
  expect_error(run_summarise(fx$m, feature_ids = "ghost_feature"),
               "must all be found")
})
