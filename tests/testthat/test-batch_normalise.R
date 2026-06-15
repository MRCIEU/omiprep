library(testthat)

# Tests for batch_normalise() (R/batch_normalise.R). It appends a dest_layer
# (initialised as a copy of the "input" layer) and divides each value by its
# per-batch, per-feature median. Critically, different run modes can map to
# DIFFERENT batch columns in the sample data (e.g. NEG features batched by
# neg_batch, POS features by pos_batch).

# ---- fixture ----------------------------------------------------------------

make_bn_fixture <- function() {
  X <- rbind(
    S1 = c(2,  1, 3, 2),
    S2 = c(4,  1, 3, 4),
    S3 = c(6,  4, 9, 6),
    S4 = c(10, 5, 9, 8),
    S5 = c(20, 5, 7, NA),   # NA cell to check missing-value handling
    S6 = c(30, 5, 7, 200)
  )
  colnames(X) <- c("F1", "F2", "F3", "F4")
  samples <- data.frame(
    sample_id = rownames(X),
    neg_batch = c("A", "A", "A", "B", "B", "B"),  # NEG features use this
    pos_batch = c("P", "P", "P", "P", "Q", "Q")   # POS features use this
  )
  features <- data.frame(
    feature_id = colnames(X),
    run_mode   = c("NEG", "NEG", "POS", "POS")
  )
  list(
    X = X,
    colmap = c(NEG = "neg_batch", POS = "pos_batch"),
    m = suppressMessages(Omiprep(X, samples, features)),
    samples = samples, features = features
  )
}

# Independent re-implementation of the documented normalisation, used to check
# the values produced by batch_normalise().
expected_bn <- function(X, samples, features, colmap, run_mode_col = "run_mode") {
  out <- X
  for (mode in names(colmap)) {
    feats <- features$feature_id[features[[run_mode_col]] == mode]
    bcol  <- colmap[[mode]]
    for (bid in unique(samples[[bcol]])) {
      s <- samples$sample_id[samples[[bcol]] == bid]
      for (f in feats) {
        med <- stats::median(X[s, f], na.rm = TRUE)
        out[s, f] <- X[s, f] / med
      }
    }
  }
  out
}

run_bn <- function(fx, ...) {
  suppressMessages(batch_normalise(
    fx$m, run_mode_col = "run_mode", run_mode_colmap = fx$colmap, ...))
}

# ---- structure --------------------------------------------------------------

test_that("batch_normalise appends a dest layer and leaves input untouched", {
  fx <- make_bn_fixture()
  o  <- run_bn(fx)

  expect_identical(dim(o@data)[3], dim(fx$m@data)[3] + 1L)
  expect_true("batch_normalised" %in% dimnames(o@data)[[3]])
  # source data preserved exactly
  expect_identical(o@data[, , "input"], fx$X)
  # metadata unchanged
  expect_identical(o@samples, fx$m@samples)
  expect_identical(o@features, fx$m@features)
})

test_that("dest_layer name is configurable", {
  fx <- make_bn_fixture()
  o  <- run_bn(fx, dest_layer = "bn_custom")
  expect_true("bn_custom" %in% dimnames(o@data)[[3]])
  expect_false("batch_normalised" %in% dimnames(o@data)[[3]])
})

# ---- correctness ------------------------------------------------------------

test_that("normalised values equal input divided by per-batch per-feature median", {
  fx <- make_bn_fixture()
  o  <- run_bn(fx)
  expect_equal(
    o@data[, , "batch_normalised"],
    expected_bn(fx$X, fx$samples, fx$features, fx$colmap)
  )
})

test_that("each mode is batched by its own mapped sample column", {
  fx <- make_bn_fixture()
  o  <- run_bn(fx)
  bn <- o@data[, , "batch_normalised"]

  # F3 is POS -> batched by pos_batch (P = S1:S4, Q = S5:S6), NOT neg_batch.
  # P median of c(3,3,9,9) = 6 -> 0.5,0.5,1.5,1.5 ; Q median c(7,7)=7 -> 1,1
  expect_equal(unname(bn[c("S1", "S2", "S3", "S4"), "F3"]),
               c(0.5, 0.5, 1.5, 1.5))
  expect_equal(unname(bn[c("S5", "S6"), "F3"]), c(1, 1))
})

test_that("post-normalisation per-batch per-feature median is 1", {
  fx <- make_bn_fixture()
  o  <- run_bn(fx)
  bn <- o@data[, , "batch_normalised"]

  for (mode in names(fx$colmap)) {
    feats <- fx$features$feature_id[fx$features$run_mode == mode]
    bcol  <- fx$colmap[[mode]]
    for (bid in unique(fx$samples[[bcol]])) {
      s <- fx$samples$sample_id[fx$samples[[bcol]] == bid]
      for (f in feats) {
        expect_equal(stats::median(bn[s, f], na.rm = TRUE), 1,
                     info = sprintf("mode=%s batch=%s feature=%s", mode, bid, f))
      }
    }
  }
})

test_that("NA values are preserved through normalisation", {
  fx <- make_bn_fixture()
  o  <- run_bn(fx)
  expect_true(is.na(o@data["S5", "F4", "batch_normalised"]))
  # and nowhere else NA was introduced
  expect_identical(which(is.na(o@data[, , "batch_normalised"])),
                   which(is.na(fx$X)))
})

# ---- validation -------------------------------------------------------------

test_that("batch_normalise rejects a run_mode_col not in the features data", {
  fx <- make_bn_fixture()
  expect_error(
    suppressMessages(batch_normalise(fx$m, run_mode_col = "not_a_column",
                                     run_mode_colmap = fx$colmap)),
    "should be one of|arg"
  )
})

test_that("batch_normalise requires colmap to cover every run mode", {
  fx <- make_bn_fixture()
  expect_error(
    suppressMessages(batch_normalise(fx$m, run_mode_col = "run_mode",
                                     run_mode_colmap = c(NEG = "neg_batch"))),
    "names must match"
  )
})

test_that("batch_normalise requires colmap values to be sample columns", {
  fx <- make_bn_fixture()
  expect_error(
    suppressMessages(batch_normalise(
      fx$m, run_mode_col = "run_mode",
      run_mode_colmap = c(NEG = "neg_batch", POS = "not_a_sample_col"))),
    "must be columns in sample data"
  )
})

test_that("batch_normalise rejects an unknown source_layer", {
  fx <- make_bn_fixture()
  expect_error(
    suppressMessages(batch_normalise(
      fx$m, run_mode_col = "run_mode", run_mode_colmap = fx$colmap,
      source_layer = "not_a_layer")),
    "should be one of|arg"
  )
})

# Regression for the previously hardcoded "input" seed: the dest layer must be
# seeded from and normalised against `source_layer`, not always "input".
test_that("source_layer (not 'input') seeds and drives the normalisation", {
  fx <- make_bn_fixture()
  # add a 'qc' layer that clearly differs from input
  qc <- fx$X + 100
  m2 <- fx$m
  m2@data <- add_layer(m2@data, qc, "qc")

  o <- suppressMessages(batch_normalise(
    m2, run_mode_col = "run_mode", run_mode_colmap = fx$colmap,
    source_layer = "qc"))

  # dest == qc divided by the per-batch per-feature median OF QC
  expect_equal(o@data[, , "batch_normalised"],
               expected_bn(qc, fx$samples, fx$features, fx$colmap))
  # both source layers preserved
  expect_identical(o@data[, , "input"], fx$X)
  expect_identical(o@data[, , "qc"], qc)

  # and it must NOT match the old behaviour: input seed / median(qc)
  old_buggy <- fx$X
  for (mode in names(fx$colmap)) {
    feats <- fx$features$feature_id[fx$features$run_mode == mode]
    bcol  <- fx$colmap[[mode]]
    for (bid in unique(fx$samples[[bcol]])) {
      s <- fx$samples$sample_id[fx$samples[[bcol]] == bid]
      for (f in feats) old_buggy[s, f] <- fx$X[s, f] / stats::median(qc[s, f], na.rm = TRUE)
    }
  }
  expect_false(isTRUE(all.equal(o@data[, , "batch_normalised"], old_buggy)))
})
