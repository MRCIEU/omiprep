library(testthat)

# Regression tests for the two-sided symmetric winsorization branch of
# quality_control() (R/quality_control.R, outlier_treatment = "winsorize").
#
# Behaviour being pinned:
#   1. Outliers are identified PER FEATURE (not the global union of all
#      outlying samples across features). The reference quantiles for a
#      feature are computed from that feature's own non-outlier samples.
#   2. Clamping is symmetric: upper-tail outliers -> winsorize_quantile,
#      lower-tail outliers -> (1 - winsorize_quantile).
#   3. The split between upper and lower outliers is by the per-feature
#      median of the reference (non-outlier) values.
#   4. Non-outlier values are left untouched in the "qc" layer.
#   5. winsorize_quantile outside [0.5, 1] errors.
#
# The run is configured to exclude NOTHING (missingness / TSA / PCA outlier
# filters all disabled) so the winsorized "qc" layer is not subsequently
# NA-masked and can be asserted on directly.

# ---- shared fixture ---------------------------------------------------------

make_winsorize_fixture <- function() {
  set.seed(1234)
  n_samples  <- 60
  n_features <- 30
  data <- matrix(rnorm(n_samples * n_features, mean = 10, sd = 1),
                 nrow = n_samples, ncol = n_features)
  rownames(data) <- paste0("Sample_", seq_len(n_samples))
  colnames(data) <- paste0("Met_", seq_len(n_features))

  # Feature 1: extreme outliers in samples 1, 2 (high) and 3 (low).
  data["Sample_1", "Met_1"] <-  30
  data["Sample_2", "Met_1"] <-  30
  data["Sample_3", "Met_1"] <- -10
  # Feature 2: DISJOINT outlier samples 30, 31 (high) and 32 (low).
  data["Sample_30", "Met_2"] <-  30
  data["Sample_31", "Met_2"] <-  30
  data["Sample_32", "Met_2"] <- -10

  # Union-vs-subset guard: give each feature's reference set some legitimately
  # high (but NON-outlier) values located at the OTHER feature's outlier
  # samples. If the code wrongly used the global union of outlier samples to
  # build the reference set, these values would be dropped and the reference
  # quantiles (hence the clamp targets) would shift measurably.
  data[c("Sample_30", "Sample_31", "Sample_32"), "Met_1"] <- 13
  data[c("Sample_1",  "Sample_2",  "Sample_3"),  "Met_2"] <- 13

  samples  <- data.frame(sample_id  = rownames(data))
  features <- data.frame(feature_id = colnames(data))

  list(data = data,
       m = Omiprep(data, samples, features))
}

# Independent re-implementation of the documented winsorize spec for one
# feature, used to compute the expected clamp targets. `restrict_to` lets a
# test simulate the (buggy) "global union" reference set.
expected_winsor <- function(data, fid, outlier_samps,
                            lower_q, upper_q, restrict_to = outlier_samps) {
  ref_vals <- data[!(rownames(data) %in% restrict_to), fid]
  qs       <- stats::quantile(ref_vals, probs = c(lower_q, upper_q), na.rm = TRUE)
  ref_med  <- stats::median(data[!(rownames(data) %in% outlier_samps), fid],
                            na.rm = TRUE)
  raw      <- data[outlier_samps, fid]
  list(
    upper_target = unname(qs[2]),
    lower_target = unname(qs[1]),
    upper_samps  = outlier_samps[raw >= ref_med],
    lower_samps  = outlier_samps[raw <  ref_med]
  )
}

run_qc_winsorize <- function(m, winsorize_quantile = 0.95) {
  quality_control(
    m,
    source_layer           = "input",
    sample_missingness     = NULL,   # disable all sample/feature dropping so the
    feature_missingness    = NULL,   # winsorized qc layer is never NA-masked
    total_sum_abundance_sd = NULL,
    pc_outlier_sd          = NULL,
    outlier_udist          = 5,
    outlier_treatment      = "winsorize",
    winsorize_quantile     = winsorize_quantile,
    max_num_pcs            = 2,
    cores                  = 1
  )
}

# ---- tests ------------------------------------------------------------------

test_that("winsorize clamps per-feature outliers to symmetric reference quantiles", {
  fx   <- make_winsorize_fixture()
  data <- fx$data
  m    <- suppressMessages(run_qc_winsorize(fx$m, winsorize_quantile = 0.95))

  qc <- m@data[, , "qc"]

  # outlier_detection should flag exactly the injected samples per feature
  expect_setequal(
    rownames(data)[outlier_detection(data, nsd = 5, by = "column")[, "Met_1"] == 1],
    c("Sample_1", "Sample_2", "Sample_3")
  )
  expect_setequal(
    rownames(data)[outlier_detection(data, nsd = 5, by = "column")[, "Met_2"] == 1],
    c("Sample_30", "Sample_31", "Sample_32")
  )

  for (fid in c("Met_1", "Met_2")) {
    outs <- if (fid == "Met_1") c("Sample_1", "Sample_2", "Sample_3")
            else                c("Sample_30", "Sample_31", "Sample_32")
    exp  <- expected_winsor(data, fid, outs, lower_q = 0.05, upper_q = 0.95)

    # unname(): matrix subsetting drops names when the index is length-1, so
    # compare values only (cell identity is fixed by the index itself).
    expect_equal(unname(qc[exp$upper_samps, fid]),
                 rep(exp$upper_target, length(exp$upper_samps)),
                 info = paste("upper clamp,", fid))
    expect_equal(unname(qc[exp$lower_samps, fid]),
                 rep(exp$lower_target, length(exp$lower_samps)),
                 info = paste("lower clamp,", fid))
  }
})

test_that("winsorize reference set is per-feature, not the global outlier union", {
  fx   <- make_winsorize_fixture()
  data <- fx$data
  m    <- suppressMessages(run_qc_winsorize(fx$m, winsorize_quantile = 0.95))
  qc   <- m@data[, , "qc"]

  outs        <- c("Sample_1", "Sample_2", "Sample_3")
  union_samps <- c(outs, "Sample_30", "Sample_31", "Sample_32")

  correct <- expected_winsor(data, "Met_1", outs,
                             lower_q = 0.05, upper_q = 0.95)
  buggy   <- expected_winsor(data, "Met_1", outs,
                             lower_q = 0.05, upper_q = 0.95,
                             restrict_to = union_samps)

  # the fixture is built so the two differ measurably
  expect_gt(abs(correct$upper_target - buggy$upper_target), 0.1)

  # actual matches the per-feature (correct) target, not the union (buggy) one
  actual_upper <- qc["Sample_1", "Met_1"]
  expect_equal(actual_upper, correct$upper_target)
  expect_false(isTRUE(all.equal(actual_upper, buggy$upper_target)))
})

test_that("winsorize leaves non-outlier values untouched", {
  fx   <- make_winsorize_fixture()
  data <- fx$data
  m    <- suppressMessages(run_qc_winsorize(fx$m, winsorize_quantile = 0.95))
  qc   <- m@data[, , "qc"]

  outlier_cells <- list(
    Met_1 = c("Sample_1", "Sample_2", "Sample_3"),
    Met_2 = c("Sample_30", "Sample_31", "Sample_32")
  )
  for (fid in colnames(data)) {
    untouched <- setdiff(rownames(data), outlier_cells[[fid]])
    expect_equal(qc[untouched, fid], data[untouched, fid],
                 info = paste("non-outlier values preserved,", fid))
  }
})

test_that("winsorize_quantile outside [0.5, 1] errors", {
  fx <- make_winsorize_fixture()
  expect_error(
    suppressMessages(run_qc_winsorize(fx$m, winsorize_quantile = 0.3)),
    "winsorize_quantile"
  )
})
