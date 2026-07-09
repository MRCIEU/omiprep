library(testthat)

# Tests for the standalone statistics helpers:
#   total_sum_abundance() (R/total_sum_abundance.R)
#   feature_skewness()    (R/feature_skewness.R)
#   cramerV()             (R/cramer_v.R)
#   variable_by_factor()  (R/variable_by_factor.R)  -- plot, smoke-tested
#   multivariate_anova()  (R/multivariate_anova.R)  -- plot, smoke-tested

# ---- total_sum_abundance ----------------------------------------------------

test_that("total_sum_abundance (no z-transform) sums rows and complete features", {
  M <- rbind(s1 = c(1, 2, 3), s2 = c(4, 5, 6), s3 = c(7, 8, NA))
  colnames(M) <- c("f1", "f2", "f3")
  out <- total_sum_abundance(M, ztransform = FALSE)

  expect_s3_class(out, "data.frame")
  expect_identical(names(out), c("sample_id", "non_na_feature_count", "tsa_total",
                                 "tsa_complete_features", "complete_feature_count"))
  expect_identical(out$sample_id, rownames(M))
  # tsa_total uses all features (NA dropped); tsa_complete uses only f1,f2
  expect_equal(out$tsa_total, c(6, 15, 15), ignore_attr = TRUE)
  expect_equal(out$tsa_complete_features, c(3, 9, 15), ignore_attr = TRUE)
  # f3 is NA for s3; 2 of 3 features (f1,f2) are complete across all samples
  expect_equal(out$non_na_feature_count, c(3, 3, 2), ignore_attr = TRUE)
  expect_equal(out$complete_feature_count, rep(2, 3), ignore_attr = TRUE)
})

test_that("total_sum_abundance equals complete-feature TSA when no data missing", {
  set.seed(11)
  M <- matrix(runif(30, 1, 100), nrow = 6)
  rownames(M) <- paste0("s", 1:6); colnames(M) <- paste0("f", 1:5)
  out <- total_sum_abundance(M, ztransform = FALSE)
  expect_equal(out$tsa_total, out$tsa_complete_features)
  expect_equal(out$tsa_total, rowSums(M), ignore_attr = TRUE)
})

test_that("total_sum_abundance z-transform yields non-negative finite sums", {
  set.seed(12)
  M <- matrix(rnorm(40), nrow = 8)
  rownames(M) <- paste0("s", 1:8); colnames(M) <- paste0("f", 1:5)
  out <- total_sum_abundance(M, ztransform = TRUE)
  # shifting each column by abs(global min) makes all values >= 0
  expect_true(all(out$tsa_total >= 0))
  expect_true(all(is.finite(out$tsa_total)))
  expect_identical(nrow(out), nrow(M))
})

# ---- feature_skewness -------------------------------------------------------

make_skew_matrix <- function() {
  cbind(
    left  = c(rep(10, 18), 1, 2),    # low-value tail  -> negative skew
    right = c(rep(1, 18), 19, 20),   # high-value tail -> positive skew
    sym   = c(1:10, 10:1)            # roughly symmetric
  )
}

test_that("feature_skewness returns per-feature skew with the right signs", {
  X  <- make_skew_matrix()
  fs <- feature_skewness(X)

  expect_s3_class(fs, "data.frame")
  expect_identical(fs$feature_id, colnames(X))
  expect_lt(fs$skew[fs$feature_id == "left"], 0)
  expect_gt(fs$skew[fs$feature_id == "right"], 0)
  # threshold NULL -> no exclusion decision
  expect_true(all(is.na(fs$exclude_by_skewness)))
})

test_that("feature_skewness flags by direction relative to threshold", {
  X <- make_skew_matrix()

  left  <- feature_skewness(X, threshold = 0.5, direction = "left")
  expect_true(left$exclude_by_skewness[left$feature_id == "left"])
  expect_false(left$exclude_by_skewness[left$feature_id == "right"])

  right <- feature_skewness(X, threshold = 0.5, direction = "right")
  expect_true(right$exclude_by_skewness[right$feature_id == "right"])
  expect_false(right$exclude_by_skewness[right$feature_id == "left"])

  both <- feature_skewness(X, threshold = 0.5, direction = "both")
  expect_true(both$exclude_by_skewness[both$feature_id == "left"])
  expect_true(both$exclude_by_skewness[both$feature_id == "right"])
  expect_false(both$exclude_by_skewness[both$feature_id == "sym"])
})

test_that("feature_skewness validates its inputs", {
  X <- make_skew_matrix()
  expect_error(feature_skewness(as.data.frame(X)), "must be a matrix")
  expect_error(feature_skewness(X, direction = "sideways"), "should be one of|arg")
  expect_error(feature_skewness(X, threshold = -1), "non-negative")
})

# ---- cramerV ----------------------------------------------------------------

test_that("cramerV is 1 for perfect association and 0 for independence", {
  expect_equal(unname(cramerV(matrix(c(10, 0, 0, 10), 2, 2))), 1)
  expect_equal(unname(cramerV(matrix(c(5, 5, 5, 5), 2, 2))), 0)
  # vector inputs: perfectly associated nominal vectors
  expect_equal(unname(cramerV(c("a", "a", "b", "b"), c("p", "p", "q", "q"))), 1)
})

test_that("cramerV returns a named scalar without ci and a data.frame with ci", {
  v <- cramerV(matrix(c(10, 2, 3, 12), 2, 2))
  expect_named(v, "Cramer V")
  expect_length(v, 1)

  ci <- cramerV(matrix(c(10, 2, 3, 12), 2, 2), ci = TRUE, R = 200)
  expect_s3_class(ci, "data.frame")
  expect_identical(names(ci), c("Cramer.V", "lower.ci", "upper.ci"))
  expect_identical(nrow(ci), 1L)
})

# ---- plot helpers (smoke / sanity) ------------------------------------------

test_that("variable_by_factor returns a ggplot for violin and box variants", {
  set.seed(20)
  x <- c(rnorm(20, 10, 2), rnorm(20, 20, 2))
  y <- factor(c(rep("A", 20), rep("B", 20)))
  expect_s3_class(variable_by_factor(x, y, violin = TRUE), "ggplot")
  expect_s3_class(variable_by_factor(x, y, violin = FALSE), "ggplot")
})

test_that("multivariate_anova returns a ggplot table figure", {
  set.seed(21)
  n  <- 80
  df <- data.frame(age = rnorm(n, 45, 8), bmi = rnorm(n, 25, 4))
  dep <- 0.3 * df$age - 0.2 * df$bmi + rnorm(n)
  expect_s3_class(multivariate_anova(dep = dep, indep_df = df), "ggplot")
})

# ---- remove_perfect_correlation (internal) ----------------------------------
# Regression: cor() of perfectly collinear data returns one ULP below 1, so the
# previous exact `== 1` test never fired. The tolerance (tol = 0.001) treats
# |r| >= 0.999 as collinear.

test_that("remove_perfect_correlation drops a perfectly collinear column", {
  df  <- data.frame(a = 1:20, b = 2 * (1:20) + 5, c = (1:20)^2 - 3 * (1:20))
  out <- remove_perfect_correlation(df)
  # a and b are perfectly collinear -> exactly one survives; c is retained
  expect_identical(ncol(out), 2L)
  expect_true("c" %in% names(out))
  expect_length(intersect(names(out), c("a", "b")), 1L)
})

test_that("remove_perfect_correlation handles perfect negative collinearity", {
  df  <- data.frame(a = 1:20, neg = -(1:20))
  out <- remove_perfect_correlation(df)
  expect_identical(ncol(out), 1L)
})

test_that("remove_perfect_correlation keeps merely-correlated columns", {
  set.seed(30)
  a <- rnorm(200)
  df <- data.frame(
    a        = a,
    moderate = a + rnorm(200, sd = 0.7),   # ~0.8 correlation, below 0.999
    indep    = rnorm(200)
  )
  out <- remove_perfect_correlation(df)
  expect_identical(names(out), names(df))   # nothing dropped
})
