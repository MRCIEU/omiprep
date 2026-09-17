#' @title Estimate Feature Skewness
#' @description
#' Estimate per-feature skewness and optionally flag features beyond a user-defined
#' skewness threshold.
#'
#' @details
#' Skewness is calculated for each feature across samples using [psych::describe()],
#' which by default returns the type 3 estimator
#' \eqn{b_1 = g_1 \cdot ((n-1)/n)^{3/2}}{b1 = g1 * ((n-1)/n)^(3/2)}, where \eqn{g_1} is the standard
#' moment coefficient of skewness and \eqn{n} is the number of non-missing values for
#' that feature.
#'
#' Missing values are dropped feature-wise, so \eqn{n} varies between features and no
#' imputation is performed. No minimum sample-size is imposed: a feature with only a
#' handful of non-missing values will still return an estimate, and such estimates are
#' unstable. Features with zero variance return `NaN` and are never flagged.
#'
#' The threshold is a cutoff on the skewness statistic itself, not a p-value, and is
#' applied to the **signed** value unless `direction = "both"`:
#'
#' \itemize{
#'   \item `direction = "left"` flags features where `skew <= -threshold`;
#'   \item `direction = "right"` flags features where `skew >= threshold`;
#'   \item `direction = "both"` flags features where `abs(skew) >= threshold`.
#' }
#'
#' Skewness is computed on the values as supplied. It is not invariant to
#' transformation, so a feature's skewness on a raw intensity scale and on a
#' log scale can differ substantially.
#'
#' A skewed distribution is a description of shape, not evidence of technical
#' failure. Genuine biological variation, subgroup structure and low-frequency
#' phenomena all produce skewed features. Treat the flag as a prompt to inspect a
#' feature rather than as grounds to remove it; see `vignette("skewness_qc")` for
#' guidance on interpretation.
#'
#' @param data matrix, a numeric matrix with samples in rows and features in columns.
#' @param threshold numeric, optional skewness threshold. If `NULL` (the default),
#' only skewness is returned and no exclusion flag is calculated.
#' @param direction character, direction of skewness to flag. One of `"left"`
#' (the default), `"right"`, or `"both"`.
#'
#' @return data.frame with columns `feature_id`, `skew`, and
#' `exclude_by_skewness` (logical; `NA` if `threshold = NULL`).
#'
#' @seealso [quality_control()] for applying this filter within the QC pipeline,
#' and [feature_describe()], which also reports per-feature skewness.
#'
#' @importFrom psych describe
#' @export
feature_skewness <- function(data, threshold = NULL, direction = "left") {
  stopifnot("data must be a matrix" = is.matrix(data))
  direction <- match.arg(direction, choices = c("left", "right", "both"))
  stopifnot("threshold must be NULL or a non-negative scalar" =
              is.null(threshold) || (is.numeric(threshold) && length(threshold) == 1 && !is.na(threshold) && threshold >= 0))

  desc <- psych::describe(data)
  out <- data.frame(
    feature_id = rownames(desc),
    skew = as.numeric(desc[, "skew"])
  )

  if (is.null(threshold)) {
    out$exclude_by_skewness <- NA
  } else {
    out$exclude_by_skewness <- switch(
      direction,
      "left" = out$skew <= -threshold,
      "right" = out$skew >= threshold,
      "both" = abs(out$skew) >= threshold
    )
  }

  out
}
