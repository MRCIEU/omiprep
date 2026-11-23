#' @title Simulate an omics dataset with known problems
#' @description
#' Generates a synthetic omics data matrix (samples-rows, features-cols), with correlated
#' feature clusters, missing values and planted quality problems (optional via the arguments).
#'
#' @details
#' 1. **Base data.** Randomly fill the matrix with independent random numbers. The first
#'    `frac_clustered` of the features are then givne one of `n_factors` latent factors and
#'    mixed with it (loading drawn uniformly from `loading_range`) to give clusters of
#'    correlated features. Each feature is then given its own mean (between 8 - 14) and
#'    own standard deviation (between 0.4 - 1.5).
#' 2. **Planted problems and biological controls** (`plant` argument), applied to the
#'    complete data. Each sample or feature gets one, or none (if `plant = NULL`), of them.
#' 3. **Missing values**, injected according to `miss_mech` and `miss_rate`.
#' 4. **Extra missing values** for planted poor features and failed samples.
#'
#' @section Planted problems (`plant`):
#' `plant` is a named list. Each `*_frac` entry is the share (0 to 1) of samples or features given
#' that problem (default 0). Entries not supplied keep their default.
#'
#' Technical problems:
#'
#' * `dilution_frac`, `dilution_shift` (default 1): diluted or concentrated samples. Every feature
#'   of the sample is shifted up or down (direction chosen at random per sample) by
#'   `dilution_shift` feature standard deviations. Label = `"tech_dilution"`.
#' * `outlier_frac`, `outlier_shift` (default 1): outlier samples. Each sample is shifted along one
#'   random direction across the features, shared by the planted samples and not present in the
#'   rest of the data, scaled by `outlier_shift` feature standard deviations. Label =
#'   `"tech_outlier"`.
#' * `failed_frac`, `failed_miss` (default `NULL`): failed samples, given extra values missing at
#'   random. `failed_miss` is the share of features set missing; if `NULL` it is drawn per sample
#'   from 0.5 to 0.9. Label = `"tech_failed"`.
#' * `poor_feat_frac`, `poor_miss` (default `NULL`), `poor_mech` (default `"random"`): poorly
#'   detected features, given extra missing values. `poor_miss` is the share of samples set
#'   missing; if `NULL` it is drawn per feature from 0.5 to 0.95. With `poor_mech = "random"` the
#'   extra values are missing at random; with `"lod"` the lowest values still observed are set
#'   missing (simulating the scenario where values are below a detection limit). Label =
#'   `"tech_poor_feature"`.
#' * `trunc_frac`, `trunc_quantile` (default 0.6): features capped at an upper limit. Values above
#'   the `trunc_quantile` quantile of the feature are set to that quantile, which gives a
#'   left-skewed feature. Label = `"tech_truncated"`.
#'
#' Biological controls (ideally should not be flagged by QC):
#'
#' * `bio_subgroup_frac`, `bio_subgroup_shift` (default 2): a subgroup of samples shifted by
#'   `bio_subgroup_shift` feature standard deviations on the features of the first latent factor
#'   only. Requires `frac_clustered > 0`. Label = `"bio_subgroup"`.
#' * `bio_tsa_frac`, `bio_tsa_shift` (default 0.5): a subgroup of samples with higher abundance,
#'   shifted by `bio_tsa_shift` feature standard deviations on a random 40% of the features.
#'   Label = `"bio_tsa"`.
#' * `rare_frac`: features with a rare effect. A random 8% of samples are shifted up by 3 feature
#'   standard deviations, which gives a right-skewed feature. Label = `"bio_rare_effect"`.
#'
#' @param n integer, the number of samples (rows).
#' @param p integer, the number of features (columns).
#' @param n_factors integer, the number of latent factors behind the correlated feature clusters.
#'   `NULL` (default) uses one factor per five clustered features.
#' @param frac_clustered numeric (0 to 1), the share of features that load on a latent factor.
#'   The remaining features are independent.
#' @param loading_range numeric vector of length 2, the range from which each clustered feature's
#'   loading on its factor is drawn. Values between 0 and 1.
#' @param miss_rate numeric (0 to <1), the target share of missing values, as the mean of the
#'   per-feature rates.
#' @param miss_mech character, how values go missing: `"lod"` (default), the lowest values of each
#'   feature; `"mcar"`, values missing completely at random; `"mixed"`, half of each feature's missing values by each mechanism;
#'   `"none"`, no missing values.
#' @param complete_frac numeric (0 to 1), the share of features that are left with no missing
#'   values from `miss_mech`.
#' @param plant named list of problems to plant, or `NULL` (default) for none. See section
#'   "Planted problems".
#' @param seed a single integer, or `NULL` (default). If supplied, [set.seed()] is called with it
#'   so the same call returns the same data. This changes the random number state of the session.
#' @param equal_rates logical. `FALSE` (default): per-feature missingness rates vary around
#'   `miss_rate` (beta-distributed, at most 0.75). Some features labelled
#'   `"clean"` will then exceed a missingness threshold by themselves (about 1 in 8 features is at
#'   least 20% missing when `miss_rate = 0.1`). `TRUE`: every feature has exactly `miss_rate`, so
#'   only planted features have high missingness.
#' @param miss_ignores_shift logical. `FALSE` (default): values below the detection limit are
#'   chosen after the dilution and outlier shifts, so a sample shifted down also gains missing
#'   values. `TRUE`: they are chosen on the unshifted values, so diluted and outlier samples have
#'   the same missingness as other samples and differ only in their values.
#'
#' @param as_omiprep logical. `FALSE` (default) returns a list; `TRUE` returns an [Omiprep()]
#'   object built from it.
#'
#' @return If `as_omiprep = FALSE`, a list with three elements that can be passed directly to
#' [Omiprep()]:
#' * `data`: numeric matrix, samples in rows and features in columns, with row and column names.
#' * `samples`: data frame with columns `sample_id`, `batch` (a random label, `"b1"` or `"b2"`,
#'   unrelated to the data) and `truth`, the label of what was planted in the sample (`"clean"`
#'   if nothing).
#' * `features`: data frame with columns `feature_id`, `platform` (constant), `truth`, the label
#'   of what was planted in the feature (`"clean"` if nothing), and `cluster_factor`, the latent
#'   factor a clustered feature loads on (`NA` for independent features).
#'
#' If `as_omiprep = TRUE`, an `Omiprep` object holding the same data, with the tables above in
#' its `samples` and `features` properties.
#'
#' @examples
#' # 200 samples by 200 features, 10% of each feature missing, with a few failed samples and
#' # poorly detected features planted
#' sim <- simulate_omics(
#'   n = 200, p = 200, miss_rate = 0.10, equal_rates = TRUE, seed = 1,
#'   plant = list(failed_frac = 0.02, poor_feat_frac = 0.05)
#' )
#' dim(sim$data)
#' table(sim$samples$truth)
#' table(sim$features$truth)
#'
#' # the same data as an Omiprep object; the truth labels are kept in the sample and feature tables
#' m <- simulate_omics(
#'   n = 200, p = 200, miss_rate = 0.10, equal_rates = TRUE, seed = 1, as_omiprep = TRUE,
#'   plant = list(failed_frac = 0.02, poor_feat_frac = 0.05)
#' )
#' table(m@samples$truth)
#'
#' # run quality control and compare what was excluded with what was planted
#' \dontrun{
#' m <- quality_control(m, source_layer = "input", cores = 1)
#' excluded_samples  <- unique(unlist(m@exclusions$samples))
#' excluded_features <- unique(unlist(m@exclusions$features))
#' table(truth = m@samples$truth,  excluded = m@samples$sample_id   %in% excluded_samples)
#' table(truth = m@features$truth, excluded = m@features$feature_id %in% excluded_features)
#' }
#'
#' @md
#' @importFrom stats rnorm runif rbeta quantile
#' @importFrom utils modifyList
#' @export
#'
simulate_omics <- function(n,
                           p,
                           n_factors          = NULL,
                           frac_clustered     = 0.3,
                           loading_range      = c(0.6, 0.9),
                           miss_rate          = 0.10,
                           miss_mech          = c("lod", "mcar", "mixed", "none"),
                           complete_frac      = 0,
                           plant              = NULL,
                           seed               = NULL,
                           equal_rates        = FALSE,
                           miss_ignores_shift = FALSE,
                           as_omiprep         = FALSE) {

  # checks ----
  miss_mech <- match.arg(miss_mech)
  is_count  <- function(x) is.numeric(x) && length(x) == 1 && !is.na(x) && x >= 1 && x == round(x)
  is_frac   <- function(x) is.numeric(x) && length(x) == 1 && !is.na(x) && x >= 0 && x <= 1
  if (!is_count(n) || n < 2) stop("`n` must be a whole number of at least 2, got: ", deparse(n))
  if (!is_count(p) || p < 2) stop("`p` must be a whole number of at least 2, got: ", deparse(p))
  if (!is.null(n_factors) && !is_count(n_factors)) stop("`n_factors` must be NULL or a positive whole number")
  if (!is_frac(frac_clustered)) stop("`frac_clustered` must be a single number between 0 and 1")
  if (!is_frac(complete_frac))  stop("`complete_frac` must be a single number between 0 and 1")
  if (!(is.numeric(loading_range) && length(loading_range) == 2 && all(loading_range >= 0 & loading_range <= 1) && loading_range[1] <= loading_range[2])) {
    stop("`loading_range` must be two increasing numbers between 0 and 1")
  }
  if (!(is_frac(miss_rate) && miss_rate < 1)) stop("`miss_rate` must be a single number from 0 to below 1")
  if (!is.null(seed) && !(is.numeric(seed) && length(seed) == 1 && !is.na(seed))) stop("`seed` must be NULL or a single number")
  stopifnot(is.logical(equal_rates), length(equal_rates) == 1, is.logical(miss_ignores_shift), length(miss_ignores_shift) == 1, is.logical(as_omiprep), length(as_omiprep) == 1)

  defaults <- .simulate_omics_plant_defaults()
  if (!is.null(plant)) {
    if (!is.list(plant) || (length(plant) > 0 && (is.null(names(plant)) || any(names(plant) == "")))) stop("`plant` must be a named list")
    unknown <- setdiff(names(plant), names(defaults))
    if (length(unknown)) {
      stop("Unknown `plant` entries: ", paste(unknown, collapse = ", "), ". Valid entries: ", paste(names(defaults), collapse = ", "))
    }
  }
  pl <- modifyList(defaults, if (is.null(plant)) list() else plant)
  frac_names <- grep("_frac$", names(pl), value = TRUE)
  if (!all(vapply(pl[frac_names], is_frac, logical(1)))) stop("Every `*_frac` entry of `plant` must be a single number between 0 and 1")
  s_frac <- c("dilution_frac", "outlier_frac", "failed_frac", "bio_subgroup_frac", "bio_tsa_frac")
  f_frac <- c("trunc_frac", "rare_frac", "poor_feat_frac")
  if (sum(round(unlist(pl[s_frac]) * n)) > n) stop("The sample fractions in `plant` add up to more than all samples")
  if (sum(round(unlist(pl[f_frac]) * p)) > p) stop("The feature fractions in `plant` add up to more than all features")
  if (!pl$poor_mech %in% c("random", "lod")) stop("`plant$poor_mech` must be \"random\" or \"lod\"")

  if (!is.null(seed)) set.seed(seed)
  sid <- sprintf("s%05d", seq_len(n))
  fid <- sprintf("f%05d", seq_len(p))


  # create base data ----
  n_clustered <- floor(p * frac_clustered)
  if (is.null(n_factors)) {
    n_factors <- max(1, round(n_clustered / 5))
  }
  latent  <- matrix(rnorm(n * n_factors), n, n_factors)
  z       <- matrix(rnorm(n * p), n, p)
  fac_of  <- rep(NA_integer_, p)
  if (n_clustered > 0) {
    fac_of[seq_len(n_clustered)] <- sample(rep_len(seq_len(n_factors), n_clustered))
    load <- runif(n_clustered, loading_range[1], loading_range[2])
    for (j in seq_len(n_clustered)) {
      z[, j] <- load[j] * latent[, fac_of[j]] + sqrt(1 - load[j]^2) * z[, j]
    }
  }


  # put on log-abundance scale (feature-specific mean and SD)
  mu  <- runif(p, 8, 14)
  sdv <- runif(p, 0.4, 1.5)
  mat <- sweep(sweep(z, 2, sdv, `*`), 2, mu, `+`)
  dimnames(mat) <- list(sid, fid)
  ts <- data.frame(sample_id = sid, truth = "clean", stringsAsFactors = FALSE)
  tf <- data.frame(feature_id = fid, truth = "clean", factor = fac_of, stringsAsFactors = FALSE)


  # planted problems (add to complete data, before missingness injected) ---

  # helpers
  pick <- function(pool, k) if (k > 0) sample(pool, k) else integer()
  free_s <- seq_len(n)
  take_s <- function(frac) { k <- round(frac * n); i <- pick(free_s, k); free_s <<- setdiff(free_s, i); i }
  free_f <- seq_len(p)
  take_f <- function(frac, pool = free_f) { k <- round(frac * p); i <- pick(pool, k); free_f <<- setdiff(free_f, i); i }
  shift_rows <- list()  # list of row indexs -> shift added by dilution / outlier planting


  # technical: dilution / concentration (global shift, in feature-SD units)
  i_dil <- take_s(pl$dilution_frac)
  if (length(i_dil)) {
    sgn <- sample(c(-1, 1), length(i_dil), replace = TRUE)
    d <- outer(sgn * pl$dilution_shift, sdv)
    mat[i_dil, ] <- mat[i_dil, ] + d
    for (r in seq_along(i_dil)) shift_rows[[as.character(i_dil[r])]] <- d[r, ]
    ts$truth[i_dil] <- "tech_dilution"
  }


  # technical: outlier samples, shifted along a new direction not shared with the rest of the data
  i_out <- take_s(pl$outlier_frac)
  if (length(i_out)) {
    w <- rnorm(p)
    for (i in i_out) {
      d <- pl$outlier_shift * w * sdv * runif(1, 0.8, 1.2)
      mat[i, ] <- mat[i, ] + d
      shift_rows[[as.character(i)]] <- d
    }
    ts$truth[i_out] <- "tech_outlier"
  }


  # technical: failed samples (extra missingness is applied after the missingness step)
  i_fail <- take_s(pl$failed_frac)
  ts$truth[i_fail] <- "tech_failed"


  # biological: subgroup shifted on the features of the first latent factor
  i_sub <- take_s(pl$bio_subgroup_frac)
  if (length(i_sub) && n_clustered > 0) {
    f1 <- which(fac_of == 1)
    mat[i_sub, f1] <- mat[i_sub, f1] + pl$bio_subgroup_shift * matrix(sdv[f1], length(i_sub), length(f1), byrow = TRUE)
    ts$truth[i_sub] <- "bio_subgroup"
  }


  # biological: higher abundance in a subgroup (40% of features)
  i_tsa <- take_s(pl$bio_tsa_frac)
  if (length(i_tsa)) {
    fsub <- sample(seq_len(p), round(0.4 * p))
    mat[i_tsa, fsub] <- mat[i_tsa, fsub] + pl$bio_tsa_shift * matrix(sdv[fsub], length(i_tsa), length(fsub), byrow = TRUE)
    ts$truth[i_tsa] <- "bio_tsa"
  }


  # features: capped at an upper limit (left-skewed)
  j_trunc <- take_f(pl$trunc_frac)
  for (j in j_trunc) {
    ul <- quantile(mat[, j], pl$trunc_quantile)
    mat[mat[, j] > ul, j] <- ul
  }
  tf$truth[j_trunc] <- "tech_truncated"


  # features: rare biological effect in 8% of samples (right-skewed)
  j_rare <- take_f(pl$rare_frac)
  for (j in j_rare) {
    r <- sample(seq_len(n), round(0.08 * n))
    mat[r, j] <- mat[r, j] + 3 * sdv[j]
  }
  tf$truth[j_rare] <- "bio_rare_effect"


  # features: poor detection (extra missingness applied after the missingness step)
  j_poor <- take_f(pl$poor_feat_frac)
  tf$truth[j_poor] <- "tech_poor_feature"


  # missingness ----
  ref <- mat
  if (miss_ignores_shift) for (r in names(shift_rows)) ref[as.integer(r), ] <- ref[as.integer(r), ] - shift_rows[[r]]
  if (miss_mech != "none" && miss_rate > 0) {
    if (equal_rates) {
      rates <- rep(miss_rate, p)
    } else {
      rates <- pmin(rbeta(p, 1.2, 1.2 / miss_rate - 1.2), 0.6)
      rates <- rates * miss_rate / mean(rates)
      rates <- pmin(rates, 0.75)
    }
    if (complete_frac > 0) rates[sample(p, round(complete_frac * p))] <- 0
    for (j in seq_len(p)) {
      k <- round(rates[j] * n); if (k == 0) next
      if (miss_mech == "mcar") {
        idx <- sample(n, k)
      } else if (miss_mech == "lod") {
        idx <- order(ref[, j])[seq_len(k)]
      } else {
        k1 <- k %/% 2
        idx <- c(order(ref[, j])[seq_len(k1)], sample(setdiff(seq_len(n), order(ref[, j])[seq_len(k1)]), k - k1))
      }
      mat[idx, j] <- NA
    }
  }


  # extra missing values for poor features and failed samples
  for (j in j_poor) {
    k <- round((if (is.null(pl$poor_miss)) runif(1, 0.5, 0.95) else pl$poor_miss) * n)
    if (identical(pl$poor_mech, "lod")) {
      obs <- which(!is.na(mat[, j]))
      idx <- obs[order(ref[obs, j])][seq_len(min(k, length(obs)))]
    } else idx <- sample(n, k)
    mat[idx, j] <- NA
  }
  for (i in i_fail) mat[i, sample(p, round((if (is.null(pl$failed_miss)) runif(1, 0.5, 0.9) else pl$failed_miss) * p))] <- NA


  # return ----
  samples  <- data.frame(sample_id = sid, batch = sample(c("b1", "b2"), n, TRUE), stringsAsFactors = FALSE)
  features <- data.frame(feature_id = fid, platform = "pos", stringsAsFactors = FALSE)
  samples$truth           <- ts$truth
  features$truth          <- tf$truth
  features$cluster_factor <- tf$factor
  if (as_omiprep) return(Omiprep(data = mat, samples = samples, features = features))
  list(data = mat, samples = samples, features = features)
}


# defaults: nothing planted
.simulate_omics_plant_defaults <- function() {
  list(
    dilution_frac     = 0, dilution_shift     = 1.0,
    outlier_frac      = 0, outlier_shift      = 1.0,
    failed_frac       = 0, failed_miss        = NULL,
    bio_subgroup_frac = 0, bio_subgroup_shift = 2,
    bio_tsa_frac      = 0, bio_tsa_shift      = 0.5,
    trunc_frac        = 0, trunc_quantile     = 0.6,
    rare_frac         = 0,
    poor_feat_frac    = 0, poor_miss          = NULL, poor_mech = "random"
  )

}
