## ── Validate correlation matrix approach for independent feature selection ────
## Checks that the new median-impute → rank → Rfast::cora() path selects the
## same independent features as the original stats::cor(method="spearman",
## use="pairwise.complete.obs") path.
##
## 350 features is used because the old pairwise Spearman path was used below
## ~300 features; 350 sits just above that, in the range where the new approach
## is expected to be in use.
##
## Usage: Rscript scripts/validate_cor_feature_selection.R
##        or source() with metaboprep loaded.
## ─────────────────────────────────────────────────────────────────────────────

library(metaboprep)
library(Rfast)

## ── 1. Synthetic data ─────────────────────────────────────────────────────────
## Reuse the same correlated-cluster structure as the timing benchmark so the
## data is realistic (30 % correlated features, 15 % missingness).

make_raw_matrix <- function(n_samples, n_features, seed = 42) {
  set.seed(seed)

  sample_ids  <- paste0("s_", seq_len(n_samples))
  feature_ids <- paste0("f_", seq_len(n_features))

  n_clustered   <- floor(n_features * 0.30)
  cluster_sizes <- sample(3:8, ceiling(n_clustered / 5), replace = TRUE)
  cluster_sizes <- cluster_sizes[cumsum(cluster_sizes) <= n_clustered]
  n_factors     <- length(cluster_sizes)

  latent <- matrix(rnorm(n_samples * n_factors), nrow = n_samples, ncol = n_factors)
  mat    <- matrix(rnorm(n_samples * n_features), nrow = n_samples, ncol = n_features)

  col <- 1L
  for (k in seq_len(n_factors)) {
    loading <- runif(cluster_sizes[k], 0.6, 0.9)
    for (j in seq_len(cluster_sizes[k])) {
      mat[, col] <- loading[j] * latent[, k] + sqrt(1 - loading[j]^2) * mat[, col]
      col <- col + 1L
    }
  }

  dimnames(mat) <- list(sample_ids, feature_ids)

  # ~15 % missingness
  n_miss <- floor(0.15 * n_samples * n_features)
  mat[sample(length(mat), n_miss)] <- NA

  mat
}

N_SAMPLES  <- 500
N_FEATURES <- 350
TREE_CUT   <- 0.5

cat(sprintf("Synthetic data: %d samples × %d features  |  tree_cut_height = %.2f\n\n",
            N_SAMPLES, N_FEATURES, TREE_CUT))

mat <- make_raw_matrix(N_SAMPLES, N_FEATURES)


## ── 2. Shared pre-processing (mirrors tree_and_independent_features) ──────────
## Both approaches start from the same cleaned matrix.

# remove zero-variance features
var0 <- which(apply(mat, 2, function(x) var(x, na.rm = TRUE) == 0))
if (length(var0) > 0) mat <- mat[, -var0]

# remove features with <80 % presence
low_pres <- which(colMeans(!is.na(mat)) < 0.8)
if (length(low_pres) > 0) mat <- mat[, -low_pres]

cat(sprintf("Features after pre-processing: %d\n\n", ncol(mat)))


## ── 3. Helper: hclust → cutree → representative feature selection ─────────────

select_independent <- function(cor_matrix, mat, tree_cut_height, feature_selection = "max_var_exp") {
  dist_matrix <- stats::as.dist(1 - abs(cor_matrix))
  stree       <- stats::hclust(dist_matrix, method = "complete")
  k           <- stats::cutree(stree, h = tree_cut_height)
  k_group     <- table(k)

  # single-feature clusters are independent by definition
  ind_k <- names(k_group[k_group == 1])
  ind   <- names(k)[k %in% ind_k]

  # pick representative from multi-feature clusters
  cluster_ids <- names(k_group[k_group > 1])

  if (feature_selection == "max_var_exp") {
    ind2 <- sapply(cluster_ids, function(cluster) {
      members      <- names(k)[k == as.integer(cluster)]
      sub_cor      <- cor_matrix[members, members]
      cor_sums     <- rowSums(abs(sub_cor), na.rm = TRUE)
      names(which.max(cor_sums))
    })
  } else {
    N    <- apply(mat, 2, function(x) sum(!is.na(x)))
    ind2 <- sapply(cluster_ids, function(cluster) {
      members <- names(k)[k == as.integer(cluster)]
      names(sort(N[members]))[1]
    })
  }

  sort(unique(c(ind, ind2)))
}


## ── 4. Old approach: pairwise Spearman via stats::cor ────────────────────────

cat("Running OLD approach  (stats::cor, method = 'spearman', use = 'pairwise.complete.obs') ...\n")
t_old <- system.time({
  cor_old <- stats::cor(mat, method = "spearman", use = "pairwise.complete.obs")
  selected_old <- select_independent(cor_old, mat, TREE_CUT)
})
cat(sprintf("  → %.2f s  |  %d independent features selected\n\n",
            t_old["elapsed"], length(selected_old)))


## ── 5. New approach: median-impute → rank → Rfast::cora ──────────────────────

cat("Running NEW approach  (median impute → rank → Rfast::cora) ...\n")
t_new <- system.time({
  mat_imp  <- apply(mat, 2, function(x) { x[is.na(x)] <- median(x, na.rm = TRUE); x })
  ranked   <- apply(mat_imp, 2, rank, ties.method = "average")
  cor_new  <- Rfast::cora(ranked)
  rownames(cor_new) <- colnames(mat)
  colnames(cor_new) <- colnames(mat)
  selected_new <- select_independent(cor_new, mat, TREE_CUT)
})
cat(sprintf("  → %.2f s  |  %d independent features selected\n\n",
            t_new["elapsed"], length(selected_new)))


## ── 6. Comparison ─────────────────────────────────────────────────────────────

only_old <- setdiff(selected_old, selected_new)
only_new <- setdiff(selected_new, selected_old)
in_both  <- intersect(selected_old, selected_new)

cat("── Feature selection agreement ───────────────────────────────────────────\n")
cat(sprintf("  Selected by both        : %d\n", length(in_both)))
cat(sprintf("  Only in OLD (dropped)   : %d\n", length(only_old)))
cat(sprintf("  Only in NEW (added)     : %d\n", length(only_new)))
cat(sprintf("  Agreement               : %.1f %%\n\n",
            100 * length(in_both) / length(union(selected_old, selected_new))))

if (length(only_old) > 0) {
  cat("Features selected ONLY by old approach:\n")
  print(only_old)
  cat("\n")
}
if (length(only_new) > 0) {
  cat("Features selected ONLY by new approach:\n")
  print(only_new)
  cat("\n")
}

## ── 7. Correlation matrix agreement ──────────────────────────────────────────
## Check how similar the two correlation matrices are, to distinguish
## near-identical results from cases where the approaches genuinely differ.

cat("── Correlation matrix comparison ─────────────────────────────────────────\n")
cor_diff  <- cor_old - cor_new
upper_idx <- upper.tri(cor_diff)
diffs     <- cor_diff[upper_idx]

cat(sprintf("  Max absolute difference : %.6f\n", max(abs(diffs))))
cat(sprintf("  Mean absolute difference: %.6f\n", mean(abs(diffs))))
cat(sprintf("  Pearson r (old vs new)  : %.8f\n\n",
            cor(cor_old[upper_idx], cor_new[upper_idx])))

## ── 8. Cluster-level comparison ───────────────────────────────────────────────
## For any discrepant features, show which cluster they belonged to under each
## approach, to help diagnose whether the difference is a near-tie broken
## differently or a genuine re-clustering.

if (length(only_old) > 0 || length(only_new) > 0) {
  discrepant <- union(only_old, only_new)

  dist_old <- stats::as.dist(1 - abs(cor_old))
  k_old    <- stats::cutree(stats::hclust(dist_old, method = "complete"), h = TREE_CUT)

  dist_new <- stats::as.dist(1 - abs(cor_new))
  k_new    <- stats::cutree(stats::hclust(dist_new, method = "complete"), h = TREE_CUT)

  cat("── Cluster membership for discrepant features ────────────────────────────\n")
  diag_df <- data.frame(
    feature    = discrepant,
    cluster_old = k_old[discrepant],
    cluster_new = k_new[discrepant],
    in_old      = discrepant %in% selected_old,
    in_new      = discrepant %in% selected_new
  )
  print(diag_df, row.names = FALSE)
  cat("\n")
}

## ── 9. Summary verdict ────────────────────────────────────────────────────────

if (length(only_old) == 0 && length(only_new) == 0) {
  cat("PASS: both approaches select identical independent features.\n")
} else {
  cat(sprintf("DIFF: %d feature(s) differ between approaches — review output above.\n",
              length(only_old) + length(only_new)))
}
