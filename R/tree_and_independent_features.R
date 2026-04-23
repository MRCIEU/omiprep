.par_lapply <- function(X, FUN, cores, ...) {
  if (.Platform$OS.type == "windows" && cores > 1) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))
    parallel::parLapply(cl, X, FUN, ...)
  } else {
    parallel::mclapply(X, FUN, mc.cores = cores, ...)
  }
}

#' @title Identify Independent Features in a Numeric Matrix
#' @description
#' This function identifies independent features using Spearman's rho correlation distances, and a dendrogram
#' tree cut step.
#' @param data matrix, the 'omics data matrix. samples in row, features in columns
#' @param tree_cut_height the tree cut height. A value of 0.2 (1-Spearman's rho) is equivalent to saying that features with a rho >= 0.8 are NOT independent.
#' @param features_exclude character, vector of feature id indicating features to exclude from the sample and PCA summary analysis but keep in the data
#' @param feature_selection character. Method for selecting a representative feature from each correlated feature cluster. 
#' @param cores number of cores available for parallelism; the default null will try find the maximum available cores - 1; set to 1 for linear, but potentially slow, computation of the correlation matrix. 
#' @param fast If \code{TRUE}, accelerates correlation computation by imputing missing values to the column minimum, pre-ranking all columns, and computing Pearson correlation on ranked data (approximating Spearman). Substantially faster than exact Spearman at large feature dimensions (\eqn{p > 5000}) but assumes missing data are missing at random. Features with high missingness will have inflated rank ties at the median (ensure these are filtered out appropriately with the missingness option). Default \code{FALSE}.
#' One of:
#' \describe{
#'   \item{\code{"max_var_exp"}}{(Default) Selects the feature with the highest sum of absolute Spearman correlations to other features in the cluster; 
#'   effectively the feature explaining the most shared variance.}
#'   \item{\code{"least_missingness"}}{Selects the feature with the fewest missing values within the cluster.}
#' }
#' @return A list with the following components:
#' \describe{
#'   \item{data}{A `data.frame` with:
#'     \itemize{
#'       \item `feature_id`: Feature (column) names from the input matrix.
#'       \item `k`: The cluster index assigned to each feature after tree cutting.
#'       \item `independent_features`: Logical indicator of whether the feature was selected as an independent (representative) feature.
#'     }
#'   }
#'   \item{tree}{A `hclust` object representing the hierarchical clustering of the features based on 1 - |Spearman's rho| distance.}
#' }
#' 
#' @importFrom stats as.dist hclust cutree
#' @importFrom parallel mclapply makeCluster stopCluster parLapply
#'
#' @export
#'
tree_and_independent_features = function(data, tree_cut_height = 0.5, features_exclude = NULL, feature_selection = "max_var_exp", cores = NULL, fast = FALSE){

  # testing
  if (FALSE) {
    tree_cut_height = 0.5
    features_exclude = NULL
    feature_selection <- "least_missingness"
    set.seed(123)
    data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
    colnames(data) <- paste0("F", 1:10)
    data[sample(length(data), 100)] <- NA
    data <- cbind(data, F_dup = data[, "F2"] + rnorm(100, sd = 0.01))
  }
  
  
  # checks 
  feature_selection <- match.arg(feature_selection, choices = c("max_var_exp", "least_missingness"))
  
  
  # find available cores for parallel processing
  if (is.null(cores)) {
    cores <- local({
      slurm <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) # guard against cluster specifying all node cores
      if (!is.na(slurm)) slurm else max(parallel::detectCores() - 1, 1)
    })
  }
  if (!(is.numeric(cores) && length(cores) == 1 && cores >= 1 && cores == as.integer(cores))) {
    stop(sprintf("Cores must be a positive integer, got: %s", deparse(cores)))
  }
  
  # remove excluded features
  if (!is.null(features_exclude)) {
    data <- data[, !colnames(data) %in% features_exclude]
  }
  
  
  # remove features with no variance
  row_var0 <- which( apply(data, 2, function(x) var(x,na.rm=T)==0) )
  if(length(row_var0) > 0){
    data <- data[, -row_var0]
  }
  
  
  # remove features with <80% presence 
  low_presence <- which(colMeans(!is.na(data)) < 0.8)
  if(length(low_presence) > 0){
    data <- data[, -low_presence]
  }
  
  
  # if fast then impute and rank once so we can avoid "pairwise.complete.obs" overhead
  if (fast) {
    col_medians  <- apply(data, 2, min, na.rm = TRUE)
    na_idx       <- which(is.na(data), arr.ind = TRUE)
    data[na_idx] <- col_medians[na_idx[, 2]]
    
    ranked <- apply(data, 2, rank, ties.method = "average")
    
    # parallelism (ML suggestion) - pearsons on imputed ranked = fast
    if (cores == 1) {
      cor_matrix <- stats::cor(ranked, method = "pearson")
    } else {
      idx    <- split(1:ncol(ranked), cut(seq_along(1:ncol(ranked)), breaks=cores, labels=FALSE))
      c_list <- .par_lapply(idx, function(j) {
        stats::cor(ranked[, j, drop=FALSE], ranked, method="pearson")
      }, cores=cores)
      cor_matrix <- do.call(rbind, c_list)
    }
  } else {
    # parallelism (ML suggestion) - spearman on potentially missing = slow due to pairwise comparisons
    if (cores == 1) {
      cor_matrix <- stats::cor(data, method = "spearman",  use = "pairwise.complete.obs")
    } else {
      idx <- split(1:ncol(data), cut(seq_along(1:ncol(data)), breaks=cores, labels=FALSE))
      c_list <- .par_lapply(idx, function(j) {
        stats::cor(data[, j, drop=FALSE], data, method = "spearman",  use = "pairwise.complete.obs")
      }, cores=cores)
      cor_matrix <- do.call(rbind, c_list)
    }
  }
  rownames(cor_matrix) <- colnames(data)
  colnames(cor_matrix) <- colnames(data)
  dist_matrix <- stats::as.dist(1 - abs(cor_matrix))
  stree       <- stats::hclust(dist_matrix, method = "complete")
  
  
  # restrict based on cut off
  k <- stats::cutree(stree, h = tree_cut_height)
  k_group <- table(k)
  

  # keep all single-feature clusters
  ind_k <- names(k_group[k_group == 1])
  ind   <- names(k)[k %in% ind_k]
  
  
  # pick 1 feature within the clusters of size > 1
  cluster_ids <- names(k_group[k_group > 1])
  
  if (feature_selection == "least_missingness") {
    
    N       <- apply(data, 2, function(x){ sum(!is.na(x)) })
    ind2 <- sapply(cluster_ids, function(x){
      w   <- which(k %in% x)
      n   <- names( k[w] )
      o   <- sort(N[n], decreasing = TRUE)
      out <- names(o)[1]
      return(out)
    })
    
  } else if (feature_selection == "max_var_exp") {
    
    # clusters with >1 feature, pick one with highest variance explained
    ind2 <- sapply(cluster_ids, function(cluster) {
      members      <- names(k)[k == as.integer(cluster)]
      sub_cor      <- cor_matrix[members, members]
      cor_sums     <- rowSums(abs(sub_cor), na.rm = TRUE)
      best_feature <- names(which.max(cor_sums))
      return(best_feature)
    })
    
  } else {
    
    stop("feature_selection parameter not recognised - (this shouldn't happen)")
    
  }

  
  independent_features <- paste( c(ind, ind2) )


  # return 
  out <- data.frame("feature_id"           = colnames(data),
                    "k"                    = k[match(colnames(data), names(k))],
                    "independent_features" = colnames(data) %in% independent_features)
  
  return(list(data = out, tree = stree))
}


