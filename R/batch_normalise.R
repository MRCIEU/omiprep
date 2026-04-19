#' @title Batch Normalisation
#' @description
#' Run batch normalisation based on the platform flag in the features data 
#' @param omiprep an object of class Omiprep
#' @param run_mode_col character, column name in features data containing the run mode
#' @param run_mode_colmap named character vector or list, c(mode = "mode col name in samples") 
#' @param source_layer character, which data layer to get the data from
#' @param dest_layer character, which data layer to put the the data in to
#' @include class_omiprep.R
#' @importFrom stats median
#' @export
batch_normalise <- new_generic(name = "batch_normalise", dispatch_args = c("omiprep"), function(omiprep, run_mode_col, run_mode_colmap, source_layer = "input", dest_layer = "batch_normalised") { S7_dispatch() })
#' @name batch_normalise
method(batch_normalise, Omiprep) <- function(omiprep, run_mode_col, run_mode_colmap, source_layer = "input", dest_layer = "batch_normalised") {

  # checks
  run_mode_col <- match.arg(run_mode_col, choices = names(omiprep@features))
  stopifnot("`run_mode_colmap` names must match unique `run_mode_col` entries" = length(setdiff(omiprep@features[[run_mode_col]], names(run_mode_colmap)))==0)
  stopifnot("`run_mode_colmap` values must be columns in sample data" = all(unname(run_mode_colmap) %in% names(omiprep@samples)))
  
  # add a copy of the raw data to the back of the matrix stack
  omiprep@data <- array(
    c(omiprep@data, omiprep@data[, , "input"]),
    dim = c(dim(omiprep@data)[1], dim(omiprep@data)[2], dim(omiprep@data)[3] + 1),
    dimnames = list(
      dimnames(omiprep@data)[[1]],
      dimnames(omiprep@data)[[2]],
      c(dimnames(omiprep@data)[[3]], dest_layer)
    )
  )

  # get the unique flags by which to batch normalise
  for (i in seq_along(run_mode_colmap)) {
    
    mode     <- names(run_mode_colmap)[i]
    mode_col <- run_mode_colmap[[i]]

    f_idx <- which(omiprep@features[[run_mode_col]] == mode)

    batch_ids <- unique(omiprep@samples[[mode_col]])

    for (bid in batch_ids) {

      s_idx <- which(omiprep@samples[[mode_col]] == bid)

      m <- apply(omiprep@data[s_idx, f_idx, source_layer], 2, function(x) stats::median(x, na.rm = TRUE))

      for (j in 1:length(m)) {
        omiprep@data[s_idx, f_idx[j], dest_layer] <- omiprep@data[s_idx, f_idx[j], dest_layer] / m[j]
      }
    }
  }

  return(omiprep)
}
