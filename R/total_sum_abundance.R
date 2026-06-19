#' @title Estimates total sum abundance
#' @description
#' This function estimates total sum abundance for numeric data in a matrix, for (1) all features and (2) all features with complete data.
#' @param data matrix, the 'omics data matrix. Samples in rows, features in columns
#' @param ztransform logical, should the feature data be z-transformed and absolute value minimum, mean shifted prior to summing the feature values. TRUE or FALSE.
#'
#' @return a data frame of estimates for (1) total sum abundance and (2) total sum abundance at complete features for each samples. The function also returns the number of features that are not NA at each sample and the number of features with complete data, across all sampmles, used to estimate total sum abundance at complete features.
#'
#' @export
#'
total_sum_abundance <- function(data, ztransform = TRUE){

  # z-transformed data frame
  if(ztransform == TRUE){
    data = apply(data, 2, function(x) {
      ( x - mean(x, na.rm = TRUE) ) / sd(x, na.rm = TRUE)
    })
    ## add absolute(minimum) value to all values
    # cat(paste0("\t\t\t- adding absolute minimum observed value to all values to make all values positive.\n") )
    data = data + abs(min(data, na.rm = TRUE))
  }

  # number of non-NA features
  non_na_feature_count = apply(data, 1, function(x){ sum(!is.na(x)) })
  
  # total sum abundance
  total_tsa = apply(data, 1, function(x){ sum(x, na.rm = TRUE) })

  # find features with complete data
  mis = apply(data, 2, function(x){ sum(is.na(x))/length(x) })
  completefeatures = which(mis == 0)

  # TSA with complete features
  completeF_tsa = apply(data[, completefeatures], 1, function(x){ sum(x, na.rm = TRUE) })

  # output data.frame
  out = data.frame("sample_id" = names(total_tsa), 
                   "non_na_feature_count" = non_na_feature_count, 
                   "tsa_total" = total_tsa, 
                   "tsa_complete_features" = completeF_tsa,
                   "complete_feature_count" = length(completefeatures))

  # return result
  return(out)
}
