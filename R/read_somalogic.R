#' @title Read and Process SomaLogic adat file
#' @description
#' This function reads and processes a commercial SomaLogic `.adat` file. It extracts
#' RFU (Relative Florecent Units) data for samples and controls, along with
#' their respective metadata and feature (protein) metadata. The function returns
#' a structured list suitable for further analysis.
#'
#' @param filepath A string specifying the path to the SomaLogic `.adat` file.
#' @param return_Omiprep logical, if TRUE (default) return a Omiprep object, if FALSE return a list.
#' 
#' @return A omiprep object or a named list with the following elements:
#' \describe{
#'   \item{data}{A matrix of RFU values for experimental samples, with `SampleId`
#'         as row names and `SeqId` columns as column names.}
#'   \item{samples}{A data frame containing metadata for experimental samples, with
#'         `sample_id` renamed from `SampleId`.}
#'   \item{features}{A data frame containing feature-level metadata, including a newly
#'         created `feature_id` column derived from `SeqId`.}
#'   \item{controls}{A matrix of RFU values for control samples, specifically
#'         "Calibrator" samples, with `SampleId` as row names and `SeqId` columns
#'         as column names. If duplicate control `SampleId` values are present,
#'         control `SampleId` is replaced with `paste0(SampleId, "_", SlideId)`.}
#'   \item{control_metadata}{A data frame containing metadata for control samples,
#'         specifically "Calibrator" samples, with `sample_id` renamed from `SampleId`.}
#' }
#'
#' @importFrom SomaDataIO read_adat is.soma_adat
#' @export
read_somalogic <- function(filepath, return_Omiprep = FALSE) {
  
  # testing ====
  if (FALSE) {
    filepath <- system.file(
      "extdata",
      "example_data10.adat",
      package = "SomaDataIO",
      mustWork = TRUE
    )
  }
  
  # checks 1 ====
  if (!grepl("(?i)\\.(adat)$", filepath)) {
    stop(
      "Expected a commercial Somalogic file with extension .adat",
      call. = FALSE
    )
  }
  
  # read data ====
  df <- SomaDataIO::read_adat(file = filepath)
  
  # checks 2 ====
  if (!SomaDataIO::is.soma_adat(df)) {
    stop("The file is not a SomaLogic 'adat' file.", call. = FALSE)
  }
  
  if (!"SampleId" %in% colnames(df)) {
    stop(
      "Column 'SampleId' not found in the dataset. Either your data is not Somalogic or you have renamed your sample column for some reason. Suggest you use an alternative approach to reading in your data (see Vignette XXX).",
      call. = FALSE
    )
  }
  
  if (!"SampleType" %in% colnames(df)) {
    stop(
      "Column 'SampleType' not found in the dataset.",
      call. = FALSE
    )
  }
  
  # make duplicate control SampleId values unique ====
  control_rows <- df$SampleType == "Calibrator"
  control_rows[is.na(control_rows)] <- FALSE
  
  control_sample_ids <- df$SampleId[control_rows]
  
  if (any(duplicated(control_sample_ids))) {
    
    # use Barcode if available, otherwise row names, otherwise 1:n for controls
    if ("Barcode" %in% colnames(df)) {
      
      control_suffix <- df$Barcode[control_rows]
      suffix_source <- "Barcode"
      
    } else {
      
      df_rownames <- rownames(df)
      
      if (!is.null(df_rownames) && length(df_rownames) == nrow(df)) {
        control_suffix <- df_rownames[control_rows]
        suffix_source <- "row names"
      } else {
        control_suffix <- seq_len(sum(control_rows))
        suffix_source <- "1:n control row number"
      }
    }
    
    message(
      paste0(
        "Duplicate control SampleId values detected. ",
        "Creating new control SampleId values using SampleId + ",
        suffix_source,
        "."
      )
    )
    
    if (any(is.na(control_suffix)) || any(!nzchar(as.character(control_suffix)))) {
      stop(
        paste0(
          "Duplicate control SampleId values were found, but some ",
          suffix_source,
          " values are missing or empty. Cannot create unique control SampleId values."
        ),
        call. = FALSE
      )
    }
    
    df$SampleId[control_rows] <- paste0(
      df$SampleId[control_rows],
      "_",
      control_suffix
    )
    
    if (any(duplicated(df$SampleId[control_rows]))) {
      stop(
        paste0(
          "Control SampleId values are still duplicated after combining SampleId with ",
          suffix_source,
          "."
        ),
        call. = FALSE
      )
    }
  }  # init ====
  data <- NULL
  controls <- NULL
  features <- NULL
  samples <- NULL
  control_metadata <- NULL
  
  # feature data for SAMPLES ====
  data <- df[df$SampleType == "Sample", ]
  rownames(data) <- data$SampleId
  data <- data[, grep("seq", colnames(data), value = TRUE)]
  data <- as.matrix(data)
  
  # feature data for CONTROLS ====
  controls <- df[df$SampleType == "Calibrator", ]
  rownames(controls) <- controls$SampleId
  controls <- controls[, grep("seq", colnames(controls), value = TRUE)]
  controls <- as.matrix(controls)
  
  # feature meta-data ====
  features <- attr(df, "Col.Meta")
  features$feature_id <- paste0("seq.", gsub("-", ".", features$SeqId))
  features <- features[, c("feature_id", setdiff(names(features), "feature_id"))]
  
  # sample meta-data ====
  row_meta <- attr(df, "row_meta")
  plain_df <- as.data.frame(df)
  
  samples <- plain_df[, row_meta, drop = FALSE]
  names(samples)[names(samples) == "SampleId"] <- "sample_id"
  samples <- samples[c("sample_id", setdiff(names(samples), "sample_id"))]
  samples <- samples[samples$SampleType == "Sample", ]
  
  control_sample_meta <- plain_df[, row_meta, drop = FALSE]
  names(control_sample_meta)[names(control_sample_meta) == "SampleId"] <- "sample_id"
  control_sample_meta <- control_sample_meta[
    c("sample_id", setdiff(names(control_sample_meta), "sample_id"))
  ]
  control_sample_meta <- control_sample_meta[
    control_sample_meta$SampleType == "Calibrator",
  ]
  
  # return ====
  if (return_Omiprep) {
    return(
      Omiprep(
        data = data,
        samples = samples,
        features = features
      )
    )
  } else {
    return(
      list(
        data = data,
        samples = samples,
        features = features,
        controls = controls,
        control_metadata = control_sample_meta
      )
    )
  }
}
