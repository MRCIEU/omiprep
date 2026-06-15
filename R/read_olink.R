#' @title Read and Process Olink NPX Data File
#' @description
#'  This function reads and processes an Olink NPX file in long format. It supports `.csv`, `.xls`, `.xlsx`, `.txt`, `.zip`, and `.parquet` formats, using Olink's own OlinkAnalyze::read_NPX() function, and returns a omiprep object or a list of matrices and metadata frames for further analysis.
#'
#' @param filepath A string specifying the path to the Olink NPX file.
#' @param return_Omiprep logical, if TRUE return a Omiprep object, if FALSE (default) return a list.
#' 
#' @returns Omiprep object or a named list with the following elements:
#' \describe{
#'   \item{data}{A matrix of NPX values with `SampleID` as rows and `OlinkID` as columns, containing only sample data.}
#'   \item{samples}{A `data.frame` of sample metadata, one row per sample. Columns that are constant within a sample are retained; the per-assay `QC_Warning` flag is summarised into `qc_n_warnings` (count of flagged assays) and `qc_any_warning` (logical).}
#'   \item{features}{A `data.frame` containing feature-level metadata for samples.}
#'   \item{controls}{A matrix of NPX values for control samples.}
#'   \item{control_metadata}{A `data.frame` containing metadata for control samples.}
#' }
#'
#' @details
#' The function checks whether the input data is in long format by verifying the presence of duplicate `SampleID` values. It also accommodates two variants of Olink files:
#' \itemize{
#'   \item Files that include a `Sample_Type` column with values `"SAMPLE"` and `"CONTROL"`.
#'   \item Files that use the `SampleID` column to label control samples (e.g., entries containing `"CONTROL"`).
#' }
#'
#' If neither format is detected, the function stops with an error indicating that the data is likely not from Olink.
#'
#' @examples
#' \dontrun{
#'   filepath <- system.file("extdata", "example_olink_data.txt", package = "omiprep")
#'   olink_data <- read_olink(filepath)
#' }
#' @importFrom OlinkAnalyze read_NPX
#' @importFrom reshape2 dcast
#' @export

read_olink <- function(filepath, return_Omiprep = FALSE) {
  
  # testing ====
  if (FALSE) {
    filepath <- system.file("extdata", "olink_v1_example.txt", package = "omiprep")
  }
  
  
  # checks 1 ====
  if (!grepl("(?i)\\.(csv|xls|xlsx|zip|txt|parquet)$", filepath)) {
    stop(paste0("Expected a commercial Olink file with extension .csv, .xls, .xlsx, .txt, .zip, or .parquet"), call. = FALSE)
  }
  
  
  # read data ====
  df <- OlinkAnalyze::read_NPX(filepath)
  if (!"SampleID" %in% colnames(df)) {
    stop("Column 'SampleID' not found in the dataset. Either your data is not Olink or you have renamed your sample column for some reason. Suggest you use an alternative approach to reading in your data (see Vignette XXX).", call. = FALSE)
  }
  
  
  # checks 2 ====
  if (!any(duplicated(df$SampleID))) {
    stop("The dataset does not appear to be in long format - no duplicate SampleID values found. Suggest you use an alternative approach to reading in your data (see Vignette XXX).", call. = FALSE)
  }
  
  
  # init ====
  data <- NULL
  controls <- NULL
  features <- NULL
  samples <- NULL
  control_metadata <- NULL
  
  
  # identify sample type column ====
  sample_type_col <- intersect(c("Sample_Type", "SampleType"), colnames(df))
  
  # feature data for SAMPLES ====
  ## olink have multiple versions, some have a Sample_Type/SampleType column (SAMPLE, CONTROL) some don't and instead have "CONTROL" in SampleID column
  if (length(sample_type_col) > 0) {
    
    # Use the first one found
    st_col <- sample_type_col[1]
    df_features_samples <- df[df[[st_col]] == "SAMPLE" & !grepl("empty well", df$SampleID, ignore.case = TRUE), ]
    
  } else if (any(grepl("CONTROL", df$SampleID, ignore.case = TRUE))) {
    
    df_features_samples <- df[!grepl("CONTROL", df$SampleID, ignore.case = TRUE) & !grepl("empty well", df$SampleID, ignore.case = TRUE), ]
    
  } else {
    
    stop("The dataset does not contain a 'Sample_Type' or 'SampleType' column or any 'CONTROL' entries in 'SampleID'. It is likely not an Olink file. Suggest you use an alternative approach to reading in your data (see Vignette XXX)", call. = FALSE)
  
  }

  data           <- reshape2::dcast(df_features_samples, SampleID ~ OlinkID, value.var = "NPX")
  rownames(data) <- data$SampleID
  data           <- as.matrix(data[ , setdiff(colnames(data), "SampleID") ])
    
  
  # feature data for CONTROLS ====
  ## olink have multiple versions, some have a Sample_Type/SampleType column (SAMPLE, CONTROL/SAMPLE_CONTROL) some don't and instead have "CONTROL" in SampleID column
  if (length(sample_type_col) > 0) {
    
    st_col <- sample_type_col[1]
    df_features_controls <- df[(df[[st_col]] %in% c("CONTROL", "SAMPLE_CONTROL") | grepl("SAMPLE_CONTROL", df[[st_col]], ignore.case = TRUE)) | grepl("empty well", df$SampleID, ignore.case = TRUE), ]
    
  } else if (any(grepl("CONTROL", df$SampleID, ignore.case = TRUE))) {
    
    df_features_controls <- df[grepl("CONTROL", df$SampleID, ignore.case = TRUE) | grepl("empty well", df$SampleID, ignore.case = TRUE), ]
  
  } else {
    
    stop("The dataset does not contain a 'Sample_Type' or 'SampleType' column or any 'CONTROL' entries in 'SampleID'. It is likely not an Olink file. Suggest you use an alternative approach to reading in your data (see Vignette XXX)", call. = FALSE)
  
  }
  
  controls           <- reshape2::dcast(df_features_controls, SampleID ~ OlinkID, value.var = "NPX")
  rownames(controls) <- controls$SampleID
  controls_matrix    <- as.matrix(controls[ , setdiff(colnames(controls), "SampleID") ])
  
  
  # feature meta-data ====
  features <- df_features_samples[
    , setdiff(colnames(df_features_samples), c("SampleID", "NPX", "QC_Warning", "Index", "PlateID")), 
    drop = FALSE
  ]
  features <- unique(features)
  names(features)[names(features) == "OlinkID"] <- "feature_id"
  
  ### control_feature_meta is not returned.
  # control_feature_meta <- df_features_controls[
  #   , setdiff(colnames(df_features_controls), c("SampleID", "NPX", "QC_Warning", "Index", "PlateID")), 
  #   drop = FALSE
  # ]
  # control_feature_meta <- unique(control_feature_meta)
  # names(control_feature_meta)[names(control_feature_meta) == "OlinkID"] <- "feature_id"
  
  
  # sample meta-data ====
  # Olink data is long (one row per sample x assay). 
  # Drop the per-assay columns, then collapse to one row per sample. 
  # The per-assay QC_Warning flag is summarised into qc_n_warnings / qc_any_warning
  # see collapse_olink_samples()
  per_assay_cols <- c("Index", "OlinkID", "UniProt", "Assay", "Assay_Warning",
                      "MissingFreq", "LOD", "NPX", "Panel", "Panel_Version",
                      "Normalization")
  samples <- collapse_olink_samples(
    df_features_samples[, setdiff(colnames(df_features_samples), per_assay_cols), drop = FALSE]
  )
  control_sample_meta <- collapse_olink_samples(
    df_features_controls[, setdiff(colnames(df_features_controls), per_assay_cols), drop = FALSE]
  )
  
  
  # return ====
  if (return_Omiprep) {
    return(Omiprep(data = data, 
                      samples = samples, 
                      features = features))
  } else {
    return(list(data = data,
                samples = samples,
                features = features,
                controls = controls_matrix,
                control_metadata = control_sample_meta))
  }

}


#' Collapse an Olink long-format metadata table to one row per sample
#'
#' Olink NPX data has one row per sample x assay. After the per-assay feature
#' columns are removed, this collapses the remaining rows to one row per
#' `SampleID`: columns that are constant within a sample are retained as-is,
#' the per-assay `QC_Warning` flag (if present) is summarised into
#' `qc_n_warnings` (count of assays flagged) and `qc_any_warning` (logical),
#' and any other column that still varies within a sample is dropped with a
#' warning because it cannot be represented at the sample level.
#'
#' @param df data.frame with a `SampleID` column and the per-assay feature
#'   columns already removed.
#' @param id_col character, the sample identifier column. Default `"SampleID"`.
#' @returns data.frame with one row per sample and a leading `sample_id` column.
#' @noRd
collapse_olink_samples <- function(df, id_col = "SampleID") {

  ids       <- as.character(df[[id_col]])
  uid       <- unique(ids)
  idx_by_id <- split(seq_len(nrow(df)), factor(ids, levels = uid))
  first_idx <- vapply(idx_by_id, function(i) i[1], integer(1))

  out <- data.frame(sample_id = uid, stringsAsFactors = FALSE)

  other_cols <- setdiff(colnames(df), id_col)

  # summarise the per-assay QC flag to the sample level, if present
  if ("QC_Warning" %in% other_cols) {
    is_warn <- !is.na(df[["QC_Warning"]]) &
      toupper(trimws(as.character(df[["QC_Warning"]]))) != "PASS"
    out[["qc_n_warnings"]]  <- vapply(idx_by_id, function(i) sum(is_warn[i]), integer(1))
    out[["qc_any_warning"]] <- out[["qc_n_warnings"]] > 0L
    other_cols <- setdiff(other_cols, "QC_Warning")
  }

  # keep columns constant within a sample; drop (with a warning) any that vary
  for (cc in other_cols) {
    vals       <- df[[cc]]
    n_distinct <- vapply(idx_by_id, function(i) length(unique(vals[i])), integer(1))
    if (all(n_distinct <= 1L)) {
      out[[cc]] <- unname(vals[first_idx])
    } else {
      warning(sprintf(
        "Olink column '%s' varies within a sample and was dropped from sample metadata.",
        cc), call. = FALSE)
    }
  }

  rownames(out) <- NULL
  out
}
