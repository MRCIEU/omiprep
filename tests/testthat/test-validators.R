library(testthat)

# Tests for the Omiprep S7 class: constructor behaviour, per-property
# validators, and the class-level validator (R/class_omiprep.R).

# ---- shared fixture ---------------------------------------------------------

make_omiprep <- function(n_samples = 5, n_features = 4) {
  data <- matrix(seq_len(n_samples * n_features),
                 nrow = n_samples, ncol = n_features)
  rownames(data) <- paste0("Sample_", seq_len(n_samples))
  colnames(data) <- paste0("Met_", seq_len(n_features))
  samples  <- data.frame(sample_id  = rownames(data))
  features <- data.frame(feature_id = colnames(data))
  Omiprep(data, samples, features)
}

# ---- constructor ------------------------------------------------------------

test_that("constructor promotes a 2D matrix to a 3D 'input' array layer", {
  m <- make_omiprep()
  expect_s3_class(m, "omiprep::Omiprep")
  expect_length(dim(m@data), 3)
  expect_identical(dimnames(m@data)[[3]], "input")
  expect_identical(dim(m@data)[1:2], c(5L, 4L))
  expect_identical(m@data[, , "input"],
                   matrix(seq_len(20), 5, 4,
                          dimnames = list(paste0("Sample_", 1:5),
                                          paste0("Met_", 1:4))))
})

test_that("constructor reorders samples and features to match data", {
  data <- matrix(seq_len(12), nrow = 3, ncol = 4)
  rownames(data) <- c("Sample_1", "Sample_2", "Sample_3")
  colnames(data) <- c("Met_1", "Met_2", "Met_3", "Met_4")
  # supply metadata in a DIFFERENT order than the data
  samples  <- data.frame(sample_id  = c("Sample_3", "Sample_1", "Sample_2"),
                         site        = c("c", "a", "b"))
  features <- data.frame(feature_id = c("Met_4", "Met_3", "Met_2", "Met_1"))
  m <- Omiprep(data, samples, features)

  expect_identical(m@samples[["sample_id"]], rownames(data))
  expect_identical(m@features[["feature_id"]], colnames(data))
  # carried column must be reordered alongside sample_id, not left stale
  expect_identical(m@samples[["site"]], c("a", "b", "c"))
})

test_that("constructor errors when data names are absent from metadata", {
  data <- matrix(seq_len(12), nrow = 3, ncol = 4)
  rownames(data) <- c("Sample_1", "Sample_2", "Sample_3")
  colnames(data) <- c("Met_1", "Met_2", "Met_3", "Met_4")
  good_features <- data.frame(feature_id = colnames(data))

  expect_error(
    Omiprep(data,
            data.frame(sample_id = c("Sample_1", "Sample_2", "WRONG")),
            good_features),
    "sample id\\(s\\) in `data` are absent from `samples\\$sample_id`"
  )
})

test_that("constructor coerces numeric ids to character", {
  # ids that read from file as numeric must not fail the validator's identical() check
  data <- matrix(seq_len(4), nrow = 2, ncol = 2)
  rownames(data) <- c("101", "102")
  colnames(data) <- c("2001", "2002")

  m <- Omiprep(data,
               data.frame(sample_id  = c(101, 102)),
               data.frame(feature_id = c(2001L, 2002L)))

  expect_type(m@samples[["sample_id"]], "character")
  expect_type(m@features[["feature_id"]], "character")
  expect_identical(m@samples[["sample_id"]], rownames(data))
  expect_identical(m@features[["feature_id"]], colnames(data))
})

test_that("constructor gives targeted errors for missing id columns and unnamed data", {
  data <- matrix(seq_len(4), nrow = 2, ncol = 2,
                 dimnames = list(c("s1", "s2"), c("f1", "f2")))
  samples  <- data.frame(sample_id  = c("s1", "s2"))
  features <- data.frame(feature_id = c("f1", "f2"))

  expect_error(Omiprep(data, data.frame(Sample_ID = c("s1", "s2")), features),
               "must contain a column named `sample_id`")
  expect_error(Omiprep(data, samples, data.frame(FeatureID = c("f1", "f2"))),
               "must contain a column named `feature_id`")

  unnamed <- matrix(seq_len(4), nrow = 2, ncol = 2)
  expect_error(Omiprep(unnamed, samples, features),
               "must have rownames \\(sample ids\\) and colnames \\(feature ids\\)")
})

# ---- property validators ----------------------------------------------------

test_that("samples property requires a sample_id column", {
  m <- make_omiprep()
  bad <- data.frame(not_sample_id = m@samples[["sample_id"]])
  expect_error(m@samples <- bad, "sample_id")
})

test_that("features property requires a feature_id column", {
  m <- make_omiprep()
  bad <- data.frame(not_feature_id = m@features[["feature_id"]])
  expect_error(m@features <- bad, "feature_id")
})

test_that("data property rejects a non-array numeric value", {
  m <- make_omiprep()
  expect_error(m@data <- c(1, 2, 3), "numeric matrix")
})

test_that("exclusions property rejects a wrongly-named structure", {
  m <- make_omiprep()
  expect_error(m@exclusions <- list(samples = list(), features = list()),
               "samples.*features|user_excluded")
})

test_that("exclusions property accepts the canonical structure", {
  m <- make_omiprep()
  good <- list(
    samples = list(user_excluded                     = character(),
                   extreme_sample_missingness        = character(),
                   user_defined_sample_missingness   = character(),
                   user_defined_sample_totalpeakarea = character(),
                   user_defined_sample_pca_outlier   = character()),
    features = list(user_excluded                    = character(),
                    extreme_feature_missingness      = character(),
                    user_defined_feature_missingness = character(),
                    user_defined_feature_skewness    = character())
  )
  expect_no_error(m@exclusions <- good)
})

# ---- class-level validator --------------------------------------------------

test_that("class validator enforces feature count equals data columns", {
  m <- make_omiprep()
  expect_error(m@features <- m@features[1:2, , drop = FALSE],
               "@features.*equal.*@data|feature_id.*identical")
})

test_that("class validator enforces sample_id identical to data rownames", {
  m <- make_omiprep()
  bad <- m@samples
  bad[["sample_id"]] <- paste0("X_", bad[["sample_id"]])
  expect_error(m@samples <- bad, "identical to the rownames")
})

test_that("class validator enforces feature_id identical to data colnames", {
  m <- make_omiprep()
  bad <- m@features
  bad[["feature_id"]] <- paste0("X_", bad[["feature_id"]])
  expect_error(m@features <- bad, "identical to the colnames")
})
