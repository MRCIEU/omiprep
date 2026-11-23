test_that("simulate_omics returns the documented structure", {
  sim <- simulate_omics(n = 60, p = 40, seed = 1)
  expect_named(sim, c("data", "samples", "features"))
  expect_named(sim$samples, c("sample_id", "batch", "truth"))
  expect_named(sim$features, c("feature_id", "platform", "truth", "cluster_factor"))
  expect_equal(dim(sim$data), c(60L, 40L))
  expect_identical(rownames(sim$data), sim$samples$sample_id)
  expect_identical(colnames(sim$data), sim$features$feature_id)
  expect_true(all(sim$samples$truth == "clean"))
  expect_true(all(sim$features$truth == "clean"))
  expect_s7_class(Omiprep(data = sim$data, samples = sim$samples, features = sim$features), Omiprep)

  m <- simulate_omics(n = 60, p = 40, seed = 1, as_omiprep = TRUE)
  expect_s7_class(m, Omiprep)
  expect_identical(m@samples$truth, sim$samples$truth)
  expect_identical(m@features$truth, sim$features$truth)
})

test_that("simulate_omics is reproducible for a given seed", {
  a <- simulate_omics(n = 50, p = 30, seed = 7, plant = list(failed_frac = 0.1, trunc_frac = 0.1))
  b <- simulate_omics(n = 50, p = 30, seed = 7, plant = list(failed_frac = 0.1, trunc_frac = 0.1))
  expect_identical(a, b)
  expect_false(identical(a$data, simulate_omics(n = 50, p = 30, seed = 8)$data))
})

test_that("missingness follows the requested rate and mechanism", {
  none <- simulate_omics(n = 200, p = 50, miss_mech = "none", seed = 1)
  expect_false(anyNA(none$data))

  eq <- simulate_omics(n = 200, p = 50, miss_rate = 0.1, miss_mech = "mcar", equal_rates = TRUE, seed = 1)
  expect_equal(unname(colMeans(is.na(eq$data))), rep(0.1, 50))

  # below a detection limit: the missing values of a feature are its lowest values
  full <- simulate_omics(n = 200, p = 50, miss_mech = "none", seed = 1)$data
  lod  <- simulate_omics(n = 200, p = 50, miss_rate = 0.1, miss_mech = "lod", equal_rates = TRUE, seed = 1)$data
  expect_true(all(vapply(seq_len(50), function(j) max(full[is.na(lod[, j]), j]) < min(lod[, j], na.rm = TRUE), logical(1))))
})

test_that("planted problems are labelled and have the intended effect", {
  sim <- simulate_omics(n = 200, p = 100, miss_rate = 0.05, seed = 2,
                        plant = list(failed_frac = 0.05, failed_miss = 0.6, dilution_frac = 0.05,
                                     poor_feat_frac = 0.1, poor_miss = 0.6, trunc_frac = 0.1, rare_frac = 0.1))
  ts <- table(sim$samples$truth); tf <- table(sim$features$truth)
  expect_equal(ts[["tech_failed"]], 10); expect_equal(ts[["tech_dilution"]], 10)
  expect_equal(tf[["tech_poor_feature"]], 10); expect_equal(tf[["tech_truncated"]], 10); expect_equal(tf[["bio_rare_effect"]], 10)

  s_mis <- rowMeans(is.na(sim$data)); f_mis <- colMeans(is.na(sim$data))
  expect_true(all(s_mis[sim$samples$truth == "tech_failed"] >= 0.5))
  expect_true(all(f_mis[sim$features$truth == "tech_poor_feature"] >= 0.5))
  expect_true(all(f_mis[sim$features$truth == "clean"] < 0.5))
})

test_that("simulate_omics rejects invalid input", {
  expect_error(simulate_omics(n = 1, p = 10), "`n`")
  expect_error(simulate_omics(n = 10, p = 10, miss_rate = 1), "miss_rate")
  expect_error(simulate_omics(n = 10, p = 10, plant = list(nonsense = 1)), "Unknown `plant` entries")
  expect_error(simulate_omics(n = 10, p = 10, plant = list(failed_frac = 0.6, dilution_frac = 0.6)), "more than all samples")
})

test_that("an empty plant list is accepted and plants nothing", {
  a <- simulate_omics(60, 40, plant = list(), seed = 2)
  b <- simulate_omics(60, 40, seed = 2)
  expect_identical(a, b)
  expect_true(all(a$samples$truth == "clean"))
})
