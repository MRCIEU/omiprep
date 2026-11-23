test_that("socket and forked back-ends give identical independent features", {
  skip_on_cran()
  set.seed(1)
  m <- matrix(rnorm(200 * 40), 200, dimnames = list(NULL, paste0("F", 1:40)))
  m[, 21:40] <- m[, 1:20] + rnorm(200 * 20, sd = 0.3)       # correlated pairs
  m[sample(length(m), 300)] <- NA
  withr::local_options(omiprep.parallel_backend = "socket")
  for (fast in c(TRUE, FALSE)) {
    sock <- tree_and_independent_features(m, cores = 2, fast = fast)$data
    options(omiprep.parallel_backend = "fork")
    fork <- tree_and_independent_features(m, cores = 2, fast = fast)$data
    one  <- tree_and_independent_features(m, cores = 1, fast = fast)$data
    options(omiprep.parallel_backend = "socket")
    expect_identical(sock, fork)
    expect_identical(sock, one)
  }
})

test_that("the shared socket cluster is started once and stopped", {
  skip_on_cran()
  withr::local_options(omiprep.parallel_backend = "socket")
  expect_true(omiprep:::.start_cluster(2))
  expect_false(omiprep:::.start_cluster(2))                  # already running: not started again
  expect_length(omiprep:::.omiprep_env$cl, 2)
  omiprep:::.stop_cluster()
  expect_null(omiprep:::.omiprep_env$cl)
  expect_false(omiprep:::.start_cluster(1))                  # one core: no cluster
})
