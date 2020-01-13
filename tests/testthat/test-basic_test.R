context("Basic tests")

X <- matrix(rnorm(500*100), nrow = 500, ncol = 100)
Y <- sample(c("0", "1"), size = 500, replace = TRUE, prob = c(0.5, 0.5))
Y <- Y == "1"
testthat::test_that("First test: logreghd", {
  testthat::expect_silent(reg_log_hd(X, Y))
})