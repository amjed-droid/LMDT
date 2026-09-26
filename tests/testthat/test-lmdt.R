test_that("lmdt_test returns valid htest object", {
  set.seed(123)
  n <- 60
  x <- rnorm(n)
  y <- 1 + 2 * x + rnorm(n)
  res <- lmdt_test(y ~ x, method = "asymptotic")
  
  expect_s3_class(res, "lmdt")
  expect_s3_class(res, "htest")
  expect_true(res$p.value >= 0 && res$p.value <= 1)
  expect_true(!is.na(res$statistic))
})

test_that("lmdt_graph outputs symmetric normalized matrices", {
  set.seed(123)
  X <- matrix(rnorm(40 * 2), ncol = 2)
  g <- lmdt_graph(X, k = 5)
  
  expect_equal(g$S, t(g$S), tolerance = 1e-10)
  expect_equal(diag(g$S), rep(0, 40))
})
