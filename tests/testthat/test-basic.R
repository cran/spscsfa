test_that("basic", {

  set.seed(1)

  n <- 100

  x <- rnorm(n)
  z <- runif(n)
  y <- cos(z) +
    sin(z) * x +
    rnorm(n, 0, 0.1) -
    exp(0.1 + 0.1 * z) * abs(rnorm(n))

  dat <- data.frame(
    output = y,
    input = x,
    covariate = z
  )

  fit <- spsc_sfa(
    output ~ input | covariate,
    data = dat,
    optim_method = "optim_BFGS",
    verbose = FALSE
  )

  expect_s3_class(fit, "spsc_sfa")
  expect_type(fit, "list")
  expect_true(is.list(fit$smooth_coefficients))
  expect_true(is.list(fit$mle_results))

})
