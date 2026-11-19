test_that("extractors reject wrong class", {
  expect_error(te(list()), "class 'spscsfa'")
  expect_error(inefficiency(list()), "class 'spscsfa'")
  expect_error(smoothcoef(list()), "class 'spscsfa'")
})

test_that("public functions exist", {
  expect_true(is.function(spscsfa))
  expect_true(is.function(te))
  expect_true(is.function(inefficiency))
  expect_true(is.function(smoothcoef))
})

test_that("the old names are aliases of the new ones", {
  expect_true(identical(spsc_sfa, spscsfa))
  expect_true(identical(smooth_coef, smoothcoef))
})

test_that("fits and summaries of the old class 'spsc_sfa' still work", {

  skip_on_cran()

  set.seed(1)
  n <- 60
  x <- rnorm(n)
  z <- runif(n)
  y <- cos(z) + sin(z) * x + rnorm(n, 0, 0.1) -
    exp(0.1 + 0.1 * z) * abs(rnorm(n))
  dat <- data.frame(y, x, z)

  # the old function name, with the new class
  fit <- suppressWarnings(spsc_sfa(y ~ x | z, data = dat, verbose = FALSE))
  expect_s3_class(fit, "spscsfa")

  # a fit saved with the old class: no new estimation is needed
  old <- fit
  class(old) <- "spsc_sfa"

  expect_equal(te(old), te(fit))
  expect_equal(inefficiency(old), inefficiency(fit))
  expect_equal(smooth_coef(old), smoothcoef(fit))
  expect_s3_class(summary(old), "summary.spscsfa")
  expect_identical(capture.output(print(old)), capture.output(print(fit)))

  # a summary saved with the old class
  s_old <- summary(fit)
  class(s_old) <- "summary.spsc_sfa"
  expect_identical(capture.output(print(s_old)),
                   capture.output(print(summary(fit))))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  curves <- plot(old, n_grid = 10)
  expect_named(curves, "z")
})
