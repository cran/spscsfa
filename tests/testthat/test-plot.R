test_that("plot method evaluates the smooth coefficients on a grid via ezdat", {

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

  fit <- spscsfa(
    output ~ input | covariate,
    data = dat,
    verbose = FALSE
  )

  # the estimation sample is kept in the fitted object
  expect_s3_class(fit$model_frame, "data.frame")
  expect_equal(nrow(fit$model_frame), n)

  grDevices::pdf(NULL)

  out <- plot(fit, n_grid = 25)

  expect_type(out, "list")
  expect_named(out, "covariate")

  d <- out$covariate
  expect_equal(nrow(d), 25L)
  expect_equal(ncol(d), 1L + ncol(smoothcoef(fit)$theta_beta))
  expect_equal(d$covariate[1L], min(z))
  expect_equal(d$covariate[25L], max(z))
  expect_equal(diff(d$covariate), rep((max(z) - min(z)) / 24, 24L))

  # same as calling spscsfa() with the grid in 'ezdat'
  ref <- spscsfa(
    output ~ input | covariate,
    data = dat,
    ezdat = data.frame(covariate = seq(min(z), max(z), length.out = 25)),
    bandwidth = fit$smooth_coefficients$bandwidth,
    verbose = FALSE
  )
  expect_equal(
    unname(as.matrix(d[, -1L, drop = FALSE])),
    unname(as.matrix(ref$smooth_coefficients$theta_beta)),
    tolerance = 1e-6
  )

  # selecting a smooth coefficient
  out1 <- plot(fit, which = 2, n_grid = 10)
  expect_equal(ncol(out1$covariate), 2L)

  # invalid arguments are rejected
  expect_error(plot(fit, which = 99), "which")
  expect_error(plot(fit, zvar = "unknown"), "zvar")
  expect_error(plot(fit, n_grid = 1), "n_grid")
  expect_error(plot(fit, mfrow = c(0, 1)), "mfrow")

  # an object without the estimation sample cannot be plotted
  fit_old <- fit
  fit_old$model_frame <- NULL
  expect_error(plot(fit_old), "model_frame")

  # no continuous z variable to plot against
  fit_disc <- fit
  fit_disc$model_frame$covariate <- factor(rep(c("a", "b"), length.out = n))
  expect_error(plot(fit_disc), "no continuous")

  grDevices::dev.off()
})

test_that("other z variables are held at their averages and factors are ignored", {

  skip_on_cran()

  set.seed(1)

  n <- 100

  x <- rnorm(n)
  z1 <- runif(n)
  z2 <- runif(n)
  z3 <- factor(rep(c("A", "B"), length.out = n))
  y <- cos(z1 + z2) +
    sin(z1 - z2) * x +
    rnorm(n, 0, 0.1) -
    exp(0.1 + 0.1 * z1) * abs(rnorm(n))

  dat <- data.frame(y, x, z1, z2, z3)

  fit <- spscsfa(y ~ x | z1 + z2 + z3, data = dat, verbose = FALSE)

  grDevices::pdf(NULL)

  out <- plot(fit, n_grid = 20)

  # no plot against the discrete z3
  expect_named(out, c("z1", "z2"))

  # z1 on its grid, z2 at its average, z3 at its most frequent category
  ref <- spscsfa(
    y ~ x | z1 + z2 + z3,
    data = dat,
    ezdat = data.frame(
      z1 = seq(min(z1), max(z1), length.out = 20),
      z2 = mean(z2),
      z3 = factor("A", levels = c("A", "B"))
    ),
    bandwidth = fit$smooth_coefficients$bandwidth,
    verbose = FALSE
  )
  expect_equal(
    unname(as.matrix(out$z1[, -1L, drop = FALSE])),
    unname(as.matrix(ref$smooth_coefficients$theta_beta)),
    tolerance = 1e-6
  )

  # selecting a z variable and a coefficient
  out2 <- plot(fit, zvar = "z2", which = 1, n_grid = 10)
  expect_named(out2, "z2")
  expect_equal(ncol(out2$z2), 2L)

  # a discrete z variable cannot be selected
  expect_error(plot(fit, zvar = "z3"), "discrete")

  grDevices::dev.off()
})
