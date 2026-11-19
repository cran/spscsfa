# Tests of the spscreg and gradients arguments of spscsfa() and of grad()

make_data <- function(n = 100, seed = 1) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  z1 <- runif(n)
  z2 <- runif(n)
  z3 <- factor(rep(c("A", "B", "C"), length.out = n))
  y <- cos(z1 + z2) + sin(z1 - z2) * x1 + sin(z1 + z2) * x2 +
    rnorm(n, 0, 0.1) - exp(0.1 + 0.1 * z1 + 0.1 * z2) * abs(rnorm(n))
  data.frame(y, x1, x2, z1, z2, z3)
}

test_that("spscreg is checked and gradients are not returned by default", {

  dat <- make_data(60)

  expect_error(
    spscsfa(y ~ x1 | z1, data = dat, spscreg = "foo", verbose = FALSE),
    "should be one of"
  )
  expect_error(
    spscsfa(y ~ x1 | z1, data = dat, gradients = NA, verbose = FALSE),
    "gradients"
  )

  skip_on_cran()

  fit <- suppressWarnings(
    spscsfa(y ~ x1 | z1, data = dat, verbose = FALSE)
  )

  # local-constant is the default and nothing about gradients is returned
  expect_equal(fit$smooth_coefficients$spscreg, "lc")
  expect_equal(fit$model_summary$spscreg, "lc")
  expect_false(fit$model_summary$gradients)
  expect_null(fit$smooth_coefficients_grad)
  expect_null(fit$inefficiency_grad)
  expect_error(grad(fit), "no gradients")

  out <- capture.output(print(fit))
  expect_false(any(grepl("GRADIENTS", out)))
  expect_true(any(grepl("local-constant", out)))
})

test_that("local-constant gradients: layout, alpha and grad()", {

  skip_on_cran()

  n <- 100
  dat <- make_data(n)

  fit <- suppressWarnings(
    spscsfa(y ~ x1 + x2 | z1 + z2 + z3, data = dat, gradients = TRUE,
             verbose = FALSE)
  )

  scg <- fit$smooth_coefficients_grad
  expect_named(scg, c("theta_beta_grad", "alpha_grad"))

  # theta_beta_grad: one matrix for each continuous z variable
  tbg <- scg$theta_beta_grad
  expect_named(tbg, c("z1", "z2"))
  for (zn in c("z1", "z2")) {
    expect_equal(dim(tbg[[zn]]), c(n, 3L))
    expect_equal(
      colnames(tbg[[zn]]),
      paste0("d(", c("Intercept (theta)", "x1", "x2"), ")/d(", zn, ")")
    )
  }

  # these are the gradients of lcscoefgrad() for the bandwidths used
  ref <- lcscoefgrad(fit$smooth_coefficients$bandwidth,
                     txdat = dat[c("x1", "x2")], tydat = dat$y,
                     tzdat = dat[c("z1", "z2", "z3")])
  expect_equal(unname(tbg$z1), unname(ref$z1), tolerance = 1e-6)
  expect_equal(unname(tbg$z2), unname(ref$z2), tolerance = 1e-6)

  # alpha_grad: one column for each continuous z variable
  ag <- scg$alpha_grad
  expect_equal(dim(ag), c(n, 2L))
  expect_equal(
    colnames(ag),
    c("d(Intercept (alpha))/d(z1)", "d(Intercept (alpha))/d(z2)")
  )

  # alpha_grad = theta_grad + E(u|z) * delta_1l
  tab <- fit$mle_results$delta_table
  eu <- as.vector(fit$mle_results$Eu_mle)
  for (zn in c("z1", "z2")) {
    delta <- tab$Estimate[tab$Coefficient == zn]
    expect_equal(
      unname(ag[, paste0("d(Intercept (alpha))/d(", zn, ")")]),
      unname(tbg[[zn]][, 1L] + eu * delta),
      tolerance = 1e-8
    )
  }

  # inefficiency_grad: one column for each inefficiency determinant, the
  # factor entering through its binary variables
  expect_equal(nrow(fit$inefficiency_grad), n)
  expect_equal(colnames(fit$inefficiency_grad), c("z1", "z2", "z3_B", "z3_C"))

  # grad() returns the three elements
  g <- grad(fit)
  expect_named(g, c("theta_beta_grad", "alpha_grad", "inefficiency_grad"))
  expect_identical(g$theta_beta_grad, tbg)
  expect_identical(g$alpha_grad, ag)
  expect_identical(g$inefficiency_grad, fit$inefficiency_grad)

  # ... for one continuous z variable
  g1 <- grad(fit, zvar = "z1")
  expect_named(g1$theta_beta_grad, "z1")
  expect_equal(ncol(g1$alpha_grad), 1L)
  expect_equal(colnames(g1$inefficiency_grad), "z1")

  # ... and for a factor: only inefficiency_grad
  expect_message(g3 <- grad(fit, zvar = "z3"), "not a continuous z variable")
  expect_null(g3$theta_beta_grad)
  expect_null(g3$alpha_grad)
  expect_equal(colnames(g3$inefficiency_grad), c("z3_B", "z3_C"))

  expect_error(grad(fit, zvar = "nope"), "unknown variable")

  # print and summary report the gradients, whatever 'verbose' is
  s <- summary(fit)
  expect_named(s$gradients,
               c("theta_beta_grad", "alpha_grad", "inefficiency_grad"))
  expect_named(s$gradients$theta_beta_grad, c("z1", "z2"))
  expect_equal(colnames(s$gradients$alpha_grad),
               c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."))

  out <- capture.output(print(fit))
  expect_true(any(grepl("GRADIENTS", out)))

  # the summary of alpha comes with the results of step 2
  expect_equal(rownames(s$alpha), "Intercept (alpha)")
  expect_equal(colnames(s$alpha),
               c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."))
  i_step2 <- grep("STEP 2", out)
  i_alpha <- grep("Intercept \\(alpha\\) Summary \\(over the training", out)
  expect_length(i_alpha, 1L)
  expect_gt(i_alpha, i_step2)
})

test_that("local-constant gradients with evaluation data", {

  skip_on_cran()

  n <- 100
  dat <- make_data(n)
  ez <- data.frame(z1 = seq(0.1, 0.9, length.out = 25), z2 = 0.5)

  fit <- suppressWarnings(
    spscsfa(y ~ x1 | z1 + z2, data = dat, ezdat = ez, gradients = TRUE,
             verbose = FALSE)
  )

  # gradients of theta and beta at the evaluation points, alpha_grad at the
  # training data, as alpha itself
  expect_equal(nrow(fit$smooth_coefficients$theta_beta), 25L)
  expect_equal(nrow(fit$smooth_coefficients_grad$theta_beta_grad$z1), 25L)
  expect_equal(nrow(fit$smooth_coefficients$alpha), n)
  expect_equal(nrow(fit$smooth_coefficients_grad$alpha_grad), n)
})

test_that("local-linear smooth coefficients and their gradients", {

  skip_on_cran()

  n <- 80
  dat <- make_data(n)

  fit <- suppressWarnings(
    spscsfa(y ~ x1 | z1 + z2, data = dat, spscreg = "ll", gradients = TRUE,
             verbose = FALSE)
  )

  expect_s3_class(fit, "spscsfa")
  expect_equal(fit$smooth_coefficients$spscreg, "ll")
  expect_s3_class(fit$smooth_coefficients$bandwidth, "llsmoothcoefficientbw")
  expect_equal(dim(fit$smooth_coefficients$theta_beta), c(n, 2L))
  expect_equal(colnames(fit$smooth_coefficients$theta_beta),
               c("Intercept (theta)", "x1"))

  # the smooth coefficients and the residuals are those of llscoef()
  bw <- fit$smooth_coefficients$bandwidth
  ref <- llscoef(bw, data = dat, betas = TRUE,
                 residuals = TRUE, gradients = TRUE)
  expect_equal(unname(fit$smooth_coefficients$theta_beta), unname(coef(ref)),
               tolerance = 1e-6)
  expect_equal(unname(fit$smooth_coefficients$residuals),
               unname(residuals(ref)), tolerance = 1e-6)

  # llscoef() gives d beta_1 / d z_1, d beta_1 / d z_2, d beta_2 / d z_1, ...
  # and spscsfa() one matrix for each continuous z variable
  tbg <- fit$smooth_coefficients_grad$theta_beta_grad
  expect_named(tbg, c("z1", "z2"))
  expect_equal(colnames(tbg$z1), c("d(Intercept (theta))/d(z1)", "d(x1)/d(z1)"))
  expect_equal(colnames(tbg$z2), c("d(Intercept (theta))/d(z2)", "d(x1)/d(z2)"))
  expect_equal(unname(tbg$z1), unname(ref$grad[, c(1L, 3L)]), tolerance = 1e-6)
  expect_equal(unname(tbg$z2), unname(ref$grad[, c(2L, 4L)]), tolerance = 1e-6)

  # no third formula part: all z variables are inefficiency determinants
  expect_equal(dim(fit$smooth_coefficients_grad$alpha_grad), c(n, 2L))
  expect_equal(colnames(fit$inefficiency_grad), c("z1", "z2"))

  out <- capture.output(print(fit))
  expect_true(any(grepl("local-linear", out)))
  expect_true(any(grepl("GRADIENTS", out)))

  # a z variable that is not an inefficiency determinant
  fit2 <- suppressWarnings(
    spscsfa(y ~ x1 | z1 + z2 | z1, data = dat, spscreg = "ll",
             bandwidth = bw, gradients = TRUE, verbose = FALSE)
  )
  expect_equal(
    unname(fit2$smooth_coefficients_grad$alpha_grad[, 2L]),
    unname(fit2$smooth_coefficients_grad$theta_beta_grad$z2[, 1L]),
    tolerance = 1e-8
  )
  expect_message(g <- grad(fit2, zvar = "z2"), "not an inefficiency determinant")
  expect_null(g$inefficiency_grad)

  # the plot method uses the local-linear estimator
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  curves <- plot(fit, n_grid = 10)
  expect_named(curves, c("z1", "z2"))
  expect_equal(nrow(curves$z1), 10L)
})

test_that("a user-provided bandwidth for the local-linear estimator", {

  skip_on_cran()

  dat <- make_data(80)

  # llscoefbw() with LSCV = FALSE returns a named vector
  bw <- llscoefbw(y ~ x1 | z1 + z2, bws = c(0.3, 0.4), data = dat,
                  LSCV = FALSE)

  fit <- suppressWarnings(
    spscsfa(y ~ x1 | z1 + z2, data = dat, spscreg = "ll", bandwidth = bw,
             verbose = FALSE)
  )

  expect_true(fit$smooth_coefficients$bandwidth_user_supplied)
  expect_equal(unname(fit$smooth_coefficients$bandwidth), c(0.3, 0.4))

  # ... and the same bandwidths without names, in the order of the formula
  fit2 <- suppressWarnings(
    spscsfa(y ~ x1 | z1 + z2, data = dat, spscreg = "ll",
             bandwidth = c(0.3, 0.4), verbose = FALSE)
  )
  expect_equal(fit2$smooth_coefficients$theta_beta,
               fit$smooth_coefficients$theta_beta)
})

test_that("gradients = TRUE when there is no continuous z variable", {

  skip_on_cran()

  dat <- make_data()

  # only a factor in the z part
  msgs <- character(0)
  fit <- withCallingHandlers(
    spscsfa(y ~ x1 + x2 | z3, data = dat, gradients = TRUE, verbose = FALSE),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl("no continuous z variable", msgs)))
  expect_s3_class(fit, "spscsfa")

  # no gradients of the smooth coefficients nor of alpha ...
  scg <- fit$smooth_coefficients_grad
  expect_length(scg$theta_beta_grad, 0L)
  expect_equal(dim(scg$alpha_grad), c(nrow(dat), 0L))

  # ... but the gradients of the inefficiency, for the binary variables
  # converted from the factor
  expect_equal(colnames(fit$inefficiency_grad), c("z3_B", "z3_C"))

  # print and summary: the note about the training observations comes with
  # the inefficiency gradients, and a line says why the others are missing
  s <- summary(fit)
  expect_null(s$gradients$theta_beta_grad)
  expect_null(s$gradients$alpha_grad)
  expect_equal(rownames(s$gradients$inefficiency_grad), c("z3_B", "z3_C"))

  out <- capture.output(print(fit))
  expect_true(any(grepl("No continuous z variable", out)))
  expect_true(any(grepl("JLMS inefficiency: gradients", out)))

  # grad() works too
  expect_message(g <- grad(fit, zvar = "z3"), "not a continuous z variable")
  expect_equal(colnames(g$inefficiency_grad), c("z3_B", "z3_C"))
})
