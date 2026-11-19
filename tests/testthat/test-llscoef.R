test_that("a bandwidth object from llscoefbw() keeps the formula: none is accepted", {

  skip_on_cran()

  set.seed(1)
  n <- 60
  x <- rnorm(n)
  z <- runif(n)
  y <- 1 + sin(2 * pi * z) * x + rnorm(n, 0, 0.1)
  dat <- data.frame(y, x, z)

  bw <- llscoefbw(y ~ x | z, data = dat)
  expect_s3_class(bw, "llsmoothcoefficientbw")

  # the formula is that of llscoefbw()
  fit <- llscoef(bw, data = dat, betas = TRUE)
  expect_s3_class(fit, "llsmoothcoefficient")
  expect_equal(dim(coef(fit)), c(n, 2L))

  # a formula given with the object is an error, by name or by position
  expect_error(llscoef(bw, formula = y ~ x | z, data = dat),
               "must not be given")
  expect_error(llscoef(bw, y ~ x | z, data = dat), "must not be given")

  # with the vector of bandwidths the formula is required, and the estimates
  # are the same
  expect_error(llscoef(bw$bw, data = dat), "formula")
  fit.b <- llscoef(bw$bw, formula = y ~ x | z, data = dat, betas = TRUE)
  expect_equal(coef(fit.b), coef(fit), tolerance = 1e-6)
})

test_that("the local-constant and local-linear bandwidth objects are not mixed", {

  set.seed(1)
  n <- 30
  dat <- data.frame(y = rnorm(n), x = rnorm(n), z = runif(n))

  # an object of npscoefbw() has class "scbandwidth"; one that records no
  # regtype is one of regtype = "lc" (the default of np): its bandwidths are for
  # the local-constant estimator
  bw.lc <- structure(list(bw = 0.3), class = "scbandwidth")
  expect_error(llscoef(bw.lc, formula = y ~ x | z, data = dat), "npscoefbw")
  # For spscsfa(spscreg = "lc"), an omitted regtype and omitted settings
  # represent the np defaults: LC, fixed bandwidths, and default kernels.
  expect_identical(spscsfa:::.lcscoef_np_bw(bw.lc, "bandwidth"), 0.3)
  bw.lc.explicit <- structure(
    list(bw = 0.3, regtype = "lc", type = "fixed", scaling = FALSE,
         ckertype = "gaussian", ckerorder = 2, ckerbound = "none",
         ukertype = "aitchisonaitken", okertype = "liracine"),
    class = "scbandwidth"
  )
  expect_identical(
    spscsfa:::.lcscoef_np_bw(bw.lc.explicit, "bandwidth"), 0.3
  )

  # nor does any other list with a component bw
  expect_error(llscoef(list(bw = 0.3), formula = y ~ x | z, data = dat),
               "numeric vector")

  # llscoefbw() rejects it too, and so does spscsfa() with the local-linear
  # estimator (an object without a regtype is one of regtype = "lc")
  expect_error(llscoefbw(y ~ x | z, bws = bw.lc, data = dat,
                         LSCV = FALSE), "npscoefbw")
  expect_error(spscsfa(y ~ x | z, data = dat, spscreg = "ll",
                       bandwidth = bw.lc, verbose = FALSE), "npscoefbw")

  # npscoefbw() selects local-linear bandwidths with regtype = "ll": the object
  # is rejected by spscsfa() with the local-constant estimator ...
  bw.np.ll <- structure(list(bw = 0.3, regtype = "ll"), class = "scbandwidth")
  expect_error(spscsfa(y ~ x | z, data = dat, spscreg = "lc",
                       bandwidth = bw.np.ll, verbose = FALSE),
               "regtype = \"ll\"")

  # ... but not with the local-linear one
  msg <- tryCatch({
    suppressWarnings(spscsfa(y ~ x | z, data = dat, spscreg = "ll",
                             bandwidth = bw.np.ll, verbose = FALSE))
    ""
  }, error = function(e) conditionMessage(e))
  expect_false(grepl("is an object returned by npscoefbw()", msg, fixed = TRUE))

  # a local-polynomial object is not accepted by either estimator
  bw.np.lp <- structure(list(bw = 0.3, regtype = "lp"), class = "scbandwidth")
  expect_error(spscsfa(y ~ x | z, data = dat, spscreg = "ll",
                       bandwidth = bw.np.lp, verbose = FALSE),
               "regtype = \"lp\"")
  expect_error(spscsfa(y ~ x | z, data = dat, spscreg = "lc",
                       bandwidth = bw.np.lp, verbose = FALSE),
               "regtype = \"lp\"")

  # llscoef() and llscoefbw() use the bandwidths of such an object (regtype =
  # "ll", fixed bandwidths, default kernels) as a vector of bandwidths
  fit.np  <- llscoef(bw.np.ll, formula = y ~ x | z, data = dat, betas = TRUE)
  fit.num <- llscoef(0.3, formula = y ~ x | z, data = dat, betas = TRUE)
  expect_equal(coef(fit.np), coef(fit.num))
  expect_equal(unname(llscoefbw(y ~ x | z, bws = bw.np.ll, data = dat,
                                LSCV = FALSE)), 0.3)

  # ... and reject the local-polynomial one
  expect_error(llscoef(bw.np.lp, formula = y ~ x | z, data = dat),
               "regtype = \"lp\"")
  expect_error(llscoefbw(y ~ x | z, bws = bw.np.lp, data = dat,
                         LSCV = FALSE), "regtype = \"lp\"")

  # A local-constant npscoefbw() object must also use fixed bandwidths
  # and all the default kernels used by spscsfa(spscreg = "lc").
  odd <- list(
    bwtype    = list(type = "generalized_nn"),
    bwscaling = list(scaling = TRUE),
    ckertype  = list(ckertype = "epanechnikov"),
    ckerorder = list(ckerorder = 4),
    ckerbound = list(ckerbound = "range"),
    ukertype  = list(ukertype = "liracine"),
    okertype  = list(okertype = "wangvanryzin")
  )

  for (nm in names(odd)) {
    bw.odd.lc <- structure(c(list(bw = 0.3, regtype = "lc"), odd[[nm]]),
                           class = "scbandwidth")
    expect_error(
      spscsfa(y ~ x | z, data = dat, spscreg = "lc",
              bandwidth = bw.odd.lc, verbose = FALSE),
      nm
    )

    bw.odd <- structure(c(list(bw = 0.3, regtype = "ll"), odd[[nm]]),
                        class = "scbandwidth")
    expect_error(llscoef(bw.odd, formula = y ~ x | z, data = dat), nm)
    expect_error(llscoefbw(y ~ x | z, bws = bw.odd, data = dat,
                           LSCV = FALSE), nm)
    expect_error(spscsfa(y ~ x | z, data = dat, spscreg = "ll",
                         bandwidth = bw.odd, verbose = FALSE), nm)
  }

  # the regtype of the object: "lc" if it records none
  expect_identical(spscsfa:::.spsc_np_regtype(bw.lc), "lc")
  expect_identical(spscsfa:::.spsc_np_regtype(bw.np.ll), "ll")

  # and spscsfa() with the local-constant estimator rejects the bandwidth
  # object of llscoefbw()
  bw.ll <- structure(list(bw = 0.3, formula = y ~ x | z),
                     class = "llsmoothcoefficientbw")
  expect_error(spscsfa(y ~ x | z, data = dat, spscreg = "lc",
                       bandwidth = bw.ll, verbose = FALSE), "llscoefbw")
})
