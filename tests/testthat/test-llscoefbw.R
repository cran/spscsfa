# Data of an example in which a single Nelder-Mead search that starts at half
# of the rule-of-thumb bandwidth stops in a local minimum of the LSCV criterion
# (near h = 0.05), while the criterion is much lower near h = 0.8.
max_search_data <- function() {
  set.seed(1)
  n <- 100
  x <- rnorm(n)
  z <- runif(n)
  rnorm(n)  # a variable that the example draws here: keeps its random numbers
  y <- cos(z) + sin(z) * x + rnorm(n, 0, 0.1) -
    exp(0.1 + 0.1 * z) * abs(rnorm(n, 0, 1))
  data.frame(y = y, x = x, z = z)
}

test_that("max.search must be a positive integer and is ignored without a search", {

  dat <- max_search_data()

  for (bad in list(0, -1, 1.5, NA, Inf, "2", c(1, 2), NULL)) {
    expect_error(llscoefbw(y ~ x | z, data = dat, max.search = bad), "max.search")
  }

  bw <- llscoefbw(y ~ x | z, bws = 0.2, data = dat, LSCV = FALSE,
                  max.search = 3)
  expect_equal(unname(bw), 0.2)
  expect_equal(names(bw), "z")
})

test_that("the scan rescues a search that starts in a local minimum", {

  skip_on_cran()

  dat <- max_search_data()
  n <- nrow(dat)
  start <- 0.5 * 1.059224 * sd(dat$z) * n^(-1 / 5)

  bw.1 <- llscoefbw(y ~ x | z, bws = start, data = dat, max.search = 1)
  grid_scan <- 2^c(-2:5)
  bw.2 <- llscoefbw(y ~ x | z, bws = start, data = dat, max.search = 2,
                    grid_scan = grid_scan)

  # never worse, and here clearly better
  expect_lte(bw.2$lscv, bw.1$lscv)
  expect_lt(bw.2$lscv, 0.9 * bw.1$lscv)
  expect_gt(bw.2$bw[["z"]], bw.1$bw[["z"]])

  # the evaluations of the scan and of the further search are counted
  expect_gt(bw.2$n_feval, bw.1$n_feval + length(grid_scan))

  # max.search = 2 is the default
  bw.d <- llscoefbw(y ~ x | z, bws = start, data = dat)
  expect_equal(bw.d$bw, bw.2$bw)
  expect_equal(bw.d$lscv, bw.2$lscv)
})

test_that("from the rule-of-thumb start the result is never worse", {

  skip_on_cran()

  dat <- max_search_data()

  bw.1 <- llscoefbw(y ~ x | z, data = dat, max.search = 1)
  bw.3 <- llscoefbw(y ~ x | z, data = dat, max.search = 3)

  expect_lte(bw.3$lscv, bw.1$lscv)

  # with max.search > 1 the scan evaluations are made and counted
  grid_scan <- eval(formals(spscsfa:::llscoefbw.formula)$grid_scan)
  expect_gte(bw.3$n_feval, bw.1$n_feval + length(grid_scan))
})

test_that("there is no scan without a continuous z variable", {

  skip_on_cran()

  set.seed(2)
  n <- 80
  x <- rnorm(n)
  f <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  y <- 1 + ifelse(f == "A", 1, -1) * x + rnorm(n, 0, 0.3)
  dat <- data.frame(y, x, f)

  bw.1 <- llscoefbw(y ~ x | f, data = dat, max.search = 1)
  bw.3 <- llscoefbw(y ~ x | f, data = dat, max.search = 3)

  expect_equal(bw.3$bw, bw.1$bw)
  expect_equal(bw.3$n_feval, bw.1$n_feval)
})

test_that("with continuous and categorical z variables the result is never worse", {

  skip_on_cran()

  set.seed(3)
  n <- 80
  x <- rnorm(n)
  z <- runif(n)
  f <- factor(sample(c("A", "B"), n, replace = TRUE))
  y <- sin(2 * pi * z) * x + (f == "B") + rnorm(n, 0, 0.2)
  dat <- data.frame(y, x, z, f)

  bw.1 <- llscoefbw(y ~ x | z + f, data = dat, max.search = 1)
  bw.3 <- llscoefbw(y ~ x | z + f, data = dat, max.search = 3)

  expect_lte(bw.3$lscv, bw.1$lscv)
  expect_named(bw.3$bw, c("z", "f"))
})
