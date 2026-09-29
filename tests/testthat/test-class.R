test_that("extractors reject wrong class", {
  expect_error(te(list()), "class 'spsc_sfa'")
  expect_error(inefficiency(list()), "class 'spsc_sfa'")
  expect_error(smooth_coef(list()), "class 'spsc_sfa'")
})

test_that("public functions exist", {
  expect_true(is.function(spsc_sfa))
  expect_true(is.function(te))
  expect_true(is.function(inefficiency))
  expect_true(is.function(smooth_coef))
})
