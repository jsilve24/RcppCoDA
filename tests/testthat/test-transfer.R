test_that("basic truths", {
  V <- ilrContrast(5)
  expect_equal(iiContrast(V,V), diag(4))
  expect_equal(aaContrast(3,3,5), diag(4))
})

