test_that("basic truths", {
  V <- ilrContrast(5)
  expect_equal(iiContrast(V,V), diag(4))
  expect_equal(aaContrast(3,3,5), diag(4))
})

test_that("transfers correct", {
  X <- matrix(abs(rnorm(10)), 5, 2)
  X <- clo(X)
  X.ilr <- ilr(X)
  X.clr <- clr(X)
  X.alr <- alr(X)
  V <- ilrContrast(5)
  expect_equal(X.clr, icContrast(V)%*%X.ilr)
  expect_equal(X.ilr, ciContrast(V)%*%X.clr)
  
  
  expect_equal(X.clr, ilr2clr(X.ilr, V))
  expect_equal(X.ilr, clr2ilr(X.clr, V))
  expect_equal(X.ilr, alr2ilr(X.alr, 5, V))
  expect_equal(X.alr, ilr2alr(X.ilr, V, 5))
})
