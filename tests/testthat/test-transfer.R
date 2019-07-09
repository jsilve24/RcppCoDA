test_that("basic truths", {
  V <- ilrContrast(5)
  expect_equal(iiContrast(V,V), diag(4))
  expect_equal(aaContrast(3,3,5), diag(4))
})

test_that("transfers of data correct", {
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

test_that("transfers of covariance are correct", {
  X <- matrix(abs(rnorm(10)), 5, 2)
  X <- clo(X)
  X <- ilr(X)
  V <- ilrContrast(nrow(X)+1)
  Sigma <- cov(t(X))
  
  foo <- ilrvar2clrvar_internal(Sigma, V)
  expect_equal(foo, t(V)%*%Sigma%*%V)
  Sigma <- array(c(Sigma, .0001*Sigma), c(4,4,2))
  Sigma.clr <- ilrvar2clrvar(Sigma, V)
  expect_equal(.0001*foo, Sigma.clr[,,2])
})
