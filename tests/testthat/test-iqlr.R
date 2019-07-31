test_that("var2iqlrContrast is correct", {
  V <- var2iqlrContrast(matrix(1:10)+.1, .4, .6)
  V.true <- diag(10)
  V.true[,4:7] <- V.true[,4:7] -.25
  expect_equal(V, V.true)
})


test_that("clrvar2iqlrvar is correct", {
  X <- matrix(abs(rnorm(10)), 5, 2)
  X <- clo(X)
  X <- ilr(X)
  V <- ilrContrast(nrow(X)+1)
  Sigma <- cov(t(X))
  Sigma <- array(c(Sigma, .0001*Sigma), c(4,4,2))
  Sigma <- ilrvar2clrvar(Sigma, V)
  Sigma.iqlr <- clrvar2iqlrvar(Sigma)
  Viqlr1 <- var2iqlrContrast(diag(Sigma[,,1]), .25, .75)
  Viqlr2 <- var2iqlrContrast(diag(Sigma[,,2]), .25, .75)
  expect_equal(Sigma.iqlr[,,1], clrvar2ilrvar(Sigma[,,1], Viqlr1))
  expect_equal(Sigma.iqlr[,,2], clrvar2ilrvar(Sigma[,,2], Viqlr2))
})
