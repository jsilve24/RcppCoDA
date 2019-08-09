test_that("transfers of data correct", {
  X <- matrix(abs(rnorm(10)), 5, 2)
  X <- clo(X)
  X.ilr <- ilr(X)
  X.clr <- clr(X)
  X.alr <- alr(X)
  V <- ilrContrast(5)
  expect_equal(X.clr, icTransfer(V)%*%X.ilr)
  expect_equal(X.ilr, ciTransfer(V)%*%X.clr)
  
  
  expect_equal(X.clr, ilr2clr(X.ilr, V))
  expect_equal(X.ilr, clr2ilr(X.clr, V))
  expect_equal(X.ilr, ilr2ilr(X.ilr, V, V))
  expect_equal(X.ilr, alr2ilr(X.alr, 5, V))
  expect_equal(X.alr, ilr2alr(X.ilr, V, 5))
  expect_equal(X.clr, alr2clr(X.alr, 5))
  expect_equal(X.alr, clr2alr(X.clr, 5))
  expect_equal(alr2alr(X.alr, 5, 3), ilr2alr(X.ilr, V, 3))

})

X <- matrix(abs(rnorm(10)), 5, 2)
X <- clo(X)
X <- ilr(X)
V <- ilrContrast(nrow(X)+1)
Sigma <- cov(t(X))

test_that("transfers of covariance are correct", {
  foo <- ilrvar2clrvar_internal(Sigma, V)
  expect_equal(foo, t(V)%*%Sigma%*%V)
  Sigma.ilr <- array(c(Sigma, .0001*Sigma), c(4,4,2))
  Sigma.clr <- ilrvar2clrvar(Sigma.ilr, V)
  expect_equal(.0001*foo, Sigma.clr[,,2])
  V.alr <- alrContrast(4, 5, FALSE)
  Sigma.alr <- Sigma.ilr; Sigma.alr[] <- 0
  Sigma.alr[,,1] <- V.alr %*% Sigma.clr[,,1] %*% t(V.alr)
  Sigma.alr[,,2] <- V.alr %*% Sigma.clr[,,2] %*% t(V.alr)
  expect_equal(Sigma.alr, clrvar2alrvar(Sigma.clr, 4))
  expect_equal(Sigma.alr, ilrvar2alrvar(Sigma.ilr, V, 4))
  expect_equal(Sigma.ilr, ilrvar2ilrvar(Sigma.ilr, V, V))
  expect_equal(Sigma.clr, alrvar2clrvar(Sigma.alr, 4))
  V.alr2 <- alrContrast(3, 5, FALSE)
  Sigma.alr2 <- Sigma.ilr; Sigma.alr2[] <- 0
  Sigma.alr2[,,1] <- V.alr2 %*% Sigma.clr[,,1] %*% t(V.alr2)
  Sigma.alr2[,,2] <- V.alr2 %*% Sigma.clr[,,2] %*% t(V.alr2)
  expect_equal(Sigma.alr2, alrvar2alrvar(Sigma.alr, 4, 3))
  expect_equal(Sigma.ilr, alrvar2ilrvar(Sigma.alr, 4, V))
})

test_that("phi statistics correct", {
  Sigma <- array(c(Sigma, .0001*Sigma), c(4,4,2))
  foo <- clrvar2phi(Sigma)
  expect_equal(foo[4,3,2], (Sigma[4,4,2] + Sigma[3,3,2] - 2*Sigma[4,3,2]) / Sigma[4,4,2])
})

test_that("variation array correct", {
  Sigma <- array(c(Sigma, .0001*Sigma), c(4,4,2))
  foo <- clrvar2vararray(Sigma)
  expect_equal(foo[4,3,2], (Sigma[4,4,2] + Sigma[3,3,2] - 2*Sigma[4,3,2]))
})

