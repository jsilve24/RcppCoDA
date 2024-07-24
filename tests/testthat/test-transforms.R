test_that("alr and alrInv are inverses and correct", {
  X <- abs(matrix(rnorm(10), 5, 2))
  X <- clo(X)
  Y <- alr(X, 3)
  Y.manual <- X[-3,]
  Y.manual <- sweep(Y.manual, 2, X[3,], FUN="/")
  Y.manual <- log(Y.manual)
  # Correctness
  expect_equal(Y, Y.manual)

  # Default ALR is D
  expect_equal(alr(X, d=nrow(X)), alr(X))
  Y <- alr(X)
  expect_equal(X, alrInv(Y))
  
  # Error when d is too large
  expect_error(alr(X, 100))
  
  # Inverse - and alrInv by association
  expect_equal(X, alrInv(Y))
})


test_that("clr and clrInv are inverses and correct", {
  X <- abs(matrix(rnorm(10), 5, 2))
  X <- clo(X)
  Y <- clr(X)
  Y.manual <- glr(X, clrContrast(5, FALSE))
  # Correctness
  expect_equal(Y, Y.manual)
  
  # Inverse - and alrInv by association
  expect_equal(X, clrInv(Y))
})

test_that("ilr and ilrInv are inverses and correct", {
  X <- abs(matrix(rnorm(10), 5, 2))
  X <- clo(X)
  Y <- ilr(X)
  Y.manual <- ilr(X, ilrContrast(5))
  # Correctness
  expect_equal(Y, Y.manual)
  
  # Inverse - and alrInv by association
  expect_equal(X, ilrInv(Y))
  
  V <- ilrContrast(5)
  expect_error(ilr(X, V[,-1]))
})
