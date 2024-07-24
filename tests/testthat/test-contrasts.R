x <- matrix(runif(300), 100, 3)
x <- clo(x)

# Standard here is that a contrast matrix is defined as d1 x d2 where
# d1 is dimension to go to and d2 is dimension currently in. 
# The exception to this is inverse log-ratio transforms where 
# d1 and d2 are always defined based on the forward not inverse function. 


test_that("Contrast Matricies Have Correct Dimension", {
  D <- 10
  d <- 4
  # ALR, CLR, and ILR
  expect_equal(dim(alrContrast(d,D,inv=TRUE)), c(D-1, D))
  expect_equal(dim(alrContrast(d,D,inv=FALSE)), c(D-1, D))
  expect_equal(dim(ilrContrast(D)), c(D-1, D))
  expect_equal(dim(clrContrast(D, inv=FALSE)), c(D, D))
  expect_equal(dim(clrContrast(D, inv=TRUE)), c(D, D))
})


test_that("ilr produces orthogonal transform",{
  V <- ilrContrast(4)
  expect_equal(V%*%t(V), diag(3))
  expect_equal(t(V)%*%V, diag(4) - 1/4)
  
  # Expect error when not a valid sign matrix
  S <- rbind(c(1,1,1,-1), c(1,-1, 0, 0), c(0, 0, 1, -1))
  expect_error(ilrContrast(S))
  
  S <- rbind(c(1,1,-1,-1), c(1,-1, 0, 0), c(0, 0, 1, -1))
  V <- ilrContrast(S)
  expect_equal(V%*%t(V), diag(3))
  expect_equal(t(V)%*%V, diag(4) - 1/4)
})



x <- matrix(runif(300), 100, 3)
x <- clo(x)

test_that("alrContrast correct", {
  B.inv.true <- cbind(diag(3), 0)
  B.true <- cbind(diag(3), -1)
  expect_equal(B.inv.true, alrContrast(4,4, TRUE))
  expect_equal(B.true, alrContrast(4,4,FALSE))
  expect_error(alrContrast(5, 4, inv=TRUE))
  
  B.true <- matrix(0, 4, 5)
  foo <- cbind(1:4, c(1,2,4,5))
  B.true[foo] <- 1
  B.true[,3] <- -1
  expect_equal(B.true, alrContrast(3, 5, FALSE))
}) 


test_that("clrContrast correct", {
  B <- matrix(-1, 4, 4) + 4*diag(4)
  B <- B/4
  expect_equal(clrContrast(4, FALSE), B)
  expect_equal(clrContrast(4, TRUE), diag(4))
})
