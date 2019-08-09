test_that("Contrast Dimensions", {
  
  # Standard here is that a contrast matrix is defined as d1 x d2 where
  # d1 is dimension to go to and d2 is dimension currently in. 
  # The exception to this is inverse log-ratio transforms where 
  # d1 and d2 are always defined based on the forward not inverse function. 
  
  D <- 10
  d <- 4
  # ALR, CLR, and ILR
  expect_equal(dim(alrContrast(d,D,inv=TRUE)), c(D-1, D))
  expect_equal(dim(alrContrast(d,D,inv=FALSE)), c(D-1, D))
  expect_equal(dim(ilrContrast(D)), c(D-1, D))
  expect_equal(dim(clrContrast(D, inv=FALSE)), c(D, D))
  expect_equal(dim(clrContrast(D, inv=TRUE)), c(D, D))
  
  # IQLR
  expect_equal(dim(var2iqlrContrast(1:D+0.1, qLow=0.25, qHigh=0.75)), c(D,D))
  expect_error(var2iqlrContrast(1:D+0.1, .1, .01)) # qLow must be greater than qHigh
  
  # Transfer Contrasts
  V1 <- ilrContrast(D)
  V2 <- ilrContrast(D)
  V2 <- V2[-1,]
  expect_equal(dim(iiTransfer(V1, V2)), c(nrow(V2), nrow(V1)))
  expect_equal(dim(icTransfer(V1)), c(D, D-1))
  expect_equal(dim(ciTransfer(V1)), c(D-1, D))
  expect_equal(dim(iaTransfer(V1, d, D)), c(D-1, D-1))
  expect_equal(dim(aiTransfer(d, V1, D)), c(D-1, D-1)) ###
  expect_equal(dim(caTransfer(d, D)), c(D-1, D))
  expect_equal(dim(acTransfer(d, D)), c(D, D-1))
  expect_equal(dim(aaTransfer(d, d, D)), c(D-1, D-1))
  
  # Expect certain equal to equality
  expect_equal(iiTransfer(V1, V1), diag(D-1))
  expect_equal(aaTransfer(d,d,D), diag(D-1))
  
  # test ilrContrast with sign matrix
  S <- rbind(c(1,1,-1,-1), c(1,-1, 0, 0), c(0, 0, 1, -1))
  expect_equal(dim(ilrContrast(S)), c(3,4))
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