#' Transfer data and covariance matricies from one coordinate system to another
#' 
#' Naming Conventions:
#' xxTransfer produces a transfer matrix for transforming between representations where 
#' i=ILR, c=CLR, a=ALR. x2x represent transformations of data between coordinate systems 
#' xvar2xvar represent transformations of covariance matricies between coordinate systems. 
#' By convention all functions take arguments as follows: all functions have the form function(a,b,c)
#' where a is the data/covariance, b is the argument for specifying the current coordinate system, and c 
#' represents the argument for specifying the desired coordinate system. 
#' @param X data to transform between representations
#' @param b index of dimension to operate on (e.g., index of dimension of parts or coords in X;
#' default is 1 meaning that compositions/log-ratios are rows) Note, for covariance matricies, 
#' b is meaningless, first two dimensions must be log-ratio coordinates. 
#' @param Sigma covariance matrix in specified transformed space
#' @param D total number of parts in Sigma 
#' @param V1 ILR contrast matrix of basis Sigma is already in
#' @param V2 ILR contrast matrix of basis Sigma is desired in
#' @param d1 alr reference element Sigma is already expressed with respec to
#' @param d2 alr reference element Sigma is to be expressed with respect to
#' @return matrix
#' @name convert_coda
#' @examples
#' # For Data
#' X <- matrix(abs(rnorm(10)), 5, 2)
#' X <- clo(X)
#' X.clr <- clr(X)
#' X.alr <- clr2alr(X, 3)
#' 
#' # For a covariance matrix starting in the alr
#' Sigma.alr <- X.alr%*%t(X.alr) # Just a random covariance matrix 
#' Sigma.ilr <- alrvar2ilrvar(Sigma.alr, 3, ilrContrast(5)) # Covert to default ILR
#' 
#' # Another way of doing this more "manually" - use contrast matricies
#' V <- caContrast(3, nrow(X)) # create contrast from clr to alr
#' X.alr <- V %*% X.clr
NULL

#' @rdname convert_coda
#' @export
ilr2ilr <- function(X, V1, V2, b=1){
  X <- vec_to_array(X)
  s <- dim(X)
  if (nrow(V1) != s[b]) stop("V1 has wrong dimensions")
  if (nrow(V2) != s[b]) stop("V2 has wrong dimensions")
  X <- array_pre(X, b)
  X <- ilr2ilr_internal(X, V1, V2)
  s[b] <- nrow(V2)
  X <- array_post(X, b, s)
  return(X)
}

#' @rdname convert_coda
#' @export
ilr2clr <- function(X, V1, b=1){
  X <- vec_to_array(X)
  s <- dim(X)
  if (nrow(V1) != s[b]) stop("V1 has wrong dimensions")
  X <- array_pre(X, b)
  X <- ilr2clr_internal(X, V1)
  s[b] <- ncol(V1)
  X <- array_post(X, b, s)
  return(X)
}

#' @rdname convert_coda
#' @export
clr2ilr <- function(X, V2, b=1){
  X <- vec_to_array(X)
  s <- dim(X)
  if (ncol(V2) != s[b]) stop("V1 has wrong dimensions")
  X <- array_pre(X, b)
  X <- clr2ilr_internal(X, V2)
  s[b] <- nrow(V2)
  X <- array_post(X, b, s)
  return(X)
}

#' @rdname convert_coda
#' @export
alr2clr <- function(X, d1, b=1){
  X <- vec_to_array(X)
  s <- dim(X)
  #if (d1 != s[b]) stop("V1 has wrong dimensions")
  X <- array_pre(X, b)
  X <- alr2clr_internal(X, d1)
  s[b] <- nrow(X)
  X <- array_post(X, b, s)
  return(X)
}

#' @rdname convert_coda
#' @export
clr2alr <- function(X, d2, b=1){
  X <- vec_to_array(X)
  s <- dim(X)
  #if (d1 != s[b]) stop("V1 has wrong dimensions")
  X <- array_pre(X, b)
  X <- clr2alr_internal(X, d2)
  s[b] <- nrow(X)
  X <- array_post(X, b, s)
  return(X)
}

#' @rdname convert_coda
#' @export
alr2alr <- function(X, d1, d2, b=1){
  X <- vec_to_array(X)
  s <- dim(X)
  #if (d1 != s[b]) stop("V1 has wrong dimensions")
  X <- array_pre(X, b)
  X <- alr2alr_internal(X, d1, d2)
  #s[b] <- nrow(X)
  X <- array_post(X, b, s)
  return(X)
}

#' @rdname convert_coda
#' @export
ilr2alr <- function(X, V1, d2, b=1){
  X <- vec_to_array(X)
  s <- dim(X)
  if (nrow(V1) != s[b]) stop("V1 has wrong dimensions")
  X <- array_pre(X, b)
  X <- ilr2alr_internal(X, V1, d2)
  s[b] <- nrow(X)
  X <- array_post(X, b, s)
  return(X)
}

#' @rdname convert_coda
#' @export
alr2ilr <- function(X, d1, V2, b=1){
  X <- vec_to_array(X)
  s <- dim(X)
  if (nrow(V2) != s[b]) stop("V1 has wrong dimensions")
  X <- array_pre(X, b)
  X <- alr2ilr_internal(X, d1, V2)
  s[b] <- nrow(X)
  X <- array_post(X, b, s)
  return(X)
}


# Covariance Functions ----------------------------------------------------


#' @rdname convert_coda
#' @export
ilrvar2ilrvar <- function(Sigma, V1, V2){
  b <- 1:2
  Sigma <- vec_to_array(Sigma)
  s <- dim(Sigma)
  if (any(nrow(V1) != s[b])) stop("V1 has wrong dimensions")
  if (ncol(V1) != ncol(V2)) stop("V2 has wrong dimensions")
  Sigma <- array_pre(Sigma, b)
  Sigma <- ilrvar2ilrvar_internal(Sigma, V1, V2)
  s[b] <- nrow(V2)
  Sigma <- array_post(Sigma, b, s)
  return(Sigma)
}

#' @rdname convert_coda
#' @export
glrvar2glrvar <- function(Sigma, V1, V2){
  ilrvar2ilrvar(Sigma, V1, V2)
}

#' @rdname convert_coda
#' @export
ilrvar2clrvar <- function(Sigma, V1){
  b <- 1:2
  Sigma <- vec_to_array(Sigma)
  s <- dim(Sigma)
  if (any(nrow(V1) != s[b])) stop("V1 has wrong dimensions")
  Sigma <- array_pre(Sigma, b)
  Sigma <- ilrvar2clrvar_internal(Sigma, V1)
  s[b] <- ncol(V1)
  Sigma <- array_post(Sigma, b, s)
  return(Sigma)
}

#' @rdname convert_coda
#' @export
clrvar2ilrvar <- function(Sigma, V2){
  b <- 1:2
  Sigma <- vec_to_array(Sigma)
  s <- dim(Sigma)
  if (any(ncol(V2) != s[b])) stop("V1 has wrong dimensions")
  Sigma <- array_pre(Sigma, b)
  Sigma <- clrvar2ilrvar_internal(Sigma, V2)
  s[b] <- nrow(V2)
  Sigma <- array_post(Sigma, b, s)
  return(Sigma)
}

#' @rdname convert_coda
#' @export
alrvar2clrvar <- function(Sigma, d1){
  b <- 1:2
  Sigma <- vec_to_array(Sigma)
  s <- dim(Sigma)
  #if (d1 != s[b]) stop("V1 has wrong dimensions")
  Sigma <- array_pre(Sigma, b)
  Sigma <- alrvar2clrvar_internal(Sigma, d1)
  s[b] <- nrow(Sigma)
  Sigma <- array_post(Sigma, b, s)
  return(Sigma)
}

#' @rdname convert_coda
#' @export
clrvar2alrvar <- function(Sigma, d2){
  b <- 1:2
  Sigma <- vec_to_array(Sigma)
  s <- dim(Sigma)
  #if (d1 != s[b]) stop("V1 has wrong dimensions")
  Sigma <- array_pre(Sigma, b)
  Sigma <- clrvar2alrvar_internal(Sigma, d2)
  s[b] <- nrow(Sigma)
  Sigma <- array_post(Sigma, b, s)
  return(Sigma)
}

#' @rdname convert_coda
#' @export
alrvar2alrvar <- function(Sigma, d1, d2){
  b <- 1:2
  Sigma <- vec_to_array(Sigma)
  s <- dim(Sigma)
  #if (d1 != s[b]) stop("V1 has wrong dimensions")
  Sigma <- array_pre(Sigma, b)
  Sigma <- alrvar2alrvar_internal(Sigma, d1, d2)
  #s[b] <- nrow(Sigma)
  Sigma <- array_post(Sigma, b, s)
  return(Sigma)
}

#' @rdname convert_coda
#' @export
ilrvar2alrvar <- function(Sigma, V1, d2){
  b <- 1:2
  Sigma <- vec_to_array(Sigma)
  s <- dim(Sigma)
  if (any(nrow(V1) != s[b])) stop("V1 has wrong dimensions")
  Sigma <- array_pre(Sigma, b)
  Sigma <- ilrvar2alrvar_internal(Sigma, V1, d2)
  s[b] <- nrow(Sigma)
  Sigma <- array_post(Sigma, b, s)
  return(Sigma)
}

#' @rdname convert_coda
#' @export
alrvar2ilrvar <- function(Sigma, d1, V2){
  b <- 1:2
  Sigma <- vec_to_array(Sigma)
  s <- dim(Sigma)
  if (any(nrow(V2) != s[b])) stop("V1 has wrong dimensions")
  Sigma <- array_pre(Sigma, b)
  Sigma <- alrvar2ilrvar_internal(Sigma, d1, V2)
  s[b] <- nrow(Sigma)
  Sigma <- array_post(Sigma, b, s)
  return(Sigma)
}


#' Calculate Variation Array from CLR Covariances
#' 
#' Assumes parts are the first two dimensions of Sigma
#' 
#' @param Sigma Covariance matrix Px(PN) where N is number of 
#'   covariance matricies in CLR space
#' @return Array (same dimension as Sigma) but elements represent variation of pariwise log-ratios
#' @export
#' @references Lovell, David, Vera Pawlowsky-Glahn, Juan Jose Egozcue, 
#'   Samuel Marguerat, and Jurg Bahler. 2015. Proportionality: A Valid Alternative 
#'   to Correlation for Relative Data. PLoS Computational Biology 11 (3). 
#'   Public Library of Science: e1004075.
clrvar2vararray <- function(Sigma){
  b <- 1:2
  Sigma <- vec_to_array(Sigma)
  s <- dim(Sigma)
  Sigma <- array_pre(Sigma, b)
  Sigma <- clrvar2vararray_internal(Sigma)
  Sigma <- array_post(Sigma, b, s)
  return(Sigma)
}


#' Calculate Phi Statistics (Proportionality) from CLR Covariances
#' 
#' Assumes parts are the first two dimensions of Sigma
#' 
#' @param Sigma Covariance matrix Px(PN) where N is number of 
#'   covariance matricies in CLR space
#' @return Array (same dimension as Sigma) but elements represent phi statistics
#' @export
#' @references Lovell, David, Vera Pawlowsky-Glahn, Juan Jose Egozcue, 
#'   Samuel Marguerat, and Jurg Bahler. 2015. Proportionality: A Valid Alternative 
#'   to Correlation for Relative Data. PLoS Computational Biology 11 (3). 
#'   Public Library of Science: e1004075.
clrvar2phi<- function(Sigma){
  b <- 1:2
  Sigma <- vec_to_array(Sigma)
  s <- dim(Sigma)
  Sigma <- array_pre(Sigma, b)
  Sigma <- clrvar2phi_internal(Sigma)
  Sigma <- array_post(Sigma, b, s)
  return(Sigma)
}


#' @export
#' @rdname iqlr
clrvar2iqlrvar <- function(Sigma, qLow=.25, qHigh=.75){
  b <- 1:2
  Sigma <- vec_to_array(Sigma)
  s <- dim(Sigma)
  Sigma <- array_pre(Sigma, b)
  Sigma <- clrvar2iqlrvar_internal(Sigma, qLow, qHigh)
  Sigma <- array_post(Sigma, b, s)
  return(Sigma)
}