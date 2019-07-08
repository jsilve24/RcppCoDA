#' If X is a vector - convert it into a matrix (with 1 column)
#' If X is an array it stays an array
#' @param X thing to test/transform
#' @return matrix
#' @export
#'
#' @examples
#' vec_to_array(c(1,2,3))
#' vec_to_array(rbind(c(1,2,3), c(1,2,3)))
vec_to_array <- function(x){
  if (is.vector(x)) {
    n <- names(x)
    x <- matrix(x, ncol = 1)
    colnames(x) <- n
  }
  x
}
  
#' Functions for setting up transforms and handing weird arrays 
#' @param X array
#' @param b part/coord index in array
#' @return matrix (if array_pre) array (if array_post)
#' @return s dimensions of array desired (for array_post)
#' @name array_munge
array_pre <- function(X, b){
  #X <- vec_to_array(X)
  if (b < 1) stop("b must be >=1")
  if (b > length(dim(X))) stop("b must be <= number of array dimensions")
  s <- dim(X)
  if (b !=1){
    perm <- rep(0, length(s))
    seen <- FALSE
    for (i in 1:length(s)){
      if (i==b) { perm[i] <- 1 ; seen <- TRUE}
      else {
        if (!seen) {perm[i] <- i+1}
        else {perm[i] <- i}  
      }
    }
    X <- aperm(X, perm)
  }
  X <- matrix(X, s[b], prod(s[-b]))
  return(X)
}

#' @rdname array_munge
array_post <- function(X, b, s){
  #X <- vec_to_array(X)
  if (b < 1) stop("d must be >=1")
  if (b > length(s)) stop("d must be <= number of array dimensions")
  stmp <- c(s[b], s[-b]) # handles if b==1 or not
  X <- array(X, dim=stmp)
  if (b != 1){
    perm <- 1:length(s)
    perm <- c(b, perm[-b])
    X <- aperm(X, perm)
  }
  return(X)
}

#' Closure operation - divide elements by sum of elements in b dimension
#' @param X vector, matrix, or array of data 
#'   (b is the dimension which is relevant to be closed)
#' @param b index of dimension to operate on 
#'   (e.g., index of dimension of parts or coords in X;
#'   default is 1 meaning that compositions/log-ratios are rows)
clo <- function(X, b=1){
  X <- vec_to_array(X)
  s <- dim(X)
  X <- array_pre(X, b)
  X <- clo_internal(X)
  X <- array_post(X, b, s)
  return(X)
}

#' Center operation - subtract from each element the mean of elements in b dimension
#' @param X vector, matrix, or array of data 
#'   (b is the dimension which is relevant to be centered)
#' @param b index of dimension to operate on 
#'   (e.g., index of dimension of parts or coords in X;
#'   default is 1 meaning that compositions/log-ratios are rows)
center <- function(X, b=1){
  X <- vec_to_array(X)
  s <- dim(X)
  X <- array_pre(X, b)
  X <- center_internal(X)
  X <- array_post(X, b, s)
  return(X)
}


#' Log-Ratio Transformation
#'
#' \code{glr} is generic log-ratio transform, code used by other
#' transforms, can be called directly. *Contrast functions produce contrast 
#' matricies (e.g., "V") that define the various transforms. 
#'
#' @param X vector, matrix, or array of data (b is the dimension which is compositionally relevant)
#' @param V transformation matrix (defines transform; P x D where D is number of parts)
#'   if NULL then uses ilr default basis (see details)
#' @param d for ALR, which component (integer position) to take as reference
#' (default is nrow(x) for alr or nrow(x)+1 for alrInv) for alrInv corresponds 
#' to column position in untransformed matrix.
#' @param b index of dimension to operate on (e.g., index of dimension of parts or coords in X;
#' default is 1 meaning that compositions/log-ratios are rows)
#' @param inv for ALR and CLR, transformation matrix is different forward and inverse
#' @param D the number of parts (e.g., number of columns in untransformed data)
#' @return matrix (converts vectors to column matricies)
#' @details The default ILR base formed by Gram-Schmidt orthogonalization of an ALR_D basis.
#' @name base_lr_transforms
#' @examples
#' #ALR Transform
#' x <- matrix(runif(30), 10, 3)
#' x <- clo(x)
#' x.alr <- alr(x, 2)
#' x <- alrInv(x.alr, 2)
#'
#' # ILR
#' x.ilr <- ilr(x)
#' x <- ilrInv(x.ilr)
#'
#' # CLR
#' x.clr <- clr(x)
#' x <- clrInv(x.clr)
#'
#' # CUSTOM - Be careful if your custom matrix is not
#' # orthogonal the inverse transform may not be given by just the transpose!
#' # For example, this is the case for the ALR
#' V <- matrix(c(1, 1, -1), 1, 3)
#' x.custom <- glr(x, V)
NULL

#' @rdname base_lr_transforms
glr <- function(X, V, b=1){
  X <- vec_to_array(X)
  s <- dim(X)
  if (ncol(V) != s[b]) stop("V has wrong dimensions")
  X <- array_pre(X, b)
  X <- glr_internal(X, V)
  s[b] <- nrow(V)
  X <- array_post(X, b, s)
  return(X)
}

#' @rdname base_lr_transforms
glrInv <- function(X, V, b=1){
  X <- vec_to_array(X)
  s <- dim(X)
  if (nrow(V) != s[b]) stop("V has wrong dimensions")
  X <- array_pre(X, b)
  X <- glrInv_internal(X, V)
  s[b] <- ncol(V)
  X <- array_post(X, b, s)
  return(X)
}

#' @rdname base_lr_transforms
alr <- function(X, d=NULL, b=1){
  X <- vec_to_array(X)
  if (is.null(d)) d <- nrow(X)
  s <- dim(X)
  if (d > s[b] |d <1) stop("d must be between 1 and dim(X)[b]+1")
  X <- array_pre(X, b)
  X <- alr_internal(X, d)
  s[b] <- s[b]-1
  X <- array_post(X, b, s)
  return(X)
}

#' @rdname base_lr_transforms
alrInv <- function(X, d=NULL, b=1){
  X <- vec_to_array(X)
  if (is.null(d)) d <- nrow(X)+1
  s <- dim(X)
  if (d > s[b]+1 | d < 1) stop("d must be between 1 and dim(X)[b]+1")
  X <- array_pre(X, b)
  X <- alrInv_internal(X, d)
  s[b] <- s[b]+1
  X <- array_post(X, b, s)
  return(X)
}

#' @rdname base_lr_transforms
clr <- function(X, b=1){
  X <- vec_to_array(X)
  s <- dim(X)
  X <- array_pre(X, b)
  X <- clr_internal(X)
  X <- array_post(X, b, s)
  return(X)
}

#' @rdname base_lr_transforms
clrInv <- function(X, b=1){
  X <- vec_to_array(X)
  s <- dim(X)
  X <- array_pre(X, b)
  X <- clrInv_internal(X)
  X <- array_post(X, b, s)
  return(X)
}

#' @rdname base_lr_transforms
ilr <- function(X, V=NULL, b=1){
  X <- vec_to_array(X)
  s <- dim(X)
  if (is.null(V)) V <- ilrContrast(s[b])
  if (ncol(V) != s[b]) stop("V has wrong dimensions")
  X <- array_pre(X, b)
  X <- ilr_internal(X, V)
  s[b] <- nrow(V)
  X <- array_post(X, b, s)
  return(X)
}

#' @rdname base_lr_transforms
ilrInv <- function(X, V=NULL, b=1){
  X <- vec_to_array(X)
  s <- dim(X)
  if (is.null(V)) V <- ilrContrast(s[b]+1)
  if (nrow(V) != s[b]) stop("V has wrong dimensions")
  X <- array_pre(X, b)
  X <- ilrInv_internal(X, V)
  s[b] <- ncol(V)
  X <- array_post(X, b, s)
  return(X)
}

