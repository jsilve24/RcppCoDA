## #' @rdname base_lr_transforms
## #' @export
## alr <- function(X, d=NULL) {
##   if(is.null(d)){
##     d <- nrow(X)
##   }
##   alr_internal(X, d)
##  }


## #' @rdname base_lr_transforms
## #' @export
## alrInv <- function(X, d=NULL) {
##   if(is.null(d)){
##     d <- nrow(X)+1
##   }
##   alrInv_internal(X, d)
##  }
