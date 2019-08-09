#' @rdname base_lr_transforms
#' @export
ilrContrast <- function(d){
  if (length(d)==1){
    return(ilrContrastDefault_internal(as.integer(d)))
  } else if (length(d) > 1){
    if (!all(rowSums(d) ==0)) stop("rowSums of your contrast matrix do not sum to 0")
    mode(d) <- "integer"
    return (ilrContrastSign_internal(d))
  }
  else{
    stop("not sure what you did but d must either be an integer or a sign matrix")
  }
}
