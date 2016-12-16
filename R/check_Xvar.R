#' @title Check covariance of exposure values was defined, if not a covariance matrix with entries equal to zero is assumed. 
#' 
#' @description Function that verifies Xvar was defined
#' 
#' @param Xvar     Input of covariance matrix of exposure values
#' 
#' @return Xvar
#' 
#' @examples 
#' #Example 1 
#' Xvar <- NA
#' check.xvar(Xvar)
#' 
#' #Example 2
#' Xvar <- 0.2
#' check.xvar(Xvar)
#' 
#' @export

check.xvar <- function(Xvar){
  if(is.na(Xvar)){
    warning("Covariance of exposure values Xvar wasn't defined, default to zero")
    n    <- dim(as.matrix(Xvar))[2]
    Xvar <- matrix(0, ncol = n, nrow = n)
  }
  return(Xvar)
}
