#' @title Check covariance of exposure values was defined, if not a covariance matrix with entries equal to zero is assumed. 
#' 
#' @description Function that verifies Xvar was defined
#' 
#' @param Xvar     Input of covariance matrix of exposure values
#' 
#' @return Xvar
#' 
#' @importFrom matrixcalc is.negative.definite
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
  
  #Convert to matrix
  Xvar <- as.matrix(Xvar)
  
  #Check nan, NA and character
  if(any(is.na(Xvar)) || any(is.nan(Xvar)) || any(is.character(Xvar))){
    
    warning("Covariance of exposure values Xvar had non-numeric arguments, defaulting everything to zero")
    n    <- ncol(Xvar)
    Xvar <- matrix(0, ncol = n, nrow = n)
    
  }
  
  #Check is squared
  if(ncol(Xvar) != nrow(Xvar)){
    stop("Covariance matrix is not square matrix")
  }
  
  #Check is positive definite
  if(is.negative.definite(.Xvar)){
    stop("Covariance matrix must be positive semi-definite.")
  }
  
  return(Xvar)
}
