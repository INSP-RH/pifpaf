#' @title Check covariance of exposure values was defined.
#' 
#' @description Function that verifies Xvar was defined 
#' 
#' @param Xvar     Input of covariance matrix of exposure values
#' 
#' @return Xvar
#' 
#' @importFrom matrixcalc is.positive.semi.definite
#' 
#' @examples 
#' 
#' #Example 1
#' Xvar <- 0.2
#' check.xvar(Xvar)
#' 
#' @export

check.xvar <- function(Xvar){
  
  #Convert to matrix
  Xvar <- as.matrix(Xvar)
  
  #Check nan, NA and character
  if(any(is.na(Xvar)) || any(is.nan(Xvar)) || any(is.character(Xvar))){
    
    stop("Covariance of exposure values Xvar had non-numeric arguments")
    
  }
  
  #Check is squared
  if(ncol(Xvar) != nrow(Xvar)){
    stop("Covariance matrix is not square matrix")
  }
  
  #Check is positive definite
  if(!is.positive.semi.definite(Xvar)){
    stop("Covariance matrix must be positive semi-definite.")
  }
  
  return(Xvar)
}
