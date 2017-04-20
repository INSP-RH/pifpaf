#' @title Check covariance of exposure values are well defined.
#'   
#' @description Function that verifies \code{Xvar} are well defined for
#'   \code{paf.confidence} and \code{pif.confidence}
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
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'   
#' @seealso \code{\link{check.confidence}}, \code{\link{check.thetas}}, 
#'   \code{\link{check.cft}}, \code{\link{check.rr}}, 
#'   \code{\link{check.exposure}}, \code{\link{check.integrals}}
#' 
#' @keywords internal
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
