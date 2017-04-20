#' @title Check Relative Risk
#'   
#' @description Function for checking that Relative Risk \code{rr} equals
#'   \code{1} when evaluated in \code{0}.
#'   
#' @param rr        Function for Relative Risk which uses parameter 
#'   \code{theta}. The order of the parameters should be \code{rr(X, theta)}.
#'   
#' @param X         Random sample (data.frame) which includes exposure and
#'   covariates. or sample mean if approximate method is selected.
#'   
#' @param thetahat  Estimator (vector or matrix) of \code{theta} for the 
#'   Relative Risk function.
#'   
#' @param tol   Tolerance for concluding numeric equality.
#'   
#' @return boolean Indicating \code{TRUE} if relative risk \code{rr} is as desired.
#'   
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'   
#' @seealso \code{\link{check.confidence}}, \code{\link{check.thetas}}, 
#'   \code{\link{check.cft}}, \code{\link{check.xvar}}, 
#'   \code{\link{check.exposure}}, \code{\link{check.integrals}}
#'   
#' @examples 
#' X  <- data.frame(rnorm(100))
#' rr <- function(X, theta){exp(X*theta)}
#' check.rr(X, 1, rr)
#' 
#' @keywords internal
#' 
#' @export

check.rr <- function(X, thetahat,  rr, tol = 1.e-8){
  
  #Boolean variable = 1
  .bool <- TRUE
  
  #Check if X is numeric
  if(is.numeric(as.matrix(X))){
  
    #Create matrix of size 0
    .X0           <- as.data.frame(matrix(0, ncol = ncol(X), nrow = 1))
    colnames(.X0) <- colnames(X)
    
    #Check condition
    if (  norm(as.matrix(rr(.X0, thetahat)) - 1, type = "2") > tol) {
      .bool <- FALSE
      warning(paste("Relative Risk by definition must equal 1 when evaluated in 0.",
                    "Are you using displaced RRs?"))
    }
    
  } 
  
  return(.bool)
  
}
