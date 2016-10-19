#' @title Check Relative Risk
#' 
#' @description Function for checking that relative risk equals 1 when evaluated in 0
#' 
#' @param rr    Relative risk function
#' 
#' @param X     Matrix of exposure that is evaluated in rr
#' 
#' @param thetahat Point estimate of theta
#' 
#' @param tol   Tolerance for equality
#' 
#' @return TRUE if relative risk is as desired
#' 
#' @examples 
#' check.rr(as.matrix(rnorm(100)), 1, function(X, theta){exp(X*theta)})
#' 
#' @export

check.rr <- function(X, thetahat,  rr, tol = 1.e-8){
  
  #Boolean variable = 1
  bool <- TRUE
  
  #Create matrix of size 0
  .X0 <- matrix(0, ncol = ncol(X), nrow = 1)
  
  #Check condition
  if (  abs(rr(.X0, thetahat) - 1) > tol) {
    bool <- FALSE
    warning(paste("Relative Risk by definition must equal 1 when evaluated in 0.",
                  "Are you using displaced RRs?"))
  }
  
  return(bool)
  
}
