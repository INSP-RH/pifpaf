#' @title Check Relative Risk
#' 
#' @description Function for checking that relative risk equals 1 when evaluated in 0.
#' 
#' @param rr        Function for Relative Risk which uses parameter 
#' \code{theta}. The order of the parameters shound be \code{rr(X, theta)}.
#' 
#' @param X         Random sample (vector or matrix) which includes exposure and
#'   covariates. or sample mean if approximate method is selected.
#'   
#' @param thetahat  Estimator (vector or matrix) of \code{theta} for the 
#'   Relative Risk function.
#' 
#' @param tol   Tolerance for concluding numeric equality.
#' 
#' @return boolean Indicating \code{TRUE} if relative risk is as desired.
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí \email{daliaf172@gmail.com}
#' 
#' @examples 
#' check.rr(as.matrix(rnorm(100)), 1, function(X, theta){exp(X*theta)})
#' 
#' @export

check.rr <- function(X, thetahat,  rr, tol = 1.e-8){
  
  #Boolean variable = 1
  .bool <- TRUE
  
  #Check if X is numeric
  if(is.numeric(X)){
  
    #Create matrix of size 0
    .X0 <- matrix(0, ncol = ncol(X), nrow = 1)
    
    #Check condition
    if (  norm(rr(.X0, thetahat) - 1, type = "2") > tol) {
      .bool <- FALSE
      warning(paste("Relative Risk by definition must equal 1 when evaluated in 0.",
                    "Are you using displaced RRs?"))
    }
    
  } else {
    warning(paste("Relative Risks were not checked as input was not numeric.", 
                  "Please set check.rr = FALSE to avoid this message "))
  }
  
  return(.bool)
  
}
