#' @title Check Integrals
#' 
#' @description Function for checking that the integrals of \code{\link{pif}} are nonnegative.
#' 
#' @param meancft Mean of relative risk \code{rr} with counterfactual
#' 
#' @param meanrr  Mean of relative risk \code{rr} without counterfactual
#' 
#' @return bool   \code{TRUE} if as desired
#' 
#' @examples 
#' check.integrals(1,0)
#' 
#' \dontrun{
#' check.integrals(0,1)
#' check.integrals(1,-1)
#' }
#' 
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'   
#' @seealso \code{\link{check.confidence}}, \code{\link{check.thetas}}, 
#'   \code{\link{check.cft}}, \code{\link{check.xvar}}, 
#'   \code{\link{check.rr}}, \code{\link{check.exposure}}
#' 
#' @keywords internal
#' 
#' @export

check.integrals <- function(meanrr, meancft){
  
  #Boolean variable = 1
  .bool <- TRUE
  
  #Check that the relative risk has returned positive > 0 values
  if (meanrr <= 0){
    .bool <- FALSE
    warning("Expected value of relative risk is <= 0!")
  }
  
  #Check that the relative risk has returned positive > 0 values
  if (meancft < 0){
    .bool <- FALSE
    warning("Expected value of relative risk under counterfactual is < 0!")
  }
  
  #Check that counterfactual is reducing the risk on average
  if (meanrr < meancft){
    .bool <- FALSE
    warning(paste0("Counterfactual is increasing the Risk. Are you sure you",
                   " are specifying it correctly?"))
  }
  
  return(.bool)
  
}
