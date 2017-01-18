#' @title Check Integrals
#' 
#' @description Function for checking that the integrals of \code{\link{pif}} make sense.
#' 
#' @param meancft Mean of relative risk with counterfactual
#' 
#' @param meanrr  Mean of relative risk without counterfactual
#' 
#' @return bool   TRUE if as desired
#' 
#' @examples 
#' check.integrals(1,0)
#' 
#' \dontrun{
#' check.integrals(0,1)
#' check.integrals(1,-1)
#' }
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí \email{daliaf172@gmail.com}

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
    warning("Counterfactual is increasing the Risk. Are you sure you are specifying it correctly?")
  }
  
  return(.bool)
  
}
