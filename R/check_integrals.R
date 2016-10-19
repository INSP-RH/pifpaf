#' @title Check Integrals
#' 
#' @description Function for checking that the integrals of PIF make sense
#' 
#' @param meancft Mean of relative risk with counterfactual
#' 
#' @param meanrr  Mean of relative risk without counterfactual
#' 
#' 
#' @return TRUE if as desired
#' 
#' @examples 
#' \dontrun{
#' check.integrals(0,1)
#' check.integrals(1,0)
#' check.integrals(1,-1)
#' }
#' 
#' @export

check.integrals <- function(meanrr, meancft){
  
  #Boolean variable = 1
  bool <- TRUE
  
  #Check that the relative risk has returned positive > 0 values
  if (meanrr <= 0){
    bool <- FALSE
    stop("Expected value of relative risk is <= 0!")
  }
  
  #Check that the relative risk has returned positive > 0 values
  if (meancft < 0){
    bool <- FALSE
    stop("Expected value of counterfactual relative risk is < 0!")
  }
  
  #Check that counterfactual is reducing the risk on average
  if (meanrr < meancft){
    bool <- FALSE
    warning("Counterfactual is increasing the Risk. Are you sure you are specifying it correctly?")
  }
  
  return(bool)
  
}
