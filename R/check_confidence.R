#' @title Check confidence is between zero and less than a hundred
#' 
#' @description Function that verifies confidence level selected is between zero and less than a hundred
#' 
#' @param confidence     Confidence level desired
#' 
#' @return TRUE if Confidence is well defined between zero and a hundred
#' 
#' @examples 
#' #Example 1 
#' confidence <- 95
#' check.confidence(confidence)
#' 
#' 
#' @export

check.confidence <- function(confidence){
  
  #Boolean variable = 1
  bool <- TRUE
  
  #Check condition
  if(confidence < 0 || confidence >= 100){
    bool <- FALSE
    stop("Confidence level incorrectly specified")
  }
    
  return(bool)
}
