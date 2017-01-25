#' @title Check confidence (of confidence interval) is between zero and a hundred
#' 
#' @description Function that verifies \code{confidence} level selected is > 0 and < 100
#' 
#' @param confidence     Confidence level desired with \code{0 < confidence < 100}
#' 
#' @return TRUE if \code{confidence} is well defined between zero and a hundred
#' 
#' @examples 
#' #Example 1 
#' confidence <- 95
#' check.confidence(confidence)
#' 
#' @author Rodrigo Zepeda Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho García Formentí \email{daliaf172@gmail.com}

#' 
#' @export

check.confidence <- function(confidence){
  
  #Boolean variable = 1
  bool <- TRUE
  
  #Check condition
  if(confidence <= 0 || confidence >= 100){
    bool <- FALSE
    stop("Confidence level incorrectly specified")
  }
    
  return(bool)
}
