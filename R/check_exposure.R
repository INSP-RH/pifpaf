#' @title Check Exposure values
#' 
#' @description Function that verifies if all exposure values are greater or equal than zero.
#' 
#' @param X     Data frame of exposure that is evaluated in relative risk \code{rr}.
#' 
#' @examples 
#' #Example 1 
#' X <- matrix(runif(500, 0,1))
#' check.exposure(X)
#' 
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'   
#' @seealso \code{\link{check.confidence}}, \code{\link{check.thetas}}, 
#'   \code{\link{check.cft}}, \code{\link{check.xvar}}, 
#'   \code{\link{check.rr}}, \code{\link{check.integrals}}
#' 
#' @keywords internal
#' 
#' @export

check.exposure <- function(X){
  
  #Boolean variable = 1
  .bool <- TRUE
  
  #Get column numbers and row numbers
  .m <- nrow(X)
  .n <- ncol(X)
  
  #Check condition for all rows
  .i <- 1
  while(.i <= .m & .bool){
    
    #Loop through all columns
    .j <- 1
    while(.j <= .n & .bool){
      
      #Check that is positive numeric
      if(is.numeric(X[.i, .j]) & as.numeric(X[.i, .j]) >= 0){
        .j <- .j + 1
      } else if (is.character(X[.i, .j]) || is.factor(X[.i, .j]) || is.logical(X[.i, .j])){
        .j <- .j + 1
      } else {
        .bool <- FALSE
        warning(paste0("Some exposure values are less than zero, ", 
                       "verify this is correct."))
      }
    }
    .i <- .i + 1
  }
  
  return(.bool)
}
