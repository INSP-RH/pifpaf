#' @title Check Exposure values
#' 
#' @description Function that verifies if all exposure values are greater or equal to zero (exposure values can not be negative)
#' 
#' @param X     Matrix of exposure that is evaluated in rr
#' 
#' @return TRUE if all exposure values are greater or equal to zero
#' 
#' @examples 
#' #Example 1 
#' X <- runif(500, 0,1)
#' check.exposure(X)
#' 
#' 
#' @export

check.exposure <- function(X){
  
  #Boolean variable = 1
  bool <- TRUE
  
  
  X <- as.matrix(X)
  m <- dim(X)[1]
  n <- dim(X)[2]
  
  #Check condition
  i <- 1
  while(i <= m && bool){
    j <- 1
    while(j <= n && bool){
      if(X[i,j] >= 0){
        j <- j+1
      }else{
        bool <- FALSE
        stop(paste("Some exposure values are less than zero, exposure can't be negative."))
      }
    }
    i <- i+1
  }
  return(bool)
}
