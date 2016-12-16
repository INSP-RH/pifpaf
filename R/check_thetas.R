#' @title Check thetas against method
#' 
#' @description Function for checking that the thetas are correctly specified according to chosen method
#' 
#' @param thetavar  Variance of theta
#' 
#' @param thetahat  Point estimate of theta
#' 
#' @param thetalow  Lower bound of theta's CI
#' 
#' @param thetaup   Upper bound of theta's CI
#' 
#' @param method    Method of CI's 
#' 
#' @return bool     Boolean variable indicating if hypothesis are matched
#' 
#' @import matrixcalc
#' @export

check.thetas <- function(thetavar, thetahat, thetalow, thetaup, method){
  
  #Boolean default true
  bool <- TRUE
  
  if(is.na(thetahat[1])){
    stop("Thetahat wasn't specified")
  }
  
  switch(method, 
         
         one2one = {
           
           #Check that thetalow and thetaup exist
           if (is.na(thetalow) || is.na(thetaup)){
             stop("You have not correctly specified the bounds of the interval of confidence of theta")
           }
          .thetalow <- as.matrix(thetalow)
          .thetaup  <- as.matrix(thetaup)
          .thetahat <- as.matrix(thetahat)
           
          if((dim(.thetahat)!=dim(.thetalow) || dim(.thetahat)!=dim(.thetaup))){
            stop("Dimensions of thetahat, thetalow, and thetaup are not the same.")
          }
           #Check that thetahat < thetaup
          correct <- TRUE
          i       <- 1
          while(correct && i <= length(.thetahat)){
            if (.thetaup[i] < .thetahat[i]){
              correct <- FALSE
              stop("thetaup < thetahat Please check that you correctly specified the interval of confidence of theta")
            }
            i <- i+1
            
          }
           
          correct <- TRUE
          i       <- 1
          while(correct && i <= length(.thetahat)){
            if (.thetalow[i] > .thetahat[i]){
              correct <- FALSE
              stop("thetalow > thetahat Please check that you correctly specified the interval of confidence of theta")
            }
            i <- i+1
            
          }
       
         },
         
         {
           
           #Check that variance exists and is non-negative
           if (is.na(thetavar[1])){
             stop("Please specify variance of theta")
           }
           
           if(length(thetahat)^2 != length(thetavar)){
             stop("Variance of theta must be of dimension nxn, where n is the length of thetahat")
           }
           
           #Check that is positive semidefinite
           if (is.square.matrix(as.matrix(thetavar)) == FALSE){
             stop("Variance and covariance matrix must be a square matrix")
           }
            
           if (is.symmetric.matrix(as.matrix(thetavar)) == FALSE){
             stop("Variance and covariance matrix must be symetric")
           }
           
           if (is.positive.semi.definite(as.matrix(thetavar)) == FALSE){
             stop("Variance must be positive semi-definite")
           }
           
           
         }
         
  )
  return(bool)
  
}
