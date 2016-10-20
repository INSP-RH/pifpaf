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
  
  switch(method, 
         
         one2one = {
           
           #Check that thetalow and thetaup exist
           if (is.na(thetalow) || is.na(thetaup)){
             stop("You have not correctly specified the bounds of the interval of confidence of theta")
           }
           
           #Check that thetahat < thetaup
           if (thetaup < thetahat){
             stop("thetaup < thetahat Please check that you correctly specified the interval of confidence of theta")
           }
           
           #Check that thetalow < thetahat
           if (thetahat < thetalow){
             stop("thetahat < thetalow. Please check that you correctly specified the interval of confidence of theta")
           }
         },
         
         {
           
           #Check that variance exists and is non-negative
           if (is.na(thetavar[1])){
             stop("Please specify variance of theta")
           }
           
           #Check that is positive semidefinite
           if (is.negative.definite(thetavar)){
             stop("Variance is not positive semi-definite")
           }
           
         }
         
  )
  
  
}