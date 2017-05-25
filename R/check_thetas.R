#' @title Check theta parameters for confidence intervals
#'   
#' @description Function for checking that the theta parameters are correctly
#'   specified according to chosen \code{pif.confidence} and
#'   \code{paf.confidence} \code{confidence_method}.
#'   
#' @param thetavar  Variance of \code{thetahat}.
#'   
#' @param thetahat  Point estimate of theta.
#'   
#' @param thetalow  Lower bound of theta's CI.
#'   
#' @param thetaup   Upper bound of theta's CI.
#'   
#' @param method    Method of CI's of  \code{pif.confidence} (resp. \code{paf.confidence})
#'   
#' @return bool     Boolean variable indicating if hypothesis are matched.
#'   
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'   
#' @seealso \code{\link{check.confidence}}, \code{\link{check.xvar}}, 
#'   \code{\link{check.cft}}, \code{\link{check.rr}}, 
#'   \code{\link{check.exposure}}, \code{\link{check.integrals}}
#'   
#' @importFrom matrixcalc is.positive.semi.definite is.symmetric.matrix
#'   is.square.matrix
#' 
#' @keywords internal
#'     
#' @export

check.thetas <- function(thetavar, thetahat, thetalow, thetaup, method){
  
  #Boolean default true
  .bool <- TRUE
  
  if(is.na(thetahat[1])){
    stop("Thetahat wasn't specified")
  }
  
  switch(method, 
         
         one2one = {
           
           #Check that thetalow and thetaup exist
           if (is.na(thetalow) || is.na(thetaup)){
             stop(paste0("The bounds are not correctly specified", 
                         "of the interval of confidence of theta"))
           }
          .thetalow <- as.vector(thetalow)
          .thetaup  <- as.vector(thetaup)
          .thetahat <- as.vector(thetahat)
           
          if((length(.thetahat)!=length(.thetalow) || length(.thetahat)!=length(.thetaup))){
            stop("Dimensions of thetahat, thetalow, and thetaup are not the same.")
          }
          
           #Check that thetahat < thetaup
          .correct <- TRUE
          .i       <- 1
          while(.correct && .i <= length(.thetahat)){
            if (.thetaup[.i] < .thetahat[.i] || .thetalow[.i] > .thetahat[.i] || 
                .thetaup[.i] < .thetahat[.i]){
              .correct <- FALSE
              stop(paste0("Thetas do not comply the inequality: 'thetaup > thetahat >", 
                          "thetalow. Verify theta's confidence interval is",
                          "correctly specified."))
            }
            .i <- .i + 1
            
          }
         },
         
         {
           
           #Check that variance exists and is non-negative
           if (is.na(thetavar[1])){
             stop("Please specify variance of theta")
           }
           
           if(length(thetahat) != nrow(as.matrix(thetavar))){
             stop("Covariance matrix dimensions must be nxn, where n is the length of thetahat")
           }
           
           #Check that is positive semidefinite
           if (is.square.matrix(as.matrix(thetavar)) == FALSE){
             stop("Covariance matrix must be a square matrix")
           }
            
           if (is.symmetric.matrix(as.matrix(thetavar)) == FALSE){
             stop("Covariance matrix must be symmetric")
           }
           
           if (is.positive.semi.definite(as.matrix(thetavar)) == FALSE){
             stop("Covariance must be positive semi-definite")
           }
           
           
         }
         
  )
  return(.bool)
  
}
