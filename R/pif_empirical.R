#'@title Point Estimate of the Potential Impact Fraction via the Empirical 
#'  Method
#'  
#'@description Function that calculates the potential impact fraction \code{\link{pif}}
#'  via the empirical method. That is: for a random sample \code{X}, a relative
#'  risk function \code{rr(X, thetahat)} with parameters \code{thetahat} the
#'  empirical estimator is given by:
#'  
#'  \deqn{PIF = 1 - \sum rr(cft(X_i); \theta)/\sum rr(X_i; \theta)}
#'  
#'@param X         Random sample (\code{data.frame}) which includes exposure and 
#'  covariates. or sample mean if approximate method is selected.
#'  
#'@param thetahat  Estimator (\code{vector}) of \code{theta} for the Relative
#'  Risk function.
#'  
#'@param rr        Function for Relative Risk which uses parameter \code{theta}.
#'  The order of the parameters shound be \code{rr(X, theta)}.
#'  
#'  **Optional**
#'  
#'@param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'  the Population Attributable Fraction \code{\link{paf}} where counterfactual 
#'  is 0 exposure.
#'  
#'@param weights   Normalized survey \code{weights} for the sample \code{X}.
#'  
#'@param check_integrals Check that counterfactual and relative risk's expected 
#'  values are well defined for this scenario.
#'  
#'@param check_exposure  Check that exposure \code{X} is positive and numeric.
#'  
#'@param check_rr        Check that Relative Risk function \code{rr} equals 
#'  \code{1} when evaluated at \code{0}.
#'  
#' @param is_paf    Boolean forcing evaluation of \code{\link{paf}}.
#'  
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'  
#'@note The empirical method converges for relative risk \code{rr} functions 
#'  that are Lipschitz, convex or concave on \code{thetahat}. For stranger 
#'  functions use \code{\link{pif.kernel}}.
#'  
#'@seealso \code{\link{pif}} which is a wrapper for all pif methods 
#'  (\code{\link{pif.empirical}}, \code{\link{pif.approximate}}, 
#'  \code{\link{pif.kernel}}).
#'  
#'  For estimation of the Population Attributable Fraction see
#'  \code{\link{paf}}.
#'  
#' @examples 
#' 
#' #Example 1: Relative risk given by exponential
#'#--------------------------------------------
#' set.seed(18427)
#' X        <- data.frame(rnorm(100,3,.5))
#' thetahat <- 0.12
#' rr       <- function(X, theta){exp(theta*X)}
#' pif.empirical(X, thetahat, rr, cft = function(X){ 0.5*X })
#' 
#' #Without counterfactual estimates PAF
#' pif.empirical(X, thetahat, rr) 
#'  
#' #Example 2: Linear relative risk
#' #--------------------------------------------
#' pif.empirical(X, thetahat, rr = function(X, theta){theta*X + 1}, 
#'                cft = function(X){ 0.5*X })
#' 
#' #Example 3: Multivariate relative risk
#' #--------------------------------------------
#' set.seed(18427)
#' X1       <- rnorm(100,4,1)
#' X2       <- rnorm(100,2,0.4)
#' X        <- data.frame(cbind(X1,X2))
#' thetahat <- c(0.12, 0.03)
#' rr       <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' 
#' #Creating a counterfactual. As rr requires a bivariate input, cft should 
#' #return a two-column matrix
#' cft  <- function(X){
#'    cbind(X[,1]/2, 1.1*X[,2])
#' }
#' pif.empirical(X, thetahat, rr, cft) 
#' 
#' @keywords internal
#'@export


pif.empirical <- function(X, thetahat, rr, 
                          cft = NA, 
                          weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))),
                          check_exposure = TRUE, check_rr = TRUE, check_integrals = TRUE,
                          is_paf = FALSE){
  
  #Set X as data frame
  .X  <- as.data.frame(X)
  
  #Check exposure values are greater than zero
  if(check_exposure){ check.exposure(.X) }
  
  #Check that rr is 1 when X = 0
  if(check_rr){ check.rr(.X, thetahat, rr) }
  
  #Check cft exists
  if(!is.function(cft)){ is_paf <- TRUE}
  
  #Estimate weighted sums
  .mux   <- weighted.mean(as.matrix(rr(.X, thetahat)), weights)
  if (is_paf){
    .mucft <- 1
  } else {
    .mucft <- weighted.mean(as.matrix(rr(cft(.X), thetahat)), weights)  
  }
  
  #Check that integrals make sense
  if(check_integrals){ check.integrals(.mux, .mucft) }
  
  #Calculate PIF
  if(is.infinite(.mux)){   warning("Expected value of Relative Risk is not finite") }
  if(is.infinite(.mucft)){ warning(paste("Expected value of Relative Risk under", 
                                         "counterfactual is not finite")) }
  
  .pif   <- 1 - .mucft/.mux

  return(.pif)
  
}


