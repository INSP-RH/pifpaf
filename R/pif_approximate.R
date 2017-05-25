#'@title Point Estimate of the Potential Impact Fraction via the approximate 
#'  method.
#'  
#'@description Function for estimating the Potential Impact Fraction
#'  \code{\link{pif}} from a cross-sectional sample of the exposure \code{X}
#'  with a known Relative Risk function \code{rr} with parameters \code{theta}
#'  when only \code{mean(X)} and \code{var(X)} are known.
#'  
#'@param X      Mean value of exposure levels from a cross-sectional random 
#'  sample. If multivariate, this should be a \code{matrix} with 1 row and as 
#'  many columns as covariates
#'  
#'@param Xvar   Variance of the exposure levels.
#'  
#'@param thetahat  Estimator (\code{vector}) of \code{theta} for the Relative 
#'  Risk function.
#'  
#'@param rr        \code{function} for Relative Risk which uses parameter 
#'  \code{theta}. The order of the parameters shound be \code{rr(X, theta)}.
#'  
#' \strong{**Optional**}
#'  
#'@param cft       Twice differentiable function \code{cft(X)} for
#'  counterfactual. Leave empty for the Population Attributable Fraction
#'  \code{\link{paf}} where counterfactual is 0 exposure.
#'  
#'@param deriv.method.args \code{method.args} for 
#'  \code{\link[numDeriv]{hessian}}.
#'  
#'@param deriv.method      \code{method} for \code{\link[numDeriv]{hessian}}. 
#'  Don't change this unless you know what you are doing.
#'  
#'@param check_integrals Check that counterfactual and relative risk's expected 
#'  values are well defined for this scenario.
#'  
#'@param check_exposure  Check that exposure \code{X} is positive and numeric.
#'  
#'@param check_rr        Check that Relative Risk function \code{rr} equals 
#'  \code{1} when evaluated at \code{0}.
#'  
#'@param is_paf    Boolean forcing evaluation of \code{\link{paf}}.
#'  
#'@author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#'@author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'  
#'  
#'@note \code{pif.approximate} method should be the last choice for the case
#'  when no information on the exposure \code{X} (except for mean and standard
#'  deviation) are given. In practice \code{\link{pif.empirical}} should be
#'  prefered.
#'  
#'@seealso \code{\link{pif}} which is a wrapper for all pif methods 
#'  (\code{\link{pif.empirical}}, \code{\link{pif.approximate}}, 
#'  \code{\link{pif.kernel}}).
#'  
#'  For estimation of the Population Attributable Fraction see
#'  \code{\link{paf}}.
#'  
#'@examples 
#'
#'#Example 1
#'#--------------------------------------------
#' X         <- data.frame(2)
#' thetahat  <- 0.12
#' Xvar      <- 0.2
#' rr        <- function(X,theta){exp(X*theta)}
#' cft       <- function(X){0.5*X}
#' pif.approximate(X, Xvar, thetahat, rr, cft)
#' 
#' #Change the derivative arguments 
#' pif.approximate(X, Xvar, thetahat, rr, cft, 
#'                deriv.method = "Richardson", 
#'                deriv.method.args = list(eps=1e-8, d=0.000001))
#' 
#' #When no counterfactual is specified paf is calculated
#' pif.approximate(X, Xvar, thetahat, rr)
#' 
#' #Example 2: Multivariate
#' #--------------------------------------------
#' X1        <- 2
#' X2        <- 1.1
#' X         <- data.frame(X1,X2)
#' Xvar      <- matrix(c(1,.4,.4,1),ncol = 2, byrow = TRUE)
#' cft       <- function(X){.25*X}
#' thetahat  <- c(0.12, 0.03)
#' rr        <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' pif.approximate(X, Xvar, thetahat, rr, cft)
#' 
#' #Example 3: More multivariate
#' #--------------------------------------------
#' X1       <- rnorm(1000,3,.5)
#' X2       <- rnorm(1000,4,1)
#' X        <- cbind(X1,X2)
#' Xmean    <- data.frame(t(colMeans(X)))
#' Xvar     <- var(X)
#' thetahat <- c(0.12, 0.17)
#' thetavar <- matrix(c(0.001, 0.00001, 0.00001, 0.004), byrow = TRUE, nrow = 2)
#' rr       <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' cft      <- function(X){cbind(sqrt(X[,1] + 0.2*X[,2]), X[,1])}
#' pif.approximate(Xmean, Xvar, thetahat, rr, cft)
#' 
#'@importFrom numDeriv hessian
#'
#'@keywords internal
#'    
#'@export


pif.approximate <- function(X, Xvar, thetahat, rr, 
                            cft = NA,
                            deriv.method.args = list(), 
                            deriv.method = c("Richardson", "complex"),
                            check_exposure = TRUE, check_rr = TRUE, check_integrals = TRUE,
                            is_paf = FALSE){
  
  #Change to matrix
  .Xvar   <- check.xvar(Xvar)
  .X      <- as.matrix(X, ncol = ncol(.Xvar)) #Matrix for hessian
  
  #Get derivative method
  .method <- as.vector(deriv.method)[1]
  
  #Check exposure values are greater than zero
  if(check_exposure){ check.exposure(.X) }
  
  #Check that rr is 1 when X = 0
  if(check_rr){ check.rr(.X, thetahat, rr) }
  
  #Check if counterfactual is given
  if(!is.function(cft)){is_paf <- TRUE}
  
  #Rewrite rr and counterfactual functions as functions of X only for 
  #numerical derivatives
  .rr_cft_fun <- function(.xmat){
    return(as.matrix(rr( cft(.xmat), thetahat)))
  } 
  
  .rr_fun_x   <- function(.xmat){
    return(as.matrix(rr(.xmat, thetahat)))
  }
  
  #Estimate weighted sums
  if (is_paf){
    .mucft <- 1
  } else {
    .hcft  <- hessian(.rr_cft_fun, .X, method = .method, method.args = deriv.method.args)
    .mucft <- .rr_cft_fun(.X) + 0.5*sum(.hcft*.Xvar)
  }
  
  .hrr   <- hessian(.rr_fun_x,   .X, method = .method, method.args = deriv.method.args)
  .mux   <- .rr_fun_x(.X)   + 0.5*sum(.hrr*.Xvar)
  
  
  #Check derivative exists
  if (is.na(.mux)){   stop("Hessian might not be defined for those values of rr")}
  if (is.na(.mucft)){ stop("Hessian might not be defined for those values of rr and cft")}
  
  #Check that Hessian approximation didn't change sign of RR
  if (.mux   <= 0){ warning("Hessian cannot approximate numerically rr(X, theta) correctly.") }
  if (.mucft <= 0){ warning("Hessian cannot approximate numerically rr(cft(X),theta) correctly.") }
  
  #Check that integrals make sense
  if(check_integrals){ check.integrals(.mux, .mucft) }
  
  #Calculate PIF
  if(is.infinite(as.numeric(.mux))){   warning(paste("Expected value of Relative Risk", 
                                                "is not finite")) }
  if(is.infinite(as.numeric(.mucft))){ warning(paste("Expected value of Relative Risk", 
                                                "under counterfactual is not finite")) }
  
  #Construct pif
  .pif   <- as.numeric(1 - .mucft/.mux)
  
  return(.pif)
  
}