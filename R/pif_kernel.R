#'@title Kernel-based estimate of Potential Impact Fraction
#'  
#'@description Function for estimating the Potential Impact Fraction \code{pif} 
#'  from a cross-sectional sample of the exposure \code{X} with a known Relative
#'  Risk function \code{rr} with parameters \code{theta} using kernels.
#'  
#'@param X         Random sample (vector or matrix) which includes exposure and 
#'  covariates. or sample mean if approximate method is selected.
#'  
#'@param thetahat  Estimator (vector or matrix) of \code{theta} for the Relative
#'  Risk function.
#'  
#'@param rr        Function for Relative Risk which uses parameter \code{theta}.
#'  The order of the parameters shound be \code{rr(X, theta)}.
#'  
#'  
#'  **Optional**
#'  
#'@param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'  the Population Attributable Fraction \code{\link{paf}} where counterfactual 
#'  is 0 exposure.
#'  
#'@param weights   Normalized survey \code{weights} for the sample \code{X}.
#'  
#'@param ktype    \code{kernel} type:  \code{"gaussian"}, \code{"epanechnikov"},
#'  \code{"rectangular"}, \code{"triangular"}, \code{"biweight"},
#'  \code{"cosine"}, \code{"optcosine"}. Additional information on kernels in
#'  \code{\link[stats]{density}}
#'  
#'@param bw        Smoothing bandwith parameter from density from 
#'  \code{\link[stats]{density}}. Default \code{"SJ"}
#'  
#'@param adjust    Adjust bandwith parameter from density from 
#'  \code{\link[stats]{density}}.
#'  
#'@param n   Number of equally spaced points at which the density is to be 
#'  estimated (see \code{\link[stats]{density}}).
#'  
#'@param check_integrals Check that counterfactual and relative risk's expected 
#'  values are well defined for this scenario
#'  
#'@param check_exposure  Check that exposure \code{X} is positive and numeric
#'  
#'@param check_rr        Check that Relative Risk function \code{rr} equals 
#'  \code{1} when evaluated at \code{0}
#'  
#'@return pif      Estimate of Potential Impact Fraction
#'  
#'@author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#'@author Dalia Camacho García Formentí \email{daliaf172@gmail.com}
#'  
#'@note In practice \code{\link{pif.empirical}} should be prefered as
#'  convergence is faster than \code{\link{pif.kernel}} for most functions. In
#'  addition, the scope of \code{\link{pif.kernel}} is limited as it does not
#'  work with multivariate exposure data \code{X}.
#'  
#'@seealso \code{\link{pif}} which is a wrapper for all pif methods 
#'  (\code{\link{pif.empirical}}, \code{\link{pif.approximate}}, 
#'  \code{\link{pif.kernel}}).
#'  
#'  For estimation of the Population Attributable Fraction see
#'  \code{\link{paf}}.
#'  
#'  For more information on kernels see \code{\link[stats]{density}}
#'  
#' @examples 
#' 
#' #Example 1: Relative risk given by exponential
#'#--------------------------------------------
#' set.seed(18427)
#' X        <- rnorm(100,3,.5)
#' thetahat <- 0.12
#' rr       <- function(X, theta){exp(theta*X)}
#' pif.kernel(X, thetahat, rr, cft = function(X){ 0.5*X })
#' 
#' #Choose a different kernel
#' pif.kernel(X, thetahat, rr, cft = function(X){ 0.5*X }, ktype = "gaussian")
#' 
#' #Specify kernel options
#' pif.kernel(X, thetahat, rr, cft = function(X){ 0.5*X }, ktype = "gaussian", 
#' bw = "nrd", adjust = 0.5, n = 1100)
#' 
#' #Without counterfactual estimates PAF
#' pif.kernel(X, thetahat, rr) 
#' 
#'  
#' #Example 2: Linear relative risk
#' #--------------------------------------------
#' pif.kernel(X, thetahat, rr = function(X, theta){theta*X + 1}, 
#'                cft = function(X){ 0.5*X })
#' 
#' #Example 3: More complex counterfactual
#' #--------------------------------------------
#' set.seed(18427)
#' X       <- rnorm(100,4,1)
#' thetahat <- c(0.12, 0.03)
#' rr       <- function(X, theta){1 + theta[1]*X + theta[2]*X^2}
#' 
#' #Creating a counterfactual. As rr requires a bivariate input, cft should 
#' #return a two-column matrix
#' cft  <- function(X){
#'    X[which(X > 4)] <- 1
#'    return(X)
#' }
#' pif.kernel(X, thetahat, rr, cft) 
#' 
#'@importFrom sfsmisc integrate.xy
#'@export


pif.kernel <- function(X, thetahat, rr, 
                       cft = function(Varx){matrix(0,ncol = ncol(as.matrix(Varx)), nrow = nrow(as.matrix(Varx)))},
                       weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                       adjust = 1, n = 512,
                       ktype = c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight","cosine", "optcosine"), 
                       bw = c("SJ", "nrd0", "nrd", "ucv", "bcv"),
                       check_exposure = TRUE, check_rr = TRUE, check_integrals = TRUE){
  
  #Set X as matrix
  .X      <- as.matrix(X)
  
  #Get values for ktype and bw
  .ktype  <- as.vector(ktype)[1]
  .bw     <- as.vector(bw)[1]
  
  #Check exposure values are greater than zero
  if(check_exposure){ check.exposure(.X) }
  
  #Check that rr is 1 when X = 0
  if(check_rr){ check.rr(.X, thetahat, rr) }
  
  #Check that the number of points in integer > 0
  .n <- max(2, ceiling(n))
  if(n < 250){
    warning("I suggest you should include more points in 'n' for better integration")
  }
  
  #Check that X has only one column
  if (ncol(.X) > 1){
    stop("X has to be one-dimensional. For multidimensional estimation use empirical method")
  }
  
  #Get the kernel density of X
  #http://stats.stackexchange.com/questions/14061/area-under-the-pdf-in-kernel-density-estimation-in-r
  .fX <- density(.X, bw = .bw, adjust = adjust, kernel = .ktype, weights = weights, n = .n)
  
  #Get x value of density as matrix
  densX <- as.matrix(.fX$x)
  densY <- as.vector(.fX$y)
  
  #Integrate
  .mux   <- integrate.xy(densX, densY*rr(densX, thetahat))
  .mucft <- integrate.xy(densX, densY*rr(cft(densX), thetahat))
  
  #Check that integrals make sense
  if(check_integrals){ check.integrals(.mux, .mucft) }
  
  #Calculate pif
  .pif          <- 1 - .mucft/.mux
  
  #Return variance
  return(.pif)
  
}


