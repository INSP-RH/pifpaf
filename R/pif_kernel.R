#'@title Kernel-based estimate of Potential Impact Fraction
#'  
#'@description Function for estimating the Potential Impact Fraction \code{\link{pif}}
#'  from a cross-sectional sample of the exposure \code{X} with a known Relative
#'  Risk function \code{rr} with parameters \code{theta} using kernels.
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
#' @param check_integrals Check that counterfactual and relative risk's expected
#'   values are well defined for this scenario.
#'   
#' @param check_exposure  Check that exposure \code{X} is positive and numeric.
#'   
#' @param check_rr        Check that Relative Risk function \code{rr} equals 
#'   \code{1} when evaluated at \code{0}.
#'  
#' @param is_paf    Boolean forcing evaluation of \code{\link{paf}}.
#'  
#'@return pif      Estimate of Potential Impact Fraction
#'  
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
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
#' X        <- data.frame(rnorm(100,3,.5))
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
#' #Example 2: Linear relative risk
#' #--------------------------------------------
#' pif.kernel(X, thetahat, rr = function(X, theta){theta*X + 1}, 
#'                cft = function(X){ 0.5*X })
#' 
#' #Example 3: More complex counterfactual
#' #--------------------------------------------
#' set.seed(18427)
#' X       <- data.frame(rnorm(100,4,1))
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
#'@keywords internal
#'@export


pif.kernel <- function(X, thetahat, rr, 
                       cft = NA,
                       weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                       adjust = 1, n = 512,
                       ktype = c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight","cosine", "optcosine"), 
                       bw = c("SJ", "nrd0", "nrd", "ucv", "bcv"),
                       check_exposure = TRUE, check_rr = TRUE, check_integrals = TRUE,
                       is_paf = FALSE){
  
  #Set X as matrix
  .X      <- as.matrix(X)
  
  #Get values for ktype and bw
  .ktype  <- as.vector(ktype)[1]
  .bw     <- as.vector(bw)[1]
  
  #Check if counterfactual is defined
  if (!is.function(cft)){ is_paf <- TRUE}
  
  #Check exposure values are greater than zero
  if(check_exposure){ check.exposure(.X) }
  
  #Check that rr is 1 when X = 0
  if(check_rr){ check.rr(.X, thetahat, rr) }
  
  #Check that the number of points in integer > 0
  .n <- max(2, ceiling(n))
  if(.n < 250){
    warning("For better integration include more points in 'n'.")
  }
  
  #Check that X has only one column
  if (ncol(.X) > 1){
    stop("X has to be one-dimensional. For multidimensional estimation use empirical method.")
  }
  
  #Get the kernel density of X
  #http://stats.stackexchange.com/questions/14061
  #/area-under-the-pdf-in-kernel-density-estimation-in-r
  .fX <- density(.X, bw = .bw, adjust = adjust, kernel = .ktype, weights = weights, n = .n)
  
  #Get x value of density as matrix
  densX <- as.matrix(.fX$x)
  densY <- as.vector(.fX$y)
  
  #Integrate the expected values of the densities.
  #Check that cft is well defined 
  .prod1   <- as.matrix(rr(densX, thetahat))
  .naprod1 <- which(is.na(.prod1))
  
  #Eliminate potential na's
  if(length(.naprod1) > 0){
    warning("Under this kernel density some Relative Risk values are NA")
    densX  <- densX[-.naprod1]
    densY  <- densY[-.naprod1]
    .prod1 <- .prod1[-.naprod1]
  }
  
  #Check expected value is finite
  if(any(is.infinite(.prod1))){   
    warning("Expected value of Relative Risk is not finite") 
    .mux <- Inf
  } else {
    #Estimate the lower integral
    .mux    <- integrate.xy(densX, densY*.prod1)
  }
  
  #Check if we are estimating a PAF or a PIF
  if (is_paf){
    .mucft   <- 1
  } else {
    .prod2   <- as.matrix(rr(cft(densX), thetahat))
    .naprod2 <- which(is.na(.prod2))
    
    #Eliminate potential na's
    if(length(.naprod2) > 0){
      warning("Under this kernel density some counterfactual values are NA")
      densX  <- densX[-.naprod2]
      densY  <- densY[-.naprod2]
      .prod2 <- .prod1[-.naprod2]
    }
    
    #Check expected value is finite
    if(any(is.infinite(.prod2))){ 
      warning("Expected value of Relative Risk under counterfactual is not finite") 
      .mucft <- Inf
    } else {
      .mucft <- integrate.xy(densX, densY*.prod2)
    }
  }  
  
  #Check that integrals make sense
  if(check_integrals){ check.integrals(.mux, .mucft) }
  
  #Get pif
  .pif          <- 1 - .mucft/.mux
  
  #Return variance
  return(.pif)
  
}


