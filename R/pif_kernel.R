#' @title Kernel-based Potential Impact Fraction
#' 
#' @description Function that calculates the potential impact fraction for a parameter theta based on kernel methods
#' 
#' @param X         Random sample (vector) of exposure
#' 
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function
#' 
#' @param rr        Function for relative risk
#' 
#' 
#' **Optional**
#' 
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for \code{PAF} where counterfactual is 0 exposure
#' 
#' @param cft.check Boolean indicating to check if the counterfactual reduces the exposure \code{X} or does not.
#' 
#' @param weights   Survey \code{weights} for the random sample \code{X}
#' 
#' @param kernel    Kernel type from  "gaussian", "epanechnikov", "rectangular","triangular", "biweight","cosine", "optcosine"
#' 
#' @param bw        Smoothing bandwith parameter from density
#' 
#' @param adjust    Adjust bandwith parameter from density
#' 
#' @param npoints   Number of points
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí 
#' 
#' @examples 
#' 
#' #Example with risk given by HR
#' set.seed(18427)
#' X <- rnorm(100,4,1)
#' thetahat <- 0.12
#' pif.kernel(X, thetahat, function(X, theta){exp(theta*X)})
#' 
#' #Same example considering counterfactual of halfing exposure
#' pif.kernel(X, thetahat, function(X, theta){exp(theta*X)}, cft = function(X){ 0.5*X })
#' 
#' 
#' @import sfsmisc
#' @export


pif.kernel <- function(X, thetahat, rr, 
                       cft = function(Varx){matrix(0,ncol = ncol(as.matrix(Varx)), nrow = nrow(as.matrix(Varx)))},
                       weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                       kernel = "epanechnikov", bw = "nrd0", adjust = 1, npoints = 1000, 
                       cft.check = TRUE){
  
  #Set X as matrix
  .X      <- as.matrix(X)
  
  #Check exposure values are greater than zero
  check.exposure(X)
  
  #Check that rr = 1 when evaluated in 0
  check.rr(.X, thetahat, rr)
  
  #Check counterfactual
  if(cft.check){ check.cft(cft, .X) }
  
  #Check that the number of points in integer > 0
  .npoints <- max(2, ceiling(npoints))
  if(npoints < 500){
    warning("You should include more points in 'npoints' for better integration")
  }
  
  #Check that X has only one column
  if (ncol(.X) > 1){
    stop("X has to be one-dimensional. For multidimensional estimation use empirical method")
  }
  
  #Get the kernel density of X
  #http://stats.stackexchange.com/questions/14061/area-under-the-pdf-in-kernel-density-estimation-in-r
  .fX <- density(X, bw = bw, adjust = adjust, kernel = kernel, weights = weights, n = .npoints)
  
  #Integrate
  .integraldown <- integrate.xy(.fX$x, .fX$y*rr(.fX$x, thetahat))
  .integralup   <- integrate.xy(.fX$x, .fX$y*rr(cft(.fX$x), thetahat))
  
  #Check that integrals make sense
  check.integrals(.integraldown, .integralup)
  
  #Calculate pif
  .pif          <- 1 - .integralup/.integraldown
  
  #Return variance
  return(.pif)
  
}


