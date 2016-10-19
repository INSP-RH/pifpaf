#' @title Potential Impact Fraction
#' 
#' @description Function that calculates the potential impact fraction via the empirical method
#' 
#' @param X         Random sample (can be vector or matrix) which includes exposure and covariates.
#' 
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function
#' 
#' @param rr        Function for relative risk
#' 
#' 
#' **Optional**
#' 
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'                  the Population Attributable Fraction \code{PAF} where counterfactual is 0 exposure
#' 
#' @param weights   Survey \code{weights} for the random sample \code{X}
#' 
#' @param method    Either \code{empirical} (default) or \code{kernel}. 
#' 
#' @param ktype    \code{kernel} type from  \code{gaussian}, \code{epanechnikov}, \code{rectangular},
#'                  \code{triangular}, \code{biweight}, \code{cosine}, \code{optcosine}
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
#' X <- rnorm(100)
#' thetahat <- 0.12
#' pif(X, thetahat, function(X, theta){exp(theta*X)})
#' 
#' #Same example with kernel method
#' pif(X, thetahat, function(X, theta){exp(theta*X)}, method = "kernel")
#' 
#' #Same example considering counterfactual of halfing exposure
#' pif(X, thetahat, function(X, theta){exp(theta*X)}, cft = function(X){ 0.5*X }, method = "empirical")
#' 
#' #Example with linear relative risk
#' pif(X, thetahat, function(X, theta){theta*X + 1}, cft = function(X){ 0.5*X })
#' 
#' #Multivariate example
#' set.seed(18427)
#' X1 <- rnorm(1000)
#' X2 <- rnorm(1000)
#' X  <- as.matrix(cbind(X1,X2))
#' thetahat <- c(0.12, 0.03)
#' rr <- function(X, theta){exp(theta[1]*X[,1] + theta[2]*X[,2])}
#' pif(X, thetahat, rr) 
#' 
#' \dontrun{
#' #Multivariate cases cannot be evaluated with kernel method
#' pif(X, thetahat, rr, method = "kernel") 
#' }
#' 
#' 
#' @import matrixStats
#' 
#' @export


pif <- function(X, thetahat, rr, 
                  cft = function(Varx){matrix(0,ncol = ncol(as.matrix(Varx)), nrow = nrow(as.matrix(Varx)))}, 
                  weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), method = c("empirical", "kernel"),
                  ktype = "epanechnikov", bw = "nrd0", adjust = 1, npoints = 1000){
  
  #Get method from vector
  .method <- as.vector(method)[1]
  
  #Check that X has at most 1 column if kernel is chosen
  .X      <- as.matrix(X)
  if(ncol(.X) > 1 & .method == "kernel"){
    .method <- "empirical"
    warning("kernel method only works with one-dimensional X. Defaulting to empirical method.")
  }
  
  switch(.method,
    
    empirical = {
      .pif <- pif.empirical(X, thetahat, rr, cft, weights)
    }, 
    
    kernel    = {
      .pif <- pif.kernel(X, thetahat, rr, cft, weights, ktype, bw, adjust, npoints)
    },
    
    {
      warning("Please specify method as either empirical or kernel. Defaulting to empirical.")
      .pif <- pif.empirical(X, thetahat, rr, cft, weights)
    }
    
  )
  
  return(.pif)
}


