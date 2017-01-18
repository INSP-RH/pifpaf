#' @title Plot of Potential Impact Fraction under different values of theta (univariate)
#' 
#' @description Function that plots the PIF under different values of an univariate parameter theta. 
#' 
#' @param X         Random sample (can be vector or matrix) which includes exposure and covariates.
#' 
#' @param thetalow  Minimum of theta (parameter of relative risk \code{rr}) for plot
#' 
#' @param thetaup   Maximum of theta (parameter of relative risk \code{rr}) for plot
#' 
#' @param rr        Function for Relative Risk which uses parameter 
#'   \code{theta}. The order of the parameters shound be \code{rr(X, theta)}.
#' 
#' **Optional**
#' 
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'   the Population Attributable Fraction \code{\link{paf}} where counterfactual
#'   is 0 exposure.
#'   
#' @param weights   Normalized survey \code{weights} for the sample \code{X}.
#'   
#' @param method    Either \code{empirical} (default), \code{kernel} or \code{approximate}.
#' 
#' @param Xvar      Variance of exposure levels (for \code{approximate} method)
#'   
#' @param deriv.method.args \code{method.args} for 
#'   \code{\link[numDeriv]{hessian}} (for \code{approximate} method).
#'   
#' @param deriv.method      \code{method} for \code{\link[numDeriv]{hessian}}. 
#'   Don't change this unless you know what you are doing (for
#'   \code{approximate} method).
#'   
#' @param ktype    \code{kernel} type:  \code{"gaussian"}, 
#'   \code{"epanechnikov"}, \code{"rectangular"}, \code{"triangular"}, 
#'   \code{"biweight"}, \code{"cosine"}, \code{"optcosine"} (for \code{kernel}
#'   method). Additional information on kernels in \code{\link[stats]{density}}
#'   
#' @param bw        Smoothing bandwith parameter from density (for \code{kernel}
#'   method) from \code{\link[stats]{density}}. Default \code{"SJ"}.
#'   
#' @param adjust    Adjust bandwith parameter from density (for \code{kernel}
#'   method) from \code{\link[stats]{density}}.
#'   
#' @param n   Number of equally spaced points at which the density (for
#'   \code{kernel} method) is to be estimated (see
#'   \code{\link[stats]{density}}).
#'   
#' @param mpoints Number of points in plot   
#' 
#' @param color Colour of plot
#' 
#' @param xlab Label of x-axis in plot
#' 
#' @param ylab Label of y-axis in plot
#' 
#' @param title Title of plot
#'   
#' @param check_integrals Check that counterfactual and relative risk's expected
#'   values are well defined for this scenario
#'   
#' @param check_exposure  Check that exposure \code{X} is positive and numeric
#'   
#' @param check_rr        Check that Relative Risk function \code{rr} equals 
#'   \code{1} when evaluated at \code{0}
#' 
#' @return pif.plot       \code{\link[ggplot2]{ggplot}} object with plot of Potential Impact Fraction 
#' as function of \code{theta}
#'   
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí \email{daliaf172@gmail.com}
#' 
#' @import ggplot2
#' 
#' @examples 
#' 
#' #Example 1: Exponential Relative Risk empirical method
#' #-----------------------------------------------------
#' set.seed(18427)
#' X <- rnorm(100, 4.2, 1.3)
#' pif.plot(X, thetalow = 0, thetaup = 2, function(X, theta){exp(theta*X)})
#' 
#' #Same example with kernel method
#' pif.plot(X, 0, 2, function(X, theta){exp(theta*X)}, method = "kernel",
#' title = "Kernel method example") 
#'  
#' #Same example for approximate method
#' Xmean <- mean(X)
#' Xvar  <- var(X)
#' pif.plot(Xmean, 0, 2, function(X, theta){exp(theta*X)}, 
#' method = "approximate", Xvar = Xvar, title = "Approximate method example")
#' 
#' #Example with counterfactual
#' pif.plot(X, 0, 2, function(X, theta){exp(theta*X)}, cft = function(X){sqrt(X)})
#' 
#' #Example for approximate method with square root counterfactual
#' pif.plot(Xmean, 0, 2, function(X, theta){exp(theta*X)},  cft = function(X){sqrt(X)},
#'  method = "approximate", Xvar = Xvar) 
#' 
#' @export

pif.plot <- function(X, thetalow, thetaup, rr,         
                     cft = function(Varx){matrix(0,ncol = ncol(as.matrix(Varx)), nrow = nrow(as.matrix(Varx)))},  #Counterfactual
                     weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                     method  = c("empirical", "kernel", "approximate"),
                     adjust = 1, n = 512, mpoints = 100,
                     Xvar    = var(X), 
                     deriv.method.args = list(), 
                     deriv.method      = c("Richardson", "complex"),
                     ktype  = c("gaussian", "epanechnikov", "rectangular", "triangular", 
                                "biweight","cosine", "optcosine"), 
                     bw     = c("SJ", "nrd0", "nrd", "ucv", "bcv"),
                     color = "darkslategrey", xlab = "Theta", ylab = "PIF",
                     title = "Potential Impact Fraction (PIF) under different values of theta",
                     check_exposure = TRUE, check_rr = TRUE, check_integrals = TRUE){
  
  #Check thetas are univariate
  if(length(thetalow) > 1 || length(thetaup) > 1){
    stop("pif.plot only works for rr with unidimensional theta")
  }
  
  #Check that we are able to plot
  if (thetalow >= thetaup){
    stop("Minimum thetalow cannot be equal or greater than maximum thetaup")
  }
    
    #Create sequence from thetalow to thetaup
    .theta <- seq(thetalow, thetaup, length.out = ceiling(mpoints))
    
    #Create data frame for saving values of theta
    .dat   <- matrix(NA, nrow = mpoints, ncol = 2)
    colnames(.dat) <- c("Theta","PIF")
    
    #Loop through values of theta for plot
    for (i in 1:mpoints){
      
      #Save theta value
      .dat[i,"Theta"] <- .theta[i]
      
      #Calculate PIF
      .dat[i,"PIF"]   <- pif(X = X, .theta[i], rr = rr, cft = cft,
                             weights =  weights, method = method, 
                             Xvar = Xvar, deriv.method.args = deriv.method.args,
                             deriv.method = deriv.method,
                             adjust = adjust, n = n,ktype  = ktype, bw     = bw,
                             check_exposure = check_exposure, check_rr = check_rr, 
                             check_integrals = check_integrals) 
      
    }
    
    #Create plot
    .thetaplot <- ggplot(as.data.frame(.dat)) + 
                  geom_path(aes(x = .dat[,"Theta"], y = .dat[,"PIF"]), color = color) +
                  xlab(xlab) + theme_bw() + ylab(ylab) +
                  ggtitle(title)
    
    return(.thetaplot)

}
