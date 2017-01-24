#' @title Plot of Potential Impact Fraction under different values of theta (univariate)
#' 
#' @description Function that plots the PAF under different values of an univariate parameter theta. 
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
#' @param weights   Normalized survey \code{weights} for the sample \code{X}.
#'   
#' @param method    Either \code{empirical} (default), \code{kernel} or \code{approximate}.
#' 
#' @param confidence Confidence level \% (default 95)
#' 
#' @param confidence_method  Either \code{linear}, \code{loglinear}, \code{bootstrap}

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
#' @param colors Colour of plot
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
#' @return paf.plot       \code{\link[ggplot2]{ggplot}} object with plot of Potential Impact Fraction 
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
#' X <- rbeta(25, 4.2, 10)
#' paf.plot(X, thetalow = 0, thetaup = 2, function(X, theta){exp(theta*X)})
#' 
#' \dontrun{
#' #Same example with kernel method
#' paf.plot(X, 0, 2, function(X, theta){exp(theta*X)}, method = "kernel",
#' title = "Kernel method example") 
#'  
#' #Same example for approximate method
#' Xmean <- mean(X)
#' Xvar  <- var(X)
#' paf.plot(Xmean, 0, 2, function(X, theta){exp(theta*X)}, 
#' method = "approximate", Xvar = Xvar, title = "Approximate method example")
#' }
#' @export

paf.plot <- function(X, thetalow, thetaup, rr,         
                     weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                     method  = c("empirical", "kernel", "approximate"),
                     adjust = 1, n = 512, mpoints = 100,
                     Xvar    = var(X), 
                     deriv.method.args = list(), 
                     deriv.method      = c("Richardson", "complex"),
                     confidence = 95,
                     confidence_method = c("bootstrap", "linear", "loglinear"),
                     ktype  = c("gaussian", "epanechnikov", "rectangular", "triangular", 
                                "biweight","cosine", "optcosine"), 
                     bw     = c("SJ", "nrd0", "nrd", "ucv", "bcv"),
                     colors = c("deepskyblue", "gray25"), xlab = "Theta", ylab = "PAF",
                     title = "Population Attributable Fraction (PAF) under different values of theta",
                     check_exposure = TRUE, check_rr = TRUE, check_integrals = TRUE){
  
  pif.plot(X = X, thetalow = thetalow, thetaup = thetaup, rr = rr, weights = weights,
           method = method, adjust = adjust, n = n, mpoints = mpoints, Xvar = Xvar,
           deriv.method.args = deriv.method.args, deriv.method = deriv.method,
           confidence = confidence, confidence_method = confidence_method,
           ktype = ktype, bw = bw, colors = colors, xlab = xlab, ylab = ylab,
           title = title, check_exposure = check_exposure, check_rr = check_rr,
           check_integrals = check_integrals, is_paf = TRUE)
  
}
