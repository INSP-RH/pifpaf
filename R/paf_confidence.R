#' @title Confidence intervals for the Population Attributable Fraction
#' 
#' @description Confidence intervals for the population attributable fraction 
#' 
#' @param X         Random sample (can be vector or matrix) which includes exposure and covariates.
#' 
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function
#' 
#' @param rr        Function for Relative Risk which uses parameter 
#'   \code{theta}. The order of the parameters shound be \code{rr(X, theta)}.
#' 
#' 
#' **Optional**
#' @param thetavar   Estimator of variance of thetahat
#' 
#' @param thetalow  Lower bound of the confidence interval (vector)
#'   
#' @param thetaup   Upper band of the confidence interval (vector)
#' 
#' @param weights    Survey \code{weights} for the random sample \code{X}.
#' 
#' @param nsim      Number of simulations for estimation of variance.
#' 
#' @param confidence Confidence level \% (default 95)
#' 
#' @param confidence_method  Either \code{inverse}, \code{one2one}, \code{linear}, \code{loglinear}, \code{bootstrap}
#' 
#' @param method    Either \code{empirical} (default), \code{kernel} or 
#'   \code{approximate}.
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
#' @param check_thetas Check that theta parameters are correctly inputed
#' 
#' @param check_exposure  Check that exposure \code{X} is positive and numeric
#' 
#' @param check_cft  Check if counterfactual function \code{cft} reduces exposure.
#' 
#' @param check_xvar Check if it is covariance matrix.
#'
#' @param check_integrals Check that counterfactual and relative risk's expected
#'   values are well defined for this scenario
#' @param check_rr        Check that Relative Risk function \code{rr} equals 
#'   \code{1} when evaluated at \code{0}
#'   
#' @param force.min Boolean indicating whether to force the \code{rr} to have a 
#'                  minimum value of 1 instead of 0 (not recommended).
#'
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí \email{daliaf172@gmail.com}
#' 
#' @examples 
#' # Example 1: univariate example
#' set.seed(82392)
#' X         <- rnorm(1000,3,.5)
#' thetahat  <- 0.5
#' thetavar  <- 0.1
#' rr        <- function(X,theta){exp(theta*X)}
#' # Default
#' paf.confidence(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr)
#' # One to one
#' paf.confidence(X = X, thetahat = thetahat, thetalow = 0.45, 
#' thetaup = 0.55, rr = rr, confidence_method = "one2one")
#' # Approximate 
#' Xmean <- mean(X)
#' Xvar  <- var(X)
#' paf.confidence(X = Xmean, thetahat = thetahat, thetavar = thetavar, 
#' rr = rr, method = "approximate", Xvar = Xvar)
#'  
#' # Example 2: multivariate example
#'  
#' X1 <- rnorm(100, 3,.5)
#' X2 <- rnorm(100,3,.5)
#' X  <- as.matrix(cbind(X1,X2))
#' thetahat <- c(0.1, 0.03)
#' thetavar <- matrix(c(0.1, 0, 0, 0.05), byrow = TRUE, nrow = 2)
#' rr       <- function(X, theta){
#'   .X <- matrix(X, ncol = 2)
#'   exp(theta[1]*.X[,1] + theta[2]*.X[,2])
#' }
#' cft <- function(X){0.5*X}
#' 
#' # Default
#' paf.confidence(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr)
#' 
#'  # Approximate 
#'  Xmean <- t(as.matrix(colMeans(X)))
#'  Xvar  <- cov(X)
#'  paf.confidence(X = Xmean, thetahat = thetahat, thetavar = thetavar, 
#'  rr = rr, method = "approximate", Xvar = Xvar)
#' # One to one
#' paf.confidence(X = X, thetahat = thetahat, thetalow = c(0.05, 0), 
#' thetaup = c(0.15, 0.08), rr = rr, confidence_method = "one2one")
#' 
#' @export

paf.confidence <- function(X, thetahat, rr,  thetavar = matrix(0, ncol = length(thetahat), nrow = length(thetahat)),
                           thetalow = thetahat, thetaup = thetahat,
                           weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                           nsim    =  100, confidence = 95,
                           confidence_method = c("bootstrap", "inverse", "one2one", "loglinear", "linear"),
                           method  = c("empirical", "kernel", "approximate"),
                           Xvar    = var(X), 
                           deriv.method.args = list(), 
                           deriv.method      = c("Richardson", "complex"),
                           adjust = 1, n = 512,
                           ktype  = c("gaussian", "epanechnikov", "rectangular", "triangular", 
                                      "biweight","cosine", "optcosine"), 
                           bw     = c("SJ", "nrd0", "nrd", "ucv", "bcv"),
                           check_exposure = TRUE, check_cft = TRUE, check_rr = TRUE,
                           check_xvar = TRUE, check_integrals = TRUE, check_thetas = TRUE,
                           force.min = FALSE){
  
  method            <- method[1]
  confidence_method <- confidence_method[1]
  
  if(confidence_method == "linear" || confidence_method == "loglinear" || confidence_method == "bootstrap" || method == "kernel"){
    pif.confidence(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr, weights = weights, 
                   nsim    =  nsim, confidence = confidence,
                   confidence_method = confidence_method, method  = method,
                   Xvar    = Xvar, deriv.method.args = deriv.method.args, 
                   deriv.method  = deriv.method, adjust = adjust, n = n,
                   ktype  = ktype,  bw = bw, check_exposure = check_exposure,
                   check_cft = check_cft, check_rr = check_rr,
                   check_xvar = check_xvar, check_integrals = check_integrals,
                   check_thetas = check_thetas, is_paf = TRUE)
  }else{
    switch (confidence_method,
            "inverse" ={
              paf.confidence.inverse(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr, 
                                     weights = weights, nsim = nsim, confidence = confidence,
                                     deriv.method.args = deriv.method.args, deriv.method = deriv.method, 
                                     force.min = force.min,  check_thetas = check_thetas,
                                     method = method, Xvar = Xvar)
            },
            "one2one" ={
              paf.confidence.one2one(X = X, thetahat = thetahat, thetalow = thetalow, thetaup = thetaup, rr = rr, 
                                     weights =  weights, confidence = confidence,
                                     check_thetas = check_thetas, deriv.method.args = deriv.method.args,
                                     deriv.method = deriv.method, method = method, Xvar = Xvar,
                                     check_exposure = check_exposure, check_rr = check_rr, 
                                     check_integrals = check_integrals)
            },
            {warning("Method of confidence interval estimation defaulted to inverse.")
              paf.confidence.inverse(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr, 
                                     weights = weights, nsim = nsim, confidence = confidence,
                                     deriv.method.args = deriv.method.args, deriv.method = deriv.method, 
                                     force.min = force.min,  check_thetas = check_thetas,
                                     method = method, Xvar = Xvar)
            }
    )
    
  }
}
