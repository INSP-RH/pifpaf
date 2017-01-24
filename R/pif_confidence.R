#' @title Confidence intervals for the Potencial Impact Fraction
#' 
#' @description Confidence intervals for the potencial impact Fraction 
#' 
#' @param X         Random sample (can be vector or matrix) which includes exposure and covariates.
#' 
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function
#' 
#' @param thetavar   Estimator of variance of thetahat
#' 
#' @param rr        Function for Relative Risk which uses parameter 
#'   \code{theta}. The order of the parameters shound be \code{rr(X, theta)}.

#' 
#' 
#' **Optional**
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'                  the Population Attributable Fraction \code{PAF} where counterfactual is 0 exposure.
#' 
#' @param weights    Survey \code{weights} for the random sample \code{X}.
#' 
#' @param nsim      Number of simulations for estimation of variance.
#' 
#' @param confidence Confidence level \% (default 95)
#' 
#' @param confidence_method  Either \code{linear}, \code{loglinear}, \code{bootstrap}
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
#' @param is_paf    Boolean forcing evaluation of paf
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
#' cft <- function(X){0.5*X}
#' # Default
#' pif.confidence(X = X, thetahat = thetahat, thetavar = thetavar, 
#' rr = rr, cft = cft)
#' 
#' # Kernel
#' pif.confidence(X = X, thetahat = thetahat, thetavar = thetavar, 
#' rr = rr, cft = cft, method = "kernel",confidence_method = "bootstrap")
#' 
#'  # Approximate 
#'  Xmean <- mean(X)
#'  Xvar  <- var(X)
#'  pif.confidence(X = Xmean, thetahat = thetahat, thetavar = thetavar, 
#'  rr = rr, cft = cft, method = "approximate", Xvar = Xvar)
#'  
#'  # Example 2: multivariate example
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
#' pif.confidence(X = X, thetahat = thetahat, thetavar = thetavar, 
#' rr = rr, cft = cft)
#' 
#'  # Approximate 
#'  Xmean <- t(as.matrix(colMeans(X)))
#'  Xvar  <- cov(X)
#'  pif.confidence(X = Xmean, thetahat = thetahat, thetavar = thetavar, 
#'  rr = rr, cft = cft, method = "approximate", Xvar = Xvar)
#' 
#' @export

pif.confidence <- function(X, thetahat, thetavar, rr,         
                   cft = function(Varx){matrix(0,ncol = ncol(as.matrix(Varx)), nrow = nrow(as.matrix(Varx)))},  #Counterfactual
                   weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                   nsim    =  100, confidence = 95,
                   confidence_method = c( "bootstrap", "linear", "loglinear"),
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
                   is_paf = FALSE){
  
  method            <- method[1]
  confidence_method <- confidence_method[1]
  switch (method,
          "empirical"   = {
            switch (confidence_method,
                    "linear"    = {
                      pif.confidence.linear(X=X, thetahat = thetahat, thetavar = thetavar, rr = rr, cft = cft, weights = weights, 
                                            confidence = confidence, check_thetas = check_thetas, is_paf = is_paf)
                    },
                    "loglinear" = {
                      pif.confidence.loglinear(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr, cft = cft, weights = weights,
                                               nsim = nsim, check_thetas = check_thetas, check_exposure = check_exposure, 
                                               check_cft = check_cft, is_paf = is_paf)
                    },
                    "bootstrap" = {
                      pif.confidence.bootstrap(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr, cft = cft, weights = weights,
                                               method = "empirical", nboost = nsim, confidence = confidence,
                                               check_exposure = check_exposure, check_rr = check_rr,
                                               check_integrals = check_integrals, check_thetas = check_thetas,
                                               is_paf = is_paf)
                    },
                    
                    {warning("Method of confidence interval estimation defaulted to loglinear.")
                      pif.confidence.loglinear(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr, cft = cft, weights = weights,
                                               nsim = nsim, check_thetas = check_thetas, check_exposure = check_exposure, 
                                               check_cft = check_cft, is_paf = is_paf)
                    }
            )
          },
          "kernel"      = {
            if(confidence_method != "bootstrap"){
              warning("Method of confidence interval estimation defaulted to bootstrap.")
            }
            pif.confidence.bootstrap(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr, cft = cft, weights = weights,
                                     method = "kernel", nboost = nsim, adjust = adjust, confidence = confidence,
                                     ktype = ktype, bw = bw, check_exposure = check_exposure, check_rr = check_rr,
                                     check_integrals = check_integrals, check_thetas = check_thetas, is_paf = is_paf)
          },
          "approximate" = {
            switch(confidence_method,
              "linear"    = {
                pif.confidence.approximate(Xmean = X, Xvar = Xvar, thetahat = thetahat, thetavar = thetavar, rr = rr,
                                           cft = cft, check_thetas = check_thetas, check_cft = check_cft,
                                           check_xvar = check_xvar, check_rr = check_rr, check_integrals = check_integrals,
                                           check_exposure = check_exposure, deriv.method.args = deriv.method.args,
                                           deriv.method = deriv.method, nsim = nsim, confidence = confidence,
                                           is_paf = is_paf)
              },
              "loglinear" = {
                pif.confidence.approximate.loglinear(Xmean = X, Xvar = Xvar, thetahat = thetahat, thetavar = thetavar, rr = rr,
                                                     cft = cft, deriv.method.args = deriv.method.args, deriv.method = deriv.method,
                                                     check_exposure = check_exposure, check_rr = check_rr, check_integrals = check_integrals,
                                                     nsim = nsim, confidence = confidence, check_thetas = check_thetas,
                                                     is_paf = is_paf)
              },
              {warning("Method of confidence interval estimation defaulted to loglinear.")
                pif.confidence.approximate.loglinear(Xmean = X, Xvar = Xvar, thetahat = thetahat, thetavar = thetavar, rr = rr,
                                                     cft = cft, deriv.method.args = deriv.method.args, deriv.method = deriv.method,
                                                     check_exposure = check_exposure, check_rr = check_rr, check_integrals = check_integrals,
                                                     nsim = nsim, confidence = confidence, check_thetas = check_thetas,
                                                     is_paf = is_paf)
              }
            )
            
          },
          {warning("Method of PAF estimation defaulted to empirical")
            switch (confidence_method,
                    "linear"    = {
                      pif.confidence.linear(X=X, thetahat = thetahat, thetavar = thetavar, rr = rr, cft = cft, weights = weights, 
                                            confidence = confidence, check_thetas = check_thetas, is_paf = is_paf)
                    },
                    "loglinear" = {
                      pif.confidence.loglinear(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr, cft = cft, weights = weights,
                                               nsim = nsim, check_thetas = check_thetas, check_exposure = check_exposure, 
                                               check_cft = check_cft, is_paf = is_paf)
                    },
                    "bootstrap" = {
                      pif.confidence.bootstrap(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr, cft = cft, weights = weights,
                                               method = "empirical", nboost = nsim, confidence = confidence,
                                               check_exposure = check_exposure, check_rr = check_rr,
                                               check_integrals = check_integrals, check_thetas = check_thetas,
                                               is_paf = is_paf)
                    },
                    
                    {warning("Method of confidence interval estimation defaulted to loglinear.")
                      pif.confidence.loglinear(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr, cft = cft, weights = weights,
                                               nsim = nsim, check_thetas = check_thetas, check_exposure = check_exposure, 
                                               check_cft = check_cft, is_paf = is_paf)
                    }
            )
          }
  )
}