#' @title Plot of Potential Impact Fraction under different values of Relative 
#'   Risk's parameter theta (univariate)
#'   
#' @description Function that plots the \code{\link{pif}} under different values
#'   of an univariate parameter theta of the relative risk function \code{rr} 
#'   which depends on the exposure \code{X} and a \code{theta} parameter 
#'   (\code{rr(X, theta)})
#'   
#' @param X         Random sample (\code{data.frame}) which includes exposure 
#'   and covariates.
#'   
#' @param thetalow  Minimum of \code{theta} (parameter of relative risk 
#'   \code{rr}) for plot.
#'   
#' @param thetaup   Maximum of \code{theta} (parameter of relative risk 
#'   \code{rr}) for plot.
#'   
#' @param rr        \code{function} for Relative Risk which uses parameter 
#'   \code{theta}. The order of the parameters shound be \code{rr(X, theta)}.
#'   
#'   **Optional**
#'   
#' @param cft       \code{function} \code{cft(X)} for counterfactual. Leave 
#'   empty for the Population Attributable Fraction \code{\link{paf}} where 
#'   counterfactual is 0 exposure.
#'   
#' @param weights   Normalized survey \code{weights} for the sample \code{X}.
#'   
#' @param method    Either \code{"empirical"} (default), \code{"kernel"} or 
#'   \code{"approximate"}.
#'   
#' @param confidence Confidence level \% (default \code{95})
#'   
#' @param confidence_method  Either \code{"linear"}, \code{"loglinear"}, 
#'   \code{"bootstrap"}
#'   
#' @param Xvar      Variance of exposure levels (for \code{"approximate"} 
#'   method)
#'   
#' @param deriv.method.args \code{method.args} for 
#'   \code{\link[numDeriv]{hessian}} (for \code{"approximate"} method).
#'   
#' @param deriv.method      \code{method} for \code{\link[numDeriv]{hessian}}. 
#'   Don't change this unless you know what you are doing (for 
#'   \code{"approximate"} method).
#'   
#' @param ktype    \code{"kernel"} type:  \code{"gaussian"}, 
#'   \code{"epanechnikov"}, \code{"rectangular"}, \code{"triangular"}, 
#'   \code{"biweight"}, \code{"cosine"}, \code{"optcosine"} (for \code{"kernel"}
#'   method). Additional information on kernels in \code{\link[stats]{density}}
#'   
#' @param bw        Smoothing bandwith parameter from density (for 
#'   \code{"kernel"} method) from \code{\link[stats]{density}}. Default 
#'   \code{"SJ"}.
#'   
#' @param adjust    Adjust bandwith parameter from density (for \code{"kernel"} 
#'   method) from \code{\link[stats]{density}}.
#'   
#' @param n   Number of equally spaced points at which the density (for 
#'   \code{"kernel"} method) is to be estimated (see 
#'   \code{\link[stats]{density}}).
#'   
#' @param mpoints Number of points in plot.
#'   
#' @param colors Colours of plot.
#'   
#' @param xlab Label of x-axis in plot.
#'   
#' @param ylab Label of y-axis in plot.
#'   
#' @param title Title of plot.
#'   
#' @param check_integrals Check that counterfactual of theoretical minimum risk 
#'   exposure and relative risk's expected values are well defined for this 
#'   scenario.
#'   
#' @param check_exposure  Check that exposure \code{X} is positive and numeric.
#'   
#' @param check_rr        Check that Relative Risk function \code{rr} equals 
#'   \code{1} when evaluated at \code{0}.
#'   
#' @param is_paf    Boolean forcing evaluation of \code{\link{paf}}
#'   
#' @return pif.plot       \code{\link[ggplot2]{ggplot}} object with plot of 
#'   Potential Impact Fraction as function of \code{theta}
#'   
#' @author Rodrigo Zepeda Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho García Formentí \email{daliaf172@gmail.com}
#'   
#' @import ggplot2
#'   
#' @examples 
#' 
#' #Example 1: Exponential Relative Risk empirical method
#' #-----------------------------------------------------
#' \dontrun{
#' set.seed(18427)
#' X <- rbeta(25, 4.2, 10)
#' pif.plot(X, thetalow = 0, thetaup = 10, rr =  function(X, theta){exp(theta*X)})
#' 
#' #Same example with kernel method
#' pif.plot(X, thetalow = 0, thetaup = 10, rr =  function(X, theta){exp(theta*X)}, method = "kernel",
#' title = "Kernel method example") 
#'  
#' #Same example for approximate method. Notice that approximate method has more uncertainty
#' Xmean <- mean(X)
#' Xvar  <- var(X)
#' pif.plot(Xmean, thetalow = 0, thetaup = 10, rr =  function(X, theta){exp(theta*X)}, 
#' method = "approximate", Xvar = Xvar, title = "Approximate method example")
#' 
#' #Example with counterfactual
#' pif.plot(X, thetalow = -10, thetaup = -5, rr = function(X, theta){exp(theta*X)}, 
#' cft = function(X){sqrt(X)})
#' 
#' #Example for approximate method with square root counterfactual
#' #Notice how the approximate represents information loss and thus the interval
#' #loses precision.
#' pif.plot(Xmean, thetalow = -10, thetaup = -5, rr = function(X, theta){exp(theta*X)},  
#' cft = function(X){sqrt(X)}, method = "approximate", Xvar = Xvar) 
#' }
#' 
#' @seealso \code{\link{pif}} for Potential Impact Fraction estimation with
#'   confidence intervals \code{\link{pif.confidence}}.
#'   
#'   \code{\link{paf.plot}} for same plot with Population Attributable Fraction.
#'   
#' @export

pif.plot <- function(X, thetalow, thetaup, rr,         
                     cft = NA,
                     weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                     method  = c("empirical", "kernel", "approximate"),
                     adjust = 1, n = 512, mpoints = 100,
                     Xvar    = var(X), 
                     deriv.method.args = list(), 
                     deriv.method      = c("Richardson", "complex"),
                     ktype  = c("gaussian", "epanechnikov", "rectangular", "triangular", 
                                "biweight","cosine", "optcosine"), 
                     bw     = c("SJ", "nrd0", "nrd", "ucv", "bcv"),
                     confidence_method = c("bootstrap", "linear", "loglinear"),
                     confidence = 95,
                     colors = c("deepskyblue", "gray25"), 
                     xlab = "Theta", ylab = "PIF",
                     title = "Potential Impact Fraction (PIF) under different values of theta",
                     check_exposure = TRUE, check_rr = TRUE, check_integrals = TRUE,
                     is_paf = FALSE){
  
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
    .dat           <- matrix(NA, nrow = mpoints, ncol = 4)
    colnames(.dat) <- c("Theta","Lower_CI", "Point_Estimate","Upper_CI")
    
    #Loop through values of theta for plot
    for (i in 1:mpoints){
      
      #Save theta value
      .dat[i,"Theta"] <- .theta[i]
      
      #Calculate PIF
      .dat[i,c("Lower_CI", "Point_Estimate","Upper_CI")]   <- 
        pif.confidence(X = X, thetahat = .theta[i], thetavar = 0, rr = rr, 
                       cft = cft, weights =  weights, 
                       method = method, Xvar = Xvar, 
                       deriv.method.args = deriv.method.args,
                       deriv.method = deriv.method,
                       confidence = confidence,
                       confidence_method = confidence_method,
                       adjust = adjust, n = n, ktype  = ktype, 
                       bw     = bw, check_exposure = check_exposure, check_rr = check_rr, 
                       check_integrals = check_integrals, is_paf = is_paf)[1:3] 
      
    }
    
    #Create plot
    .thetaplot <- ggplot(as.data.frame(.dat), aes(x = .dat[,"Theta"])) + 
                  geom_errorbar(aes(ymin = .dat[,"Lower_CI"], ymax = .dat[,"Upper_CI"], 
                                    color = "Pointwise Confidence")) +
                  geom_point(aes(y = .dat[,"Point_Estimate"], color = "Point Estimate")) +
                  xlab(xlab) + theme_bw() + ylab(ylab) +
                  ggtitle(title) + 
                  scale_color_manual(name = "", 
                                     values = c( "Point Estimate" = colors[1], 
                                                 "Pointwise Confidence" = colors[2])) 
    
    return(.thetaplot)

}
