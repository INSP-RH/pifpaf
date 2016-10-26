#' @title Plot of Potential Impact Fraction under different values of theta (univariate)
#' 
#' @description Function that plots the PIF under different values. 
#' 
#' @param X         Random sample (can be vector or matrix) which includes exposure and covariates.
#' 
#' @param thetamin  Minimum of theta for plot
#' 
#' @param thetamax  Maximum of theta for plot
#' 
#' @param rr        Function for relative risk
#' 
#' **Optional**
#' 
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for \code{PAF} where counterfactual is 0 exposure
#' 
#' @param weights   Survey \code{weights} for the random sample \code{X}
#' 
#' @param npoints   Number of points in plot (default 100).
#' 
#' @param method    Either \code{empirical} (default), \code{kernel} or {approximate}. 
#' 
#' @param Xvar      Variance of exposure levels.
#' 
#' @param ktype    \code{kernel} type from  \code{gaussian}, \code{epanechnikov}, \code{rectangular},
#'                  \code{triangular}, \code{biweight}, \code{cosine}, \code{optcosine}
#' 
#' @param bw        Smoothing bandwith parameter from density
#' 
#' @param adjust    Adjust bandwith parameter from density
#' 
#' @param mpoints   Number of points for kernel
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí
#' 
#' @import ggplot2
#' 
#' @examples 
#' 
#' #Example with risk given by HR
#' set.seed(18427)
#' library(ggplot2)
#' X <- rlnorm(100)
#' plotpif(X, 0, 2, function(X, theta){exp(theta*X)})+ ggtitle("PAF under different values of theta \n Empirical method")
#' 
#' #Example with kernel method
#' plotpif(X, 0, 2, function(X, theta){exp(theta*X)}, method = "kernel")+ ggtitle("PAF under different values of theta \n Kernel method")
#' 
#' #Example for approximate method
#' Xmean <- mean(X)
#' Xvar  <- var(X)
#' plotpif(X, 0, 2, function(X, theta){exp(theta*X)}, method = "approximate", Xvar = Xvar)+ ggtitle("PAF under different values of theta \n Approximate method")
#' 
#' #Example with square root counterfactual
#' plotpif(X, 0, 2, function(X, theta){exp(theta*X)}, cft = function(X){sqrt(X)})+ ggtitle("PIF under different values of theta \n square root counterfactual  \n Empirical method")
#' 
#' #Example for approximate method with square root counterfactual
#' plotpif(X, 0, 2, function(X, theta){exp(theta*X)},  cft = function(X){sqrt(X)}, method = "approximate", Xvar = Xvar) + ggtitle("PIF under different values of theta \n square root counterfactual  \n Approximate method")
#' 
#' @export

plotpif <- function(X, thetamin, thetamax, rr, 
                    cft = function(Varx){matrix(0,ncol = ncol(as.matrix(Varx)), nrow = nrow(as.matrix(Varx)))}, 
                    weights = rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), npoints = 100, method = c("empirical", "kernel", "approximate"),
                    Xvar = 1,
                    ktype = "epanechnikov", bw = "nrd0", adjust = 1, mpoints = 1000){
  
  #Check that we are able to plot
  if (thetamin < thetamax){
    
    #Create sequence from thetamin to thetamax
    .theta <- seq(thetamin, thetamax, length.out = ceiling(npoints))
    
    #Create data frame for saving values of theta
    .dat   <- matrix(NA, nrow = npoints, ncol = 2)
    colnames(.dat) <- c("Theta","PIF")
    
    #Loop through values of theta for plot
    for (i in 1:npoints){
      
      #Save theta value
      .dat[i,"Theta"] <- .theta[i]
      
      #Calculate PIF
      .dat[i,"PIF"]   <- pif(X, .theta[i], rr, cft, weights, method = method,
                             Xvar = Xvar,
                             ktype = ktype, bw = bw, adjust = adjust, npoints = mpoints) 
      
    }
    
    #Create plot
    return(
      ggplot(as.data.frame(.dat)) + 
      geom_path(aes(x = .dat[,"Theta"], y = .dat[,"PIF"]), color = "darkslategrey") +
      xlab("Theta") +
      theme_bw() + 
      ylab("PIF") +
      ggtitle("Potential Impact Fraction (PIF) under different values of theta")
    )
    
  } else {
    
    warning("Minimum thetamin  cannot be equal or greater than maximum thetamax")
    return(NA)
    
  }

}
