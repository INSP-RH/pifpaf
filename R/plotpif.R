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
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' 
#' @import ggplot2
#' 
#' @examples 
#' 
#' #Example with risk given by HR
#' set.seed(18427)
#' X <- rlnorm(100)
#' plotpif(X, 0, 2, function(X, theta){exp(theta*X)})
#' 
#' #Example with square root counterfactual
#' plotpif(X, 0, 2, function(X, theta){exp(theta*X)}, cft = function(X){sqrt(X)})
#' 
#' @export

plotpif <- function(X, thetamin, thetamax, rr, 
                    cft = function(Varx){matrix(0,ncol = ncol(as.matrix(Varx)), nrow = nrow(as.matrix(Varx)))}, 
                    weights = rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), npoints = 100){
  
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
      .dat[i,"PIF"]   <- pif(X, .theta[i], rr, cft, weights) 
      
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
