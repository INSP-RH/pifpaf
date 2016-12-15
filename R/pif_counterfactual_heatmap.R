#' @title  Sensitivity Analysis of PIF with counterfactual
#' 
#' @description Provides a graphic sensitivity analysis by varying the parameters of a bivariate counterfactual
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
#' @param weights   Survey \code{weights} for the random sample \code{X}
#' 
#' @param cft       Counterfactual function with parameters \code{a} and \code{b}. Default \code{aX + b}
#' 
#' @param mina      Minimum for \code{a} for the counterfactual 
#' 
#' @param minb      Minimum for \code{b} for the counterfactual 
#' 
#' @param maxa      Maximum for \code{a} for the counterfactual 
#' 
#' @param maxb      Maximum for \code{b} for the counterfactual 
#' 
#' @param nmesh     Number of points in mesh (default \code{30})
#' 
#' @param title     Title for the plot
#' 
#' @param xlab      Label for the X-axis of the plot (corresponding to "a")
#' 
#' @param ylab      Label for the Y-axis of the plot (corresponding to "b")
#' 
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí
#' 
#' @examples 
#' 
#' #Example 1
#' X  <- rnorm(25,3)                        #Sample
#' rr <- function(X,theta){exp(X*theta)}    #Relative risk
#' theta <- 0.01                            #Estimate of theta
#' pif.counterfactual.heatmap(X, theta = theta, rr = rr)
#' 
#' #Example 2
#' X     <- rbeta(100, 1, 0.2)
#' theta <- c(0.12, 1)
#' rr    <- function(X,theta){X*theta[1] + theta[2]}
#' cft   <- function(X, a, b){sin(a*X)*b}
#' pif.counterfactual.heatmap(X, theta = theta, rr = rr, cft = cft, 
#'      nmesh = 15, palette = rainbow(30), method = "kernel",
#'      title = "PIF with counterfactual cft(X) = sin(a*X)*b")
#' 
#' #You can also plot univariate counterfactuals
#' X     <- rgamma(100, 1, 0.2)
#' theta <- c(0.12, 1)
#' rr    <- function(X,theta){X*theta[1] + theta[2]}
#' cft   <- function(X, a, b){sqrt(a*X)}   #Leave two variables in it
#' pif.counterfactual.heatmap(X, theta = theta, rr = rr, 
#' cft = cft, minb = 0, maxb = 0, title ="Univariate counterfactual", ylab = "")
#' 
#' @import ggplot2
#' @export


pif.counterfactual.heatmap <-function(X, thetahat, rr, 
                                 weights = rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                                 method = "empirical", legendtitle = "PIF",
                                 mina = 0.01, maxa = 0.99, minb = -1, maxb = 0, nmesh = 10,
                                 title = "Potential Impact Fraction (PIF) with counterfactual \n f(X)= aX+b",
                                 xlab = "a", ylab = "b", cft = function(X, a, b){a*X + b},
                                 palette = heat.colors(nmesh)){
  
  #Check that mina < maxa and minb <= maxb
  if (mina >= maxa || minb > maxb){stop("Error: min values of a or b larger than max values")}
  
  #Create a mesh of all possible a and b values
  M <- expand.grid(a   = seq(mina, maxa, length.out = nmesh), 
                   b   = seq(minb, maxb, length.out = nmesh),
                   pif = rep(NA, nmesh))
  
  #Calculate the PIF for the specific counterfacvtual
  for(i in 1:nrow(M)){
     M[i,3] <- pif(X, thetahat = thetahat, rr = rr, weights = weights, 
                     cft = function(X){cft(X, M$a[i], M$b[i])} )
  }
  
  #Create Heatmap
  plotobject <- ggplot(M, aes(x = a, y = b, fill = pif)) + geom_tile() + 
    xlab(xlab) + ylab(ylab) + ggtitle(title) + theme_classic() + 
    scale_fill_gradientn(legendtitle, colours = palette)
  
  #Check that minb == maxb to delete its axis
  if(minb == maxb){
    plotobject <- plotobject + theme(axis.text.y = element_blank(), 
                                     axis.ticks.y = element_blank(),
                                     axis.line.y = element_blank())
  }
  
  return(plotobject)
}

