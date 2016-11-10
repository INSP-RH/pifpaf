#' @title 3D Sensitivity Analysis of PIF with affine counterfactual
#' 
#' @description Provides a graphic sensitivity analysis by varying the parameters of the counterfactual \code{a*X + b}
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
#' @param mina      Minimum for \code{a} for the counterfactual 
#' 
#' @param minb      Minimum for \code{b} for the counterfactual 
#' 
#' @param maxa      Maximum for \code{a} for the counterfactual 
#' 
#' @param maxb      Maximum for \code{b} for the counterfactual 
#' 
#' @param title     Title for the plot
#' 
#' @param xlab      Label for the X-axis of the plot (corresponding to "a")
#' 
#' @param ylab      Label for the Y-axis of the plot (corresponding to "b")
#' 
#' @param zlab      Label for the Z-axis of the plot (corresponding to "PIF")
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
#' pif.counterfactual.3D(X, theta = theta, rr = rr)
#' 
#' #Example 2
#' X     <- rbeta(100, 1, 0.2)
#' theta <- c(0.12, 1)
#' rr    <- function(X,theta){X*theta[1] + theta[2]}
#' pif.counterfactual.3D(X, theta = theta, rr = rr)
#' 
#' 
#' @import plot3D
#' @export


pif.counterfactual.3D <-function(X, thetahat, rr, 
                                 weights = rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                                 mina = 0, maxa = 1, minb = -1, maxb = 0, 
                                 title = "Potential Impact Fraction (PIF) with counterfactual \n f(X)= aX+b",
                                 xlab = "a", ylab = "b", zlab ="PIF"){
  
  #Check limits were correctly specified
  if( mina >= maxa){
    warning("The minimum value of a is greater or equal than the maximum value of a, verify the values are correct.")
  }
  
  if( minb >= maxb){
    warning("The minimum value of b is greater or equal than the maximum value of b, verify the values are correct.")
  }
  
  #Create a mesh of all possible a and b values
  M <- mesh(seq(mina, maxa, length.out = 30), seq(minb, maxb, length.out = 30))
  a <- M$x
  b <- M$y
  z <- a+b
  
  #Calculate the PIF for the specific counterfacvtual
  for(i in 1:30){
    for(j in 1:30){
      z[i,j] <- pif(X,thetahat = thetahat, rr = rr, weights = weights, cft = function(X){a[i,i]*X+b[j,j]})
    }
  }
  
  #Create a 3D surface
  surf3D(x = a, y = b, z = z, bty="b2", main = title, xlab = xlab, ylab = ylab, zlab = zlab)
}

