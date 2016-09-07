#' @title Multivariate convexity checker
#' 
#' @description Checks a function \code{f} for any indication of it not being convex. 
#'              **Note that passing this test does not imply convexity.**
#'  
#' @param f     Function in which convexity will be checked. Note that \code{f} must 
#'              be a function of two vectors \code{theta} and \code{X} with \code{f = f(X,theta)}
#'              
#' @param xmin  Minimum plaussible value of each entry of \code{X} can be found. That is: 
#'              \code{xmin[i]} <= \code{X[i]} <= \code{xmax[i]}
#' 
#' @param xmax  Maximum plaussible value of each entry of \code{X} can be found. That is: 
#'              \code{xmin[i]} <= \code{X[i]} <= \code{xmax[i]}
#' 
#' @param tmin  Minimum plaussible value of \code{theta} can be found. That is: 
#'              \code{tmin[i]} <= \code{theta[i]} <= \code{tmax[i]}
#'              
#' @param tmax  Maximum plaussible value of \code{theta} can be found. That is: 
#'              \code{tmin[i]} <= \code{theta[i]} <= \code{tmax[i]}           
#'              
#' **Optional**
#'              
#' @param maxval Number of values to sample to discard a function as not being convex. 
#' 
#' @param tol    Tolerance of the method (\code{var} assumed \code{0} if \code{var < tol})
#' 
#' @importFrom stats runif 
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' 
#' @examples
#' #Check function is convex
#' ismulticonvex(function(X,theta){X*exp(theta[1] + 0.12*theta[2])}, 0, 1, c(0,0.12), c(0.1, 1))
#' 
#' #Check function is convex
#' ismulticonvex(function(X,theta){sin(X[1]^(theta[1])) *
#'                 cos(X[2]*theta[2])}, c(0,0), c(1,1), c(0,0), c(1, 1))
#' 
#' #Check function is concave. (Recall function f is concave if -f is convex)
#' ismulticonvex(function(X,theta){- sqrt(X[1]*X[2]*X[3]*theta[1])}, c(0,0,0), c(1,1,1), 1, 2)
#' 
#' @export

ismulticonvex <- function(f, xmin, xmax, tmin, tmax, maxval = 1000, tol = 1.e-8){
  
  #Check that they represent intervals with min < max
  if ( !all(xmin < xmax) || !all(tmin < tmax) || 
       length(xmin) != length(xmax) || length(tmin) != length(tmax)){
    
    warning("Either x values and/or t values are incorrect. Check that xmin < xmax and tmin < tmax for all entries.")
    stop()
    
  } else if (maxval < 1) {
    
    warning("Not enough maxval values for simulation")
    stop()
    
  } else {
    
    #Initial conditions for loop
    n     <- 1       # number of loop
    cvx   <- TRUE    # convexity conclusion
    
    #Check maxval is integer
    maxval <- ceiling(maxval)
    
    #Get n alphas < 1 for checking convexity
    alpha  <- runif(maxval)
    
    #Get length of theta and length of x
    xlen   <- length(xmin)
    tlen   <- length(tmin)
    
    #Create the vectors
    xval   <- vector("numeric", length = xlen)
    theta1 <- vector("numeric", length = tlen)
    theta2 <- vector("numeric", length = tlen)
    
    #Loop through each value of maxval looking for a counterexample
    while( n <= maxval & cvx){
      
      #Simulate random x and t vectors
      for (i in 1:xlen){
        xval[i] <- runif(1, xmin[i], xmax[i])
      }
      
      for (i in 1:tlen){
        theta1[i] <- runif(1, tmin[i], tmax[i])
        theta2[i] <- runif(1, tmin[i], tmax[i])
      }
      
      #Check convexity condition
      cvx <- (alpha[n] * f(xval, theta1) + (1-alpha[n]) * f(xval, theta2) - 
                f(xval, alpha[n] * theta1 + (1-alpha[n]) * theta2 )) > - tol
      
      #Update n
      n   <- n + 1

    } #Close while
    
    return(cvx)
    
  } #Close else
} #Close function