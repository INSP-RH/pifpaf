#' @title Potential Impact Fraction
#' 
#' @description Function that calculates the potential impact fraction via the empirical method
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
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for \code{PAF} where counterfactual is 0 exposure
#' 
#' @param weights   Survey \code{weights} for the random sample \code{X}
#' 
#' @param eval.cvx  Boolean for checking if relative risk function passes convexity or concavity test.
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
#' @param maxval Number of values to sample to discard a function as not being convex. 
#' 
#' @param tol    Tolerance of the method (\code{var} assumed \code{0} if \code{var < tol})              
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' 
#' @examples 
#' 
#' #Example with risk given by HR
#' set.seed(18427)
#' X <- rlnorm(100)
#' thetahat <- 0.12
#' pif(X, thetahat, function(X, theta){exp(theta*X)})
#' 
#' #Same example considering counterfactual of halfing exposure
#' pif(X, thetahat, function(X, theta){exp(theta*X)}, cft = function(X){ 0.5*X })
#' 
#' 
#' #Example with linear regression
#' 
#' 
#' #Example with  function of risk that is not convex nor concave
#' set.seed(523)
#' X <- rgamma(100, shape = 1)
#' thetahat <- 1.81
#' pif(X, thetahat, function(X, theta){X*(theta-1.82)^3},  eval.cvx = TRUE)
#' 
#' #Notice that deactivating the convexity checker returns results but they might be incorrect
#' pif(X, thetahat, function(X, theta){X*(theta-1.82)^3}, eval.cvx = FALSE)
#' 
#' #Notice that specifying theta's domain limits the revision of convexity to the area
#' #where the function actually is convex
#' pif(X, thetahat, function(X, theta){X*(theta-1.82)^3}, eval.cvx = TRUE, tmin = 1, tmax =1.819)
#' 
#' 
#' @export



pif <- function(X, thetahat, rr, 
                cft = function(Varx){matrix(0,ncol = ncol(as.matrix(Varx)), nrow = nrow(as.matrix(Varx)))}, 
                weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), eval.cvx = TRUE, xmin = NA, 
                xmax = NA, tmin = NA, tmax = NA, maxval = 1000, tol = 1.e-8){
  
  #Set X as matrix
  X <- as.matrix(X)
  
  #Check that function is convex
  if (eval.cvx){ 
    
    #Get minima and maxima if they do not exist for X
    if (is.na(xmin)) { xmin <- 0.99*colMins(X) }
    if (is.na(xmax)) { xmax <- 1.01*colMaxs(X) }  
    
    #Get minima and maxima for t if they do not exist
    if (is.na(tmax)) { tmax <- 1.25*thetahat   }
    if (is.na(tmin)) { tmin <- 0.75*thetahat   }
    
    #Check if function is convex
    cvx   <- ismulticonvex(function(X, theta){ rr(X, theta)}, xmin, xmax, tmin, tmax, maxval, tol)
    
    #Check if function is concave
    if (!cvx){
      ccv <- ismulticonvex(function(X, theta){-rr(X, theta)}, xmin, xmax, tmin, tmax, maxval, tol)
    } else {
      ccv <- FALSE
    }
    
    if(!cvx & !ccv){
      warning(paste0("Function for Relative Risks RR(X;theta) does not appear convex nor concave in theta",
                     "If you are 100 % sure your function is convex please make eval.cvx = FALSE"))
      stop()
    }
  
  #If user disables option of checking, assume function is convex.
  } else {
    
    cvx <- TRUE
    ccv <- TRUE
    
  }
  
  #Evaluate if function is convex or concave
  if (cvx || ccv){
    
    #Estimate weighted sums
    mux   <- sum(weights*rr(X,thetahat))
    mucft <- sum(weights*rr(cft(X),thetahat))
    #mux <- median(rr(X,thetahat))
    #mucft <- median(rr(cft(X),thetahat))
    
    #Calculate PIF
    pif   <- 1 - mucft/mux
    
  } else {
    
    pif <- NA
    
  }
  
  return(pif)
  
}


