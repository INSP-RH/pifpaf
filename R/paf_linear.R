#' @title Population Attributable Fraction with Linear Relative Risk Function
#'   
#' @description Function that calculates the Population Attributable Fraction
#'   \code{\link{paf}} with linear Relative Risk function \code{rr} given by 
#'   \deqn{
#'   rr(X; \theta) = \theta_1 + \sum\limits_{i=1}^{n} \theta_{i+1} X_i.
#'   }{
#'   rr(X, theta) = theta[1] + theta[2]*X[,1] + theta[3]*X[,2] + ... + theta[n+1]*X[,n].
#'   }
#'   
#' @param X         Random sample (\code{data.frame}) which includes exposure 
#'   and covariates or sample \code{mean} if \code{"approximate"} method is 
#'   selected.
#'   
#' @param thetahat  Asymptotically consistent or Fisher consistent estimator (\code{vector}) of \code{theta} for the Relative 
#'   Risk function \code{rr}.
#'   
#'   \strong{**Optional**}
#'   
#' @param weights   Normalized survey \code{weights} for the sample \code{X}.
#'   
#' @param method    Either \code{"empirical"} (default), \code{"kernel"} or 
#'   \code{"approximate"}. For details on estimation methods see 
#'   \code{\link{pif}}.
#'   
#' @param Xvar      Variance of exposure levels (for \code{"approximate"} 
#'   method).
#'   
#' @param deriv.method.args \code{method.args} for 
#'   \code{\link[numDeriv]{hessian}} (for \code{"approximate"} method).
#'   
#' @param deriv.method      \code{method} for \code{\link[numDeriv]{hessian}}. 
#'   Don't change this unless you know what you are doing (for 
#'   \code{"approximate"} method).
#'   
#' @param ktype    \code{kernel} type:  \code{"gaussian"}, 
#'   \code{"epanechnikov"}, \code{"rectangular"}, \code{"triangular"}, 
#'   \code{"biweight"}, \code{"cosine"}, \code{"optcosine"} (for \code{"kernel"}
#'   method). Additional information on kernels in \code{\link[stats]{density}}.
#'   
#' @param bw        Smoothing bandwith parameter (for 
#'   \code{"kernel"} method) from \code{\link[stats]{density}}. Default 
#'   \code{"SJ"}.
#'   
#' @param adjust    Adjust bandwith parameter (for \code{"kernel"} 
#'   method) from \code{\link[stats]{density}}.
#'   
#' @param n   Number of equally spaced points at which the density (for 
#'   \code{"kernel"} method) is to be estimated (see 
#'   \code{\link[stats]{density}}).
#'   
#' @param check_integrals \code{boolean}  Check that counterfactual \code{cft} 
#'   and relative risk's \code{rr} expected values are well defined for this 
#'   scenario.
#'   
#' @param check_exposure  \code{boolean}  Check that exposure \code{X} is 
#'   positive and numeric.
#'   
#' @param check_rr        \code{boolean} Check that Relative Risk function
#'   \code{rr} equals \code{1} when evaluated at \code{0}.
#'   
#' @return paf      Estimate of Population Attributable Fraction with linear
#'   relative risk.
#'   
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'   
#' @note The \code{"approximate"} method should be the last choice. In practice 
#'   \code{"empirical"} should be preferred as convergence is faster than 
#'   \code{"kernel"} for most functions. In addition, the scope of
#'   \code{"kernel"} is limited as it does not work with multivariate exposure
#'   data \code{X}.
#'   
#' @note \code{\link{paf.linear}} is a wrapper for \code{\link{paf}} with linear
#'   relative risk.
#'   
#'   
#' @examples 
#' 
#' #Example 1: Univariate relative risk
#' #----------------------------------------
#' set.seed(18427)
#' X <- data.frame(Exposure = rnorm(100,3,.5))
#' thetahat <- c(1, 0.12)   #Linear risk given by 1 + 0.12*X
#' paf.linear(X, thetahat)
#' 
#' #This is the same as doing:
#' paf(X, thetahat, rr = function(X, theta){X*theta[2] + theta[1]})
#'
#' #Same example with kernel method
#' paf.linear(X, thetahat,  method = "kernel")
#' 
#' #Same example with approximate method
#' Xmean <- data.frame(mean(X[,"Exposure"]))
#' Xvar  <- var(X)
#' paf.linear(Xmean, thetahat, method = "approximate", Xvar = Xvar)
#' 
#' #Example 2: Multivariate relative risk
#' #----------------------------------------
#' X     <- data.frame(Exposure = rnorm(100,2,.7), Covariate = rnorm(100,4,1))
#' theta <- c(1, 0.3,0.1)
#' paf.linear(X, theta)   #Linear risk given by 1 + 0.3*X1 + 0.1*X2 
#' 
#' #Example 3: Polynomial relative risk
#' #----------------------------------------
#' X     <- runif(100)
#' X2    <- X^2
#' X3    <- X^3
#' matX  <- data.frame(X,X2,X3)
#' theta <- c(1, 0.3,0.1, 0.4)
#' paf.linear(matX,theta) #Polynomial risk: 1 + 0.3*X + 0.1*X^2 + 0.4*X^3
#' 
#' @seealso  
#' 
#' See \code{\link{paf}} for Population Attributable Fraction (with
#'   arbitrary relative risk), and \code{\link{pif}} for Potential Impact Fraction
#'   estimation.
#'   
#'   See \code{\link{paf.exponential}} for PAF with ready-to-use exponential
#'   relative risk function.
#'   
#'   For more information on kernels see \code{\link[stats]{density}}.
#'   
#' @references Vander Hoorn, S., Ezzati, M., Rodgers, A., Lopez, A. D., & 
#'   Murray, C. J. (2004). \emph{Estimating attributable burden of disease from 
#'   exposure and hazard data. Comparative quantification of health risks: 
#'   global and regional burden of disease attributable to selected major risk 
#'   factors}. Geneva: World Health Organization, 2129-40.
#'   
#' @export

paf.linear <- function(X, thetahat,
                       method = "empirical",
                       weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                       Xvar    = var(X), 
                       deriv.method.args = list(), 
                       deriv.method      = c("Richardson", "complex"),
                       adjust = 1, n = 512,
                       ktype  = "gaussian", 
                       bw     = "SJ",
                       check_exposure = TRUE, check_rr = TRUE, check_integrals = TRUE){
  
  #Convert exposure to matrix
  .X    <- as.data.frame(X)
  
  #Check that there are parameters for every covariate
  if(ncol(.X) != (length(thetahat) - 1)){
    stop(paste0("The amount of parameters in theta must be equal to the number ", 
                "of exposure values and covariates + 1 in each observation"))
  }
  
  #Create function for linear relative risk
  .rr <- function(.myX, .mytheta){
    
    #Sum everyone
    .expsol <- .mytheta[1]
    for (.i in 1:ncol(.myX)){
      .expsol <- .expsol + .mytheta[.i + 1]*.myX[,.i]
    }
    
    return(.expsol)
  }
  
  #Estimate Population attributable fraction
  .paf <- paf(X = X, thetahat = thetahat,   rr = .rr,         
              method  = method, weights = weights, 
              Xvar    = Xvar, deriv.method.args = deriv.method.args, 
              deriv.method = deriv.method, adjust = adjust, n = n,
              ktype  = ktype, bw     = bw,
              check_exposure = check_exposure, check_rr = check_rr,
              check_integrals = check_integrals)
  
  return(.paf)
}


