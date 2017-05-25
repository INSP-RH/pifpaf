#' @title Pivotal Boostrap Confidence Intervals for Potential Impact Fraction
#'   
#' @description  Estimates a 1 - alpha pivotal confidence interval for the
#'   potential impact fraction \code{pif} using a boostrap approximation.
#'   
#' @param X         Random sample (vector or matrix) which includes exposure and
#'   covariates.
#'   
#' @param thetahat  Maximum Likelihood estimator (vector or matrix) of
#'   \code{theta} for the Relative Risk function.
#'   
#' @param thetavar   Estimator of variance of \code{thetahat}
#'   
#' @param rr        Function for Relative Risk which uses parameter 
#'   \code{theta}. The order of the parameters shound be \code{rr(X, theta)}.
#'   
#'   
#'   **Optional**
#'   
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'   the Population Attributable Fraction \code{\link{paf}} where counterfactual
#'   is 0 exposure.
#'   
#' @param weights   Normalized survey \code{weights} for the sample \code{X}.
#'   
#' @param method    Either \code{empirical} (default), or \code{kernel}
#'   
#' @param nboost    Number of samples in Bootstrap
#'   
#' @param confidence Concidence level (0 to 100) default = \code{95} \%
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
#' @param check_thetas Checks that theta parameters are correctly inputed
#'   
#' @param check_integrals Check that counterfactual and relative risk's expected
#'   values are well defined for this scenario
#'   
#' @param check_exposure  Check that exposure \code{X} is positive and numeric
#'   
#' @param check_rr        Check that Relative Risk function \code{rr} equals 
#'   \code{1} when evaluated at \code{0}
#'   
#' @param is_paf   Boolean forcing evaluation of \code{paf}
#'   
#' @return pif      Estimate of Potential Impact Fraction
#'   
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'  
#'   
#' @importFrom MASS mvrnorm
#' @importFrom stats var quantile
#'   
#' @examples 
#' #Example 1: Exponential Relative Risk
#' #----------------------------------------
#' set.seed(18427)
#' X        <- rnorm(100,5,1)
#' thetahat <- 0.4
#' thetavar <- 0.1
#' pif.confidence.bootstrap(X, thetahat, thetavar, function(X, theta){exp(theta*X)}, 
#'                          nboost = 100) #nboost small only for example purposes
#' 
#' #This also works with kernel method
#' pif.confidence.bootstrap(X, thetahat, thetavar, function(X, theta){exp(theta*X)}, 
#'                          nboost = 100, method = "kernel") 
#' 
#' #Example 2: Multivariate example
#' #----------------------------------------
#' \dontrun{
#' set.seed(18427)
#' X1 <- rnorm(100, 1, 0.05)
#' X2 <- rnorm(100, 1, 0.05)
#' X  <- as.matrix(cbind(X1,X2))
#' thetahat <- c(2, 0.03)
#' thetavar <- matrix(c(0.1, 0, 0, 0.05), byrow = TRUE, nrow = 2)
#' rr        <- function(X, theta){
#'   .X <- as.matrix(X, ncol = 2)
#'   exp(theta[1]*.X[,1] + theta[2]*.X[,2])
#' }
#' cft <- function(X){0.5*X}#' cft <- function(X){0.95*X}
#' pif.confidence.bootstrap(X, thetahat, thetavar, rr, cft) 
#' }
#' @keywords internal
#' @export

pif.confidence.bootstrap <- function(X, thetahat, thetavar, rr,         
                                     cft = NA,
                                     weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                                     method  = c("empirical", "kernel"),
                                     nboost = 10000, 
                                     adjust = 1, n = 512,
                                     confidence = 95,
                                     ktype  = c("gaussian", "epanechnikov", "rectangular", "triangular", 
                                                "biweight","cosine", "optcosine"), 
                                     bw     = c("SJ", "nrd0", "nrd", "ucv", "bcv"),
                                     check_exposure = TRUE, check_rr = TRUE, check_integrals = TRUE,
                                     check_thetas = TRUE, is_paf = FALSE){
  
  
  #Get confidence
  check.confidence(confidence)
  .alpha <- max(0, 1 - confidence/100)
  
  #Check that nboost is integer
  .nboost <- max(10, ceiling(nboost))
  
  #Check 
  .thetavar <- as.matrix(thetavar)
  if(check_thetas){ check.thetas(.thetavar, thetahat, NA, NA, "linear") }
  
  #Save X as matrix
  .X       <- as.data.frame(X)
  
  #Get original pif
  .pif     <- pif(X = .X, thetahat = thetahat, rr = rr, cft = cft, method = method, weights = weights, 
                  adjust = adjust, n = n, ktype = ktype, bw = bw,
                  check_exposure = check_exposure, check_rr = check_rr, check_integrals = check_integrals,
                  is_paf = is_paf) 
  
  #Boostrap pif will be saved here
  .bpif   <- rep(NA, .nboost)
  
  #Get vector with numbers indicating which X's to sample
  .index  <- 1:nrow(.X)
  
  #Get new thetas using asymptotic normality
  .newtheta <- mvrnorm(.nboost, thetahat, .thetavar, empirical = TRUE)
  
  #Loop through nboost simulations resampling according to probability normalized weights
  for (i in 1:.nboost){
    
    #Get resampled values
    .whichX  <- sample(.index, nrow(.X), replace = TRUE, prob = weights) #Resampling 
    .Xboost  <- as.data.frame(.X[.whichX, ])                                           #Rows = individuals
    colnames(.Xboost) <- colnames(.X)
    .wboost  <- weights[.whichX]/sum(weights[.whichX])                   #Renormalize weights ????
    
    #Get pif
    .bpif[i] <- pif(X = .Xboost, thetahat = .newtheta[i,], rr = rr, cft = cft, 
                    method = method, weights = .wboost,  Xvar = NA,
                    adjust = adjust, n = n, ktype = ktype, bw = bw, 
                    check_exposure = FALSE, check_rr = FALSE,   #Set FALSE to avoid recalculating 
                    check_integrals = FALSE, is_paf = is_paf)                             
  }
  
  #Get the confidence interval and variance
  if (any(is.na(.bpif))){
    .quantiles        <- c(NaN, NaN)
  } else {
    .quantiles        <- quantile(.bpif, probs = c(1 - .alpha/2, .alpha/2))  
  }
  
  names(.quantiles) <- c()        #Delete names
  .variance         <- var(.bpif)
  
  #Create vector of ci
  .ci <- c("Lower_CI" = 2*.pif - .quantiles[1], "Point_Estimate"     = .pif, 
          "Upper_CI" =  2*.pif - .quantiles[2], "Estimated_Variance" = .variance)
  
  .ci["Upper_CI"] <- min(.ci["Upper_CI"], 1)
  
  #Return ci
  return(.ci)
}