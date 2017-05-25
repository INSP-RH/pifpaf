#' @title Confidence Intervals for the Population Attributable Fraction
#'   
#' @description Function that estimates confidence intervals for the Population 
#'   Attributable Fraction \code{\link{paf}} from a cross-sectional sample of 
#'   the exposure \code{X} with a known Relative Risk function \code{rr} with meta-analytical
#'   parameter \code{theta}, where the Population Attributable Fraction is given
#'   by: \deqn{ PAF = 
#'   \frac{E_X\left[rr(X;\theta)\right]-1}{E_X\left[rr(X;\theta)\right]}. }{ PAF 
#'   = mean(rr(X; theta) - 1)/mean(rr(X; theta)).}
#'   
#' @param X         Random sample (\code{data.frame}) which includes exposure 
#'   and covariates or sample \code{mean} if \code{"approximate"} method is 
#'   selected.
#'   
#' @param thetahat  Asymptotically consistent of Fisher consistent estimator (\code{vector})
#'  of \code{theta} for the Relative Risk function. \code{thetahat} should be 
#'  asymptotically normal with mean \code{theta} and variance \code{var_of_theta}.
#'   
#' @param rr        \code{function} for Relative Risk which uses parameter 
#'   \code{theta}. The order of the parameters should be \code{rr(X, theta)}.
#'   
#'   
#'  \strong{ **Optional**}
#'   
#' @param thetavar   Estimator of variance \code{var_of_theta} of asymptotic 
#'   normality of \code{thetahat}.
#'   
#' @param thetalow  (\code{vector}) lower bound of the confidence interval of 
#'   \code{theta}.
#'   
#' @param thetaup   (\code{vector}) upper bound of the confidence interval of 
#'   \code{theta}.
#'   
#' @param weights   Normalized survey \code{weights} for the sample \code{X}.
#'   
#' @param nsim      Number of simulations for estimation of variance.
#'   
#' @param confidence Confidence level \% (default \code{95}). If 
#'   \code{confidence_method} \code{"one2one"} is selected, \code{confidence} 
#'   should be at most the one from \code{theta}'s confidence interval 
#'   (\code{confidence_theta}\%).
#'   
#' @param confidence_method  Either \code{bootstrap} (default) \code{inverse}, 
#'   \code{one2one}, \code{linear}, \code{loglinear}. See details for additional
#'   explanation.
#'   
#' @param confidence_theta Confidence level \% of \code{theta} corresponding to
#' the interval [\code{thetalow}, \code{thetaup}] (default: \code{99}\%).
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
#' @param check_thetas \code{boolean} Check that theta associated parameters are
#'   correctly inputed for the model.
#'   
#' @param check_exposure  \code{boolean}  Check that exposure \code{X} is 
#'   positive and numeric.
#'   
#' @param check_cft  \code{boolean}  Check that counterfactual function 
#'   \code{cft} reduces exposure.
#'   
#' @param check_xvar \code{boolean} Check \code{Xvar} is a covariance matrix.
#'   
#' @param check_integrals \code{boolean}  Check that counterfactual \code{cft} 
#'   and relative risk's \code{rr} expected values are well defined for this 
#'   scenario.
#'   
#' @param check_rr         \code{boolean} Check that Relative Risk function 
#' \code{rr} equals \code{1} when evaluated at \code{0}.
#'   
#' @param force.min Boolean indicating whether to force the \code{rr} to have a 
#'   minimum value of 1 instead of 0 (not recommended). This works only for 
#'   \code{confidence_method} \code{"inverse"}.
#'   
#' @return pafvec Vector with lower (\code{"Lower_CI"}), and upper 
#'   (\code{"Upper_CI"}) confidence bounds for the \code{\link{paf}} as well as
#'   point estimate \code{"Point_Estimate"} and estimated variance or variance
#'   of \code{log(paf)} (if \code{confidence_method} is \code{"loglinear"}).
#'   
#' @note \code{\link{paf.confidence}} is a wrapper for
#'   \code{\link{pif.confidence}} with counterfactual of theoretical
#'   minimum risk exposure (\code{rr = 1}) .
#'   
#' @note For more information on kernels see \code{\link[stats]{density}}.
#'   
#' @note Do not use the \code{$} operator when using \code{"approximate"}
#'   \code{method}.
#'   
#' @details The \code{confidence_method} estimates confidence intervals with
#' different methods. A bootstrap approximation is conducted by 
#' \code{"bootstrap"}. The Delta Method is applied to \code{\link{paf}}
#' or \code{log(paf)} when choosing \code{"linear"} and \code{"loglinear"}
#' respectively. The \code{"inverse"} method estimates confidence intervals
#' for the Relative Risk function \code{rr} and applies the transformation
#' \code{1 - 1/rr}. Finally, \code{"one2one"} works with functions for which
#' the expected value over \code{X} of the relative risk is injective in 
#' \code{theta}. 
#' 
#' Additional information on confidence method estimations can be found
#' in the package's vignette: \code{browseVignettes("pifpaf")}.
#'   
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'  
#' 
#' @seealso  \code{\link{pif.confidence}} for confidence interval estimation of
#'   \code{\link{pif}}, and  \code{\link{paf}} for only point estimates.
#'   
#'   Sensitivity analysis plots can be done with \code{\link{paf.plot}}, and 
#'   \code{\link{paf.sensitivity}}.
#'   
#' @examples 
#' 
#' #Example 1: Exponential Relative Risk
#' #--------------------------------------------
#' set.seed(18427)
#' X        <- data.frame(Exposure = rnorm(100,3,1))
#' thetahat <- 0.32
#' thetavar <- 0.02
#' rr       <- function(X, theta){exp(theta*X)}
#' 
#' #Using bootstrap method
#' paf.confidence(X, thetahat, rr, thetavar)
#' 
#' \dontrun{
#' #Same example with loglinear method
#' paf.confidence(X, thetahat, rr, thetavar, confidence_method = "loglinear")
#' 
#' #Same example with linear method (usually the widest and least precise)
#' paf.confidence(X, thetahat, rr, thetavar, confidence_method = "linear")
#' 
#' #Same example with inverse method 
#' paf.confidence(X, thetahat, rr, thetavar, confidence_method = "inverse")
#' 
#' #Same example with one2one method 
#' #assume 99% ci of theta is [0.27, 0.35]
#' paf.confidence(X, thetahat, rr, thetalow = 0.27, thetaup = 0.35, 
#' confidence_method = "one2one", confidence_theta = 99)
#' 
#' #Example 2: Linear Relative Risk with weighted sample
#' #--------------------------------------------
#' set.seed(18427)
#' X                   <- data.frame(Exposure = rbeta(100,3,1))
#' weights             <- runif(100)
#' normalized_weights  <- weights/sum(weights)
#' thetahat            <- 0.17
#' thetavar            <- 0.01
#' rr                  <- function(X, theta){theta*X^2 + 1}
#' paf.confidence(X, thetahat, rr, thetavar, weights = normalized_weights)
#' 
#' #Change the confidence level and paf method
#' paf.confidence(X, thetahat, rr,  thetavar, weights = normalized_weights, 
#'      method = "kernel", confidence = 90)
#' 
#' 
#' #Example 3: Multivariate Linear Relative Risk
#' #--------------------------------------------
#' set.seed(18427)
#' X1       <- rnorm(100,4,1)
#' X2       <- rnorm(100,2,0.4)
#' thetahat <- c(0.12, 0.03)
#' thetavar <- diag(c(0.01, 0.02))
#' 
#' #But the approximate method crashes due to operator
#' Xmean <- data.frame(Exposure = mean(X1), 
#'                     Covariate = mean(X2))
#' Xvar  <- var(cbind(X1, X2))
#' 
#' #When creating relative risks avoid using the $ operator
#' #as it doesn't work under approximate method of PAF
#' rr_not    <- function(X, theta){
#'                exp(theta[1]*X$Exposure + theta[2]*X$Covariate)
#'              }
#' rr_better <- function(X, theta){
#'                exp(theta[1]*X[,"Exposure"] + theta[2]*X[,"Covariate"])
#'              }
#'              
#' paf.confidence(Xmean, thetahat, rr_better, thetavar,
#'                method = "approximate", Xvar = Xvar)
#' }
#' \dontrun{
#' #Warning: $ operator in rr definitions don't work in approximate
#' paf.confidence(Xmean, thetahat, rr_not, thetavar,
#'                method = "approximate", Xvar = Xvar)
#' }
#' 
#' \dontrun{
#' #Example 4: Categorical Relative Risk & Exposure
#' #--------------------------------------------
#' set.seed(18427)
#' mysample  <- sample(c("Normal","Overweight","Obese"), 100, 
#'                    replace = TRUE, prob = c(0.4, 0.1, 0.5))
#' X        <- data.frame(Exposure = mysample)
#' 
#' thetahat <- c(1, 1.2, 1.5)
#' thetavar <- diag(c(0.1, 0.2, 0.3))
#' 
#' #Categorical relative risk function
#' rr <- function(X, theta){
#' 
#'    #Create return vector with default risk of 1
#'    r_risk <- rep(1, nrow(X))
#'    
#'    #Assign categorical relative risk
#'    r_risk[which(X[,"Exposure"] == "Normal")]      <- thetahat[1]
#'    r_risk[which(X[,"Exposure"] == "Overweight")]  <- thetahat[2]
#'    r_risk[which(X[,"Exposure"] == "Obese")]       <- thetahat[3]
#'    
#'    return(r_risk)
#' }
#' 
#' paf.confidence(X, thetahat, rr, thetavar, check_rr = FALSE)
#' 
#' 
#' #Example 5: Continuous Exposure and Categorical Relative Risk
#' #------------------------------------------------------------------
#' set.seed(18427)
#' 
#' #Assume we have BMI from a sample
#' BMI          <- data.frame(Exposure = rlnorm(100, 3.1, sdlog = 0.1))
#' 
#' #Theoretical minimum risk exposure is at 20kg/m^2 in borderline "Normal" category
#' BMI_adjusted <- BMI - 20
#' 
#' thetahat <- c(Malnourished = 2.2, Normal = 1, Overweight = 1.8, 
#'               Obese = 2.5)
#' thetavar <- diag(c(0.1, 0.2, 0.2, 0.1))
#' rr       <- function(X, theta){
#'      
#'      #Create return vector with default risk of 1
#'      r_risk <- rep(1, nrow(X))
#'    
#'      #Assign categorical relative risk
#'      r_risk[which(X[,"Exposure"] < 0)]             <- theta[1] #Malnourished
#'      r_risk[intersect(which(X[,"Exposure"] >= 0), 
#'                       which(X[,"Exposure"] < 5))]  <- theta[2] #Normal
#'      r_risk[intersect(which(X[,"Exposure"] >= 5), 
#'                       which(X[,"Exposure"] < 10))] <- theta[3] #Overweight
#'      r_risk[which(X[,"Exposure"] >= 10)]           <- theta[4] #Obese
#'    
#'    return(r_risk)
#' }
#' 
#' paf.confidence(BMI_adjusted, thetahat, rr, thetavar, check_exposure = FALSE)
#' 
#' #Example 6: Bivariate exposure and rr ("classical PAF")
#' #------------------------------------------------------------------
#' set.seed(18427)
#' mysample  <- sample(c("Exposed","Unexposed"), 1000, 
#'                 replace = TRUE, prob = c(0.1, 0.9))
#' X         <- data.frame(Exposure = mysample)
#' theta     <- c("Exposed" = 2.5, "Unexposed" = 1.2)  
#' thetavar  <- matrix(c(0.04, 0.02, 0.02, 0.03), ncol = 2)
#' rr        <- function(X, theta){
#'    
#'    #Create relative risk function
#'    r_risk <- rep(1, nrow(X))
#'    
#'    #Assign values of relative risk
#'    r_risk[which(X[,"Exposure"] == "Unexposed")] <- theta["Unexposed"]
#'    r_risk[which(X[,"Exposure"] == "Exposed")]   <- theta["Exposed"]
#'    
#'    return(r_risk)
#' }    
#' 
#' paf.confidence(X, theta, rr, thetavar)
#' 
#' #Example 7: Continuous exposure, several covariates
#' #------------------------------------------------------------------
#' X <- data.frame(Exposure = rbeta(100, 2, 3),
#'                 Age      = runif(100, 20, 100),
#'                 Sex      = sample(c("M","F"), 100, replace = TRUE),
#'                 BMI      = rlnorm(100, 3.2, 0.2))
#' thetahat <- c(-0.1, 0.05, 0.2, -0.4, 0.3, 0.1)
#' 
#' #Create variance of theta
#' almostvar <- matrix(runif(6^2), ncol = 6)
#' thetavar  <- t(almostvar) %*% almostvar
#' rr <- function(X, theta){
#'      #Create risk vector
#'      Risk    <- rep(1, nrow(X))
#'      
#'      #Identify subpopulations
#'      males   <- which(X[,"Sex"] == "M")
#'      females <- which(X[,"Sex"] == "F")
#'      
#'      #Calculate population specific rr
#'      Risk[males] <- theta[1]*X[males,"Exposure"] + 
#'                                       theta[2]*X[males,"Age"]^2 + 
#'                                       theta[3]*X[males,"BMI"]/2 
#'                                      
#'      Risk[females] <- theta[4]*X[females,"Exposure"] + 
#'                                       theta[5]*X[females,"Age"]^2 + 
#'                                       theta[6]*X[females,"BMI"]/2 
#'                                      
#'     return(Risk)
#' }
#' 
#' paf.confidence(X, thetahat, rr, thetavar)
#' }
#' @export

paf.confidence <- function(X, thetahat, rr,  thetavar = NA,
                           thetalow = NA, thetaup = NA,
                           method  = "empirical",
                           confidence_method = "bootstrap",
                           confidence = 95,
                           confidence_theta = 99,
                           nsim    =  1000, 
                           weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                           Xvar    = var(X), 
                           deriv.method.args = list(), 
                           deriv.method      = "Richardson",
                           adjust = 1, n = 512,
                           ktype  = "gaussian", 
                           bw     = "SJ",
                           check_exposure = TRUE, check_cft = TRUE, check_rr = TRUE,
                           check_xvar = TRUE, check_integrals = TRUE, check_thetas = TRUE,
                           force.min = FALSE){
  
  method            <- method[1]
  confidence_method <- confidence_method[1]
  

  switch (confidence_method,
            "inverse" ={
              if(any(is.na(thetavar))){stop("Please specify thetavar, variance of thetahat")}
              paf.confidence.inverse(X = X, thetahat = thetahat, rr = rr, thetavar = thetavar, 
                                     weights = weights, method = method,  
                                     nsim = nsim, confidence = confidence,
                                     deriv.method.args = deriv.method.args, deriv.method = deriv.method, 
                                     force.min = force.min,  check_thetas = check_thetas,
                                     Xvar = Xvar)
            },
            "one2one" ={
              if(any(is.na(thetalow)) || any(is.na(thetaup))){stop("Please specify thetalow and thetaup bounds of thetahat's CI")}
              paf.confidence.one2one(X = X, thetahat = thetahat, rr = rr, thetalow = thetalow, thetaup = thetaup, 
                                     weights =  weights, confidence = confidence, confidence_theta = confidence_theta,
                                     check_thetas = check_thetas, deriv.method.args = deriv.method.args,
                                     deriv.method = deriv.method, method = method, Xvar = Xvar,
                                     check_exposure = check_exposure, check_rr = check_rr, 
                                     check_integrals = check_integrals)
            },
            {
              if(any(is.na(thetavar))){stop("Please specify thetavar, variance of thetahat")}
              pif.confidence(X = X, thetahat = thetahat, rr = rr, thetavar = thetavar, 
                            cft=NA, method  = method, confidence_method = confidence_method, 
                            confidence = confidence, nsim    =  nsim, weights = weights, 
                            Xvar    = Xvar, deriv.method.args = deriv.method.args, 
                            deriv.method  = deriv.method, adjust = adjust, n = n,
                            ktype  = ktype,  bw = bw, check_exposure = check_exposure,
                            check_cft = check_cft, check_rr = check_rr,
                            check_xvar = check_xvar, check_integrals = check_integrals,
                            check_thetas = check_thetas, is_paf = TRUE)
            }
    )
}
