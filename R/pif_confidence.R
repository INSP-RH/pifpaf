#' @title Confidence Intervals for the Potential Impact Fraction
#' 
#' @description Function that estimates confidence intervals for the Potential 
#'   Impact Fraction \code{\link{pif}} from a cross-sectional sample of 
#'   the exposure \code{X} with a known Relative Risk function \code{rr} with 
#'   parameter \code{theta}, where the Potential Impact Fraction is given
#'   by: \deqn{ PIF = 
#'   \frac{E_X\left[rr(X;\theta)\right] - 
#'   E_X\left[rr\big(\textrm{cft}(X);\theta\big)\right]} 
#'   {E_X\left[rr(X;\theta)\right]}. }{ PIF = (mean(rr(X; theta)) - 
#'   mean(rr(cft(X);theta)))/mean(rr(X; theta)) .}
#' 
#' @param X         Random sample (\code{data.frame}) which includes exposure 
#'   and covariates or sample \code{mean} if \code{"approximate"} method is 
#'   selected.
#'   
#' @param thetahat  Asymptotically consistent or Fisher consistent estimator
#'  (\code{vector}) of \code{theta} for the Relative 
#'   Risk function. \code{thetahat} should be asymptotically normal \code{(N(theta, var_of_theta))}
#'   with mean  \code{theta} and estimated variance \code{var_of_theta}.
#'   
#' @param rr        \code{function} for Relative Risk which uses parameter 
#'   \code{theta}. The order of the parameters shound be \code{rr(X, theta)}.
#'   
#' @param thetavar   Estimator of variance \code{var_of_theta}.
#' 
#' 
#' \strong{**Optional**}
#' 
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'   the Population Attributable Fraction \code{\link{paf}} where 
#'   counterfactual is that of the theoretical minimum risk exposure \code{X0} 
#'   such that \code{rr(X0,theta) = 1}.
#' 
#' @param weights   Normalized survey \code{weights} for the sample \code{X}.
#'   
#' @param nsim      Number of simulations for estimation of variance.
#'   
#' @param confidence Confidence level \% (default \code{95}). 
#'   
#' @param confidence_method  Either \code{bootstrap} (default), \code{linear}, 
#' \code{loglinear}. See \code{\link{paf}} details for additional information.
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
#' @param check_xvar \code{boolean} Check \code{Xvar} is covariance matrix.
#'   
#' @param check_integrals \code{boolean}  Check that counterfactual \code{cft} 
#'   and relative risk's \code{rr} expected values are well defined for this 
#'   scenario.
#'   
#' @param check_rr        \code{boolean} Check that Relative Risk function \code{rr} equals 
#'   \code{1} when evaluated at \code{0}.
#'   
#' @param is_paf    Boolean forcing evaluation of \code{\link{paf}}. This forces
#'   the \code{\link{pif}} function to ignore the inputed counterfactual and set 
#'   it to the theoretical minimum risk value of \code{rr = 1}.
#' 
#' @return pifvec Vector with lower (\code{"Lower_CI"}), and upper 
#'   (\code{"Upper_CI"}) confidence bounds for the \code{\link{pif}} as well as
#'   point estimate \code{"Point_Estimate"} and estimated variance 
#'   of \code{log(pif)} (if \code{confidence_method} is \code{"loglinear"}).
#'   
#' @note For more information on kernels see \code{\link[stats]{density}}.
#'   
#' @note Do not use the \code{$} operator when using \code{"approximate"}
#'   \code{method}.
#'   
#'   
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#' 
#' @seealso  
#' 
#' See \code{\link{paf.confidence}} for confidence interval estimation of
#'   \code{\link{paf}}, and  \code{\link{pif}} for only point estimate.
#'   
#'   Sensitivity analysis plots can be done with \code{\link{paf.plot}}, and 
#'   \code{\link{paf.sensitivity}}
#'   
#' @examples 
#' #Example 1: Exponential Relative Risk
#' #--------------------------------------------
#' set.seed(18427)
#' X        <- data.frame(Exposure = rnorm(100,3,1))
#' thetahat <- 0.12
#' thetavar <- 0.02
#' rr       <- function(X, theta){exp(theta*X)}
#' 
#' 
#' #Counterfactual of halving exposure
#' cft   <- function(X){ 0.5*X }
#' 
#' #Using bootstrap method
#' pif.confidence(X, thetahat, rr, thetavar, cft)
#' 
#' \dontrun{
#' #Same example with loglinear method
#' pif.confidence(X, thetahat, rr, thetavar, cft, confidence_method = "loglinear")
#' 
#' #Same example with linear method (usually the widest and least precise)
#' pif.confidence(X, thetahat, rr, thetavar, cft, confidence_method = "linear")
#' 
#' 
#' #Example 2: Linear Relative Risk
#' #--------------------------------------------
#' set.seed(18427)
#' X        <- data.frame(Exposure = rbeta(100,3,1))
#' thetahat <- 0.17
#' thetavar <- 0.01
#' rr       <- function(X, theta){theta*X + 1}
#' cft      <- function(X){ 0.5*X }
#' weights             <- runif(100)
#' normalized_weights  <- weights/sum(weights)
#' pif.confidence(X, thetahat, rr, thetavar, cft, weights = normalized_weights)
#' 
#' #Same example with more complex counterfactual that reduces 
#' #only the values > 0.75 are halved
#' cft       <- function(X){
#' 
#'    #Indentify the ones with "a lot" of exposure:
#'    where_excess_exposure    <- which(X[,"Exposure"] > 0.75)             
#'    
#'    #Halve their exposure
#'    X[where_excess_exposure, "Exposure"] <- 
#'             X[where_excess_exposure, "Exposure"]/2  
#'    return(X)
#' }
#' pif.confidence(X, thetahat, rr, thetavar, cft, weights = normalized_weights)
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
#' #Creating a counterfactual. 
#' cft  <- function(X){
#'    Y               <- X
#'    Y[,"Exposure"]  <- 0.5*X[,"Exposure"]
#'    Y[,"Covariate"] <- 1.1*X[,"Covariate"] + 1
#'    return(Y)
#' }
#' 
#' pif.confidence(Xmean, thetahat, rr_better, thetavar, cft, 
#' method = "approximate", Xvar = Xvar) 
#' }
#' 
#' \dontrun{
#' #Warning: $ operator in rr definitions don't work in approximate
#' pif.confidence(Xmean, thetahat, rr_not, thetavar, cft, method = "approximate", Xvar = Xvar)
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
#' 
#' #Counterfactual of reducing all obesity to normality
#' cft <- function(X){
#'    X[which(X[,"Exposure"] == "Obese"),] <- "Normal"
#'    return(X)
#' }
#' 
#' pif.confidence(X, thetahat, rr, thetavar, cft, check_rr = FALSE)
#' 
#' 
#' #Example 5: Categorical Relative Risk & continuous exposure
#' #----------------------------------------------------------
#' set.seed(18427)
#' BMI      <- data.frame(Exposure = rlnorm(100, 3.1, sdlog = 0.1))
#' 
#' #Theoretical minimum risk exposure is 20kg/m^2 in borderline "Normal" category
#' BMI_adjusted <- BMI - 20
#' 
#' thetahat <- c(Malnourished = 2.2, Normal = 1, Overweight = 1.8, 
#'               Obese = 2.5)
#' thetavar <- diag(c(0.1, 0.2, 0.2, 0.1))
#' 
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
#' #Counterfactual of everyone in normal range
#' cft <- function(bmi){
#'      bmi           <- data.frame(rep(2.5, nrow(bmi)), ncol = 1)
#'      colnames(bmi) <- c("Exposure")
#'      return(bmi)
#' }
#' 
#' pif.confidence(BMI_adjusted, thetahat, rr, thetavar, cft, 
#'                 check_exposure = FALSE, method = "empirical")
#' 
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
#' #Counterfactual of reducing the exposure in half of the individuals
#' cft <- function(X){
#' 
#'    #Find out which ones are exposed
#'    Xexp  <- which(X[,"Exposure"] == "Exposed")
#'    
#'    #Use only half of the exposed randomly
#'    reduc <- sample(Xexp, length(Xexp)/2)
#'    
#'    #Unexpose those individuals
#'    X[reduc, "Exposure"] <- "Unexposed"
#'    
#'    return(X)
#' }
#' 
#' pif.confidence(X, theta, rr, thetavar, cft)
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
#' 
#' thetavar <- t(almostvar) %*% almostvar
#'
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
#' #Counterfactual of reducing BMI
#' cft <- function(X){
#'     excess_bmi           <- which(X[,"BMI"] > 25)
#'     X[excess_bmi,"BMI"]  <- 25
#'     return(X)
#' }
#' 
#' pif.confidence(X, thetahat, rr, thetavar, cft)
#' }
#' @export

pif.confidence <- function(X, thetahat, rr, thetavar, 
                   cft = NA,
                   method  = "empirical",
                   confidence_method = "bootstrap",
                   confidence = 95,
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
                   is_paf = FALSE){
  
  
  #Choose method
  method            <- method[1]
  confidence_method <- confidence_method[1]
  
  switch (method,
          "kernel"      = {
            pif.confidence.bootstrap(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr, cft = cft, weights = weights,
                                     method = "kernel", nboost = nsim, n=n, adjust = adjust, confidence = confidence,
                                     ktype = ktype, bw = bw, check_exposure = check_exposure, check_rr = check_rr,
                                     check_integrals = check_integrals, check_thetas = check_thetas, is_paf = is_paf)
          },
          "approximate" = {
            switch(confidence_method,
              "linear"    = {
                pif.confidence.approximate(Xmean = X, Xvar = Xvar, thetahat = thetahat, thetavar = thetavar, rr = rr,
                                           cft = cft, check_thetas = check_thetas, check_cft = check_cft,
                                           check_xvar = check_xvar, check_rr = check_rr, check_integrals = check_integrals,
                                           check_exposure = check_exposure, deriv.method.args = deriv.method.args,
                                           deriv.method = deriv.method, nsim = nsim, confidence = confidence,
                                           is_paf = is_paf)
              },
              {
                pif.confidence.approximate.loglinear(Xmean = X, Xvar = Xvar, thetahat = thetahat, thetavar = thetavar, rr = rr,
                                                     cft = cft, deriv.method.args = deriv.method.args, deriv.method = deriv.method,
                                                     check_exposure = check_exposure, check_rr = check_rr, check_integrals = check_integrals,
                                                     nsim = nsim, confidence = confidence, check_thetas = check_thetas,
                                                     is_paf = is_paf)
              }
            )
            
          },
          {
            switch (confidence_method,
                    "linear"    = {
                      pif.confidence.linear(X = X, thetahat = thetahat, rr = rr, thetavar = thetavar,  cft = cft, weights = weights, 
                                            confidence = confidence, nsim=nsim, check_thetas = check_thetas,
                                            check_exposure = check_exposure, check_rr = check_rr,
                                            check_integrals = check_integrals,  is_paf = is_paf)
                    },
                    "loglinear" = {
                      pif.confidence.loglinear(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr, cft = cft, weights = weights,
                                               nsim = nsim, confidence=confidence, check_thetas = check_thetas, 
                                               check_exposure = check_exposure, 
                                               check_cft = check_cft, is_paf = is_paf)
                    },
                    {
                      pif.confidence.bootstrap(X = X, thetahat = thetahat, thetavar = thetavar, rr = rr, cft = cft, weights = weights,
                                               method = "empirical", nboost = nsim, confidence = confidence,
                                               check_exposure = check_exposure, check_rr = check_rr,
                                               check_integrals = check_integrals, check_thetas = check_thetas,
                                               is_paf = is_paf)
                    }
            )
          }
  )
}