#' @title Potential Impact Fraction
#'   
#' @description Function for estimating the Potential Impact Fraction \code{pif}
#'   from a cross-sectional sample of the exposure \code{X} with known Relative 
#'   Risk function \code{rr} with parameter \code{theta}, where the Potential 
#'   Impact Fraction is given by: \deqn{ PIF = 
#'   \frac{E_X\left[rr(X;\theta)\right] - 
#'   E_X\left[rr\big(\textrm{cft}(X);\theta\big)\right]} 
#'   {E_X\left[rr(X;\theta)\right]}. }{ PIF = (mean(rr(X; theta)) - 
#'   mean(rr(cft(X);theta)))/mean(rr(X; theta)) .}
#'   
#' @param X         Random sample (\code{data.frame}) which includes exposure 
#'   and covariates or sample \code{mean} if \code{"approximate"} method is 
#'   selected.
#'   
#' @param thetahat  Asymptotically consistent or Fisher consistent
#'  estimator (\code{vector}) of \code{theta} for the Relative 
#'   Risk function.
#'   
#' @param rr        \code{function} for Relative Risk which uses parameter 
#'   \code{theta}. The order of the parameters should be \code{rr(X, theta)}.
#'   
#'   
#'  \strong{**Optional**}
#'   
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'   the Population Attributable Fraction \code{\link{paf}} where 
#'   counterfactual is that of a theoretical minimum risk exposure \code{X0} 
#'   such that \code{rr(X0,theta) = 1}.
#'   
#' @param weights   Normalized survey \code{weights} for the sample \code{X}.
#'   
#' @param method    Either \code{"empirical"} (default), \code{"kernel"} or 
#'   \code{"approximate"}.
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
#' @param check_integrals Check that counterfactual of theoretical minimum risk 
#'   exposure and relative risk's expected values are well defined for this 
#'   scenario.
#'   
#' @param check_exposure  Check exposure \code{X} is positive and numeric.
#'   
#' @param check_rr        Check that Relative Risk function \code{rr} equals 
#'   \code{1} when evaluated at \code{0}.
#'   
#' @param is_paf    Boolean forcing evaluation of \code{\link{paf}}. This forces
#'   the \code{pif} function ignore the inputed counterfactual and set the 
#'   relative risk to the theoretical minimum risk value of \code{1}.
#'   
#' @return pif      Estimate of Potential Impact Fraction.
#'   
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'   
#' @note For more information on kernels see \code{\link[stats]{density}}.
#'   
#' @note Do not use the \code{$} operator when using \code{"approximate"}
#'   \code{method}.
#'   
#' @details The \code{"empirical"} method estimates the \code{pif} by \deqn{ PIF
#'   = 1 - \frac{\sum\limits_{i=1}^{n}w_i rr\big(cft(X_i); 
#'   \theta\big)}{\sum\limits_{i=1}^{n} w_i rr(X_i; \theta)}. }{ PIF = 1 - 
#'   weighted.mean(rr(cft(X), theta), weights)/ weighted.mean(rr(cft(X), theta),
#'   weights). }
#'   
#'   The \code{"kernel"} method approximates the \code{\link[stats]{density}} of
#'   the exposure \code{X} and estimates its expected value from that 
#'   approximation: 
#'   \deqn{ PIF = 1 - 
#'         \frac{\int\limits_{-\infty}^{\infty} rr\big(cft(X);\theta \big) \hat{f}(x) dx
#'   }{
#'         \int\limits_{-\infty}^{\infty} rr\big(cft(X);\theta \big) \hat{f}(x) dx}.
#'   }{ 
#'   PIF = 1 - integrate(rr(cft(X), theta)*f(x))/integrate(rr(X, theta)*f(x)).
#'   }
#'   
#'   The \code{"approximate"} method conducts a Laplace approximation of the \code{pif}.
#'   Additional information on the methods is dicussed in the package's vignette:
#'   \code{browseVignettes("pifpaf")}.
#'   
#'   In practice \code{"approximate"} method should be the last choice. 
#'   Simulations have shown that \code{"empirical"}'s convergence is faster than
#'   \code{"kernel"} for most functions. In addition, the scope of 
#'   \code{"kernel"} is limited as it does not work with multivariate exposure 
#'   data \code{X}.
#'   
#' @examples 
#' 
#' #Example 1: Exponential Relative Risk
#' #--------------------------------------------
#' set.seed(18427)
#' X        <- data.frame(Exposure = rnorm(100,3,1))
#' thetahat <- 0.12
#' rr       <- function(X, theta){exp(theta*X)}
#' 
#' #Without specifying counterfactual pif matches paf
#' pif(X, thetahat, rr)
#' paf(X, thetahat, rr)
#' 
#' #Same example with kernel method
#' pif(X, thetahat, rr, method = "kernel")
#' 
#' #Same example with approximate method
#' Xmean <- data.frame(Exposure = mean(X[,"Exposure"]))
#' Xvar  <- var(X[,"Exposure"])
#' pif(Xmean, thetahat, rr, method = "approximate", Xvar = Xvar)
#' 
#' #Same example considering counterfactual of halving exposure
#' cft   <- function(X){ 0.5*X }
#' pif(X, thetahat, rr, cft, method = "empirical")
#' 
#' #Example 2: Linear Relative Risk
#' #--------------------------------------------
#' set.seed(18427)
#' X        <- data.frame(Exposure = rbeta(100,3,1))
#' thetahat <- 0.12
#' rr       <- function(X, theta){theta*X + 1}
#' cft      <-  function(X){ 0.5*X }
#' weights             <- runif(100)
#' normalized_weights  <- weights/sum(weights)
#' pif(X, thetahat, rr, cft, weights = normalized_weights)
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
#' pif(X, thetahat, rr, cft, weights = normalized_weights)
#' 
#' 
#' #Example 3: Multivariate Linear Relative Risk
#' #--------------------------------------------
#' set.seed(18427)
#' X1       <- rnorm(100,4,1)
#' X2       <- rnorm(100,2,0.4)
#' X        <- data.frame(Exposure = X1, Covariate = X2)
#' thetahat <- c(0.12, 0.03)
#' 
#' #When creating relative risks and counterfactuals avoid using $ operator
#' #as it doesn't work under approximate method
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
#' pif(X, thetahat, rr_better, cft) 
#' 
#' #Same multivariate example for approximate method calculating 
#' #mean and variance
#' Xmean <- data.frame(Exposure = mean(X$Exposure), 
#'                    Covariate = mean(X$Covariate))
#' Xvar  <- var(X)
#' pif(Xmean, thetahat, rr_better, method = "approximate", Xvar = Xvar)
#' 
#' \dontrun{
#' #The one with $ operators doesn't work:
#' pif(Xmean, thetahat, rr_not, method = "approximate", Xvar = Xvar)
#' }
#' \dontrun{
#' #Warning: Multivariate cases cannot be evaluated with kernel method
#' pif(X, thetahat, rr_better, method = "kernel") 
#' }
#' 
#' #Example 4: Categorical Relative Risk & Exposure
#' #--------------------------------------------
#' set.seed(18427)
#' mysample  <- sample(c("Normal","Overweight","Obese"), 100, 
#'                    replace = TRUE, prob = c(0.4, 0.1, 0.5))
#' X        <- data.frame(Exposure = mysample)
#' 
#' thetahat <- c(1, 1.2, 1.5)
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
#' pif(X, thetahat, rr, check_rr = FALSE)
#' 
#' #Counterfactual of reducing all obesity to normality
#' cft <- function(X){
#'    X[which(X[,"Exposure"] == "Obese"),] <- "Normal"
#'    return(X)
#' }
#' 
#' pif(X, thetahat, rr, cft, check_rr = FALSE)
#' 
#' #Example 5: Categorical Relative Risk & continuous exposure
#' #----------------------------------------------------------
#' set.seed(18427)
#' BMI      <- data.frame(Exposure = rlnorm(100, 3.1, sdlog = 0.1))
#' 
#' #Theoretical minimum risk exposure is at 20kg/m^2 in borderline "Normal" category
#' BMI_adjusted <- BMI - 20
#' 
#' thetahat <- c(Malnourished = 2.2, Normal = 1, Overweight = 1.8, 
#'               Obese = 2.5)
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
#' pif(BMI_adjusted, thetahat, rr, cft, 
#'     check_exposure = FALSE, method = "empirical")
#' 
#' 
#' #Example 6: Bivariate exposure and rr ("classical PAF")
#' #------------------------------------------------------------------
#' set.seed(18427)
#' mysample  <- sample(c("Exposed","Unexposed"), 1000, 
#'                 replace = TRUE, prob = c(0.1, 0.9))
#' X         <- data.frame(Exposure = mysample)
#' theta     <- c("Exposed" = 2.5, "Unexposed" = 1.2)   
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
#' pif(X, theta, rr, cft)
#' 
#' #Example 7: Continuous exposure, several covariates
#' #------------------------------------------------------------------
#' X <- data.frame(Exposure = rbeta(100, 2, 3),
#'                 Age      = runif(100, 20, 100),
#'                 Sex      = sample(c("M","F"), 100, replace = TRUE),
#'                 BMI      = rlnorm(100, 3.2, 0.2))
#' thetahat <- c(-0.1, 0.05, 0.2, -0.4, 0.3, 0.1)
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
#' pif(X, thetahat, rr, cft)
#' 
#' @seealso  See \code{\link{pif.confidence}} for confidence interval estimation, 
#'   and \code{\link{paf}} for Population Attributable Fraction estimation.
#'   
#'   Sensitivity analysis plots can be done with \code{\link{pif.plot}}, 
#'   \code{\link{pif.sensitivity}}, and \code{\link{pif.heatmap}}.
#'   
#' @references Vander Hoorn, S., Ezzati, M., Rodgers, A., Lopez, A. D., & 
#'   Murray, C. J. (2004). \emph{Estimating attributable burden of disease from 
#'   exposure and hazard data. Comparative quantification of health risks: 
#'   global and regional burden of disease attributable to selected major risk 
#'   factors}. Geneva: World Health Organization, 2129-40.
#'   
#' @export


pif <- function(X, thetahat, rr,         
                cft     = NA,
                method  = "empirical",
                weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                Xvar    = var(X), 
                deriv.method.args = list(), 
                deriv.method      = "Richardson",
                adjust  = 1, n = 512,
                ktype   = "gaussian", 
                bw      = "SJ",
                check_exposure = TRUE, 
                check_integrals = TRUE,
                check_rr = TRUE, 
                is_paf = FALSE){
  
  #Get method from vector
  .method <- as.vector(method)[1]
  
  #Check if counterfactual is na then estimate paf
  if (!is.function(cft)){ is_paf <- TRUE}
  
  #Check X is data.frame
  if (!is.data.frame(X)){
    warning("Exposure X should be a data.frame.")
  }
  
  switch(.method,
         empirical   = {
           .pif <- pif.empirical(X = X, thetahat = thetahat, rr = rr, cft = cft, 
                                 weights = weights, check_exposure = check_exposure, 
                                 check_rr = check_rr, check_integrals = check_integrals,
                                 is_paf = is_paf)
           
         }, 
         kernel      = {
           .pif <- pif.kernel(X = X, thetahat = thetahat, rr = rr, cft = cft, 
                              weights = weights, adjust = adjust, n = n,
                              ktype = ktype, bw = bw, 
                              check_exposure = check_exposure,
                              check_rr = check_rr, check_integrals = check_integrals,
                              is_paf = is_paf)
         }, 
         approximate = {
           .pif <- pif.approximate(X = X, Xvar = Xvar, thetahat = thetahat, rr = rr, 
                                   cft = cft, deriv.method.args = deriv.method.args,
                                   deriv.method = deriv.method, check_exposure = check_exposure,
                                   check_rr = check_rr, check_integrals = check_integrals,
                                   is_paf = is_paf)
         },{
           stop("Please specify method as either empirical, kernel or approximate")
         }
         
  )
  return(.pif)
}


