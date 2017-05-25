#' @title Population Attributable Fraction
#'   
#' @description Function for estimating the Population Attributable Fraction 
#'   \code{paf} from a cross-sectional sample of the exposure \code{X} with a 
#'   known Relative Risk function \code{rr} with meta-analytical parameter  \code{theta}, where 
#'   the Population Attributable Fraction is given by: \deqn{ PAF = 
#'   \frac{E_X\left[rr(X;\theta)\right]-1}{E_X\left[rr(X;\theta)\right]} .}{ PAF 
#'   = mean(rr(X; theta) - 1)/mean(rr(X; theta)).}
#'   
#' @param X         Random sample (\code{data.frame}) which includes exposure 
#'   and covariates or sample \code{mean} if \code{"approximate"} method is 
#'   selected.
#'   
#' @param thetahat  Asymptotically consistent or Fisher consistent
#'  estimator (\code{vector}) of \code{theta} for the Relative 
#'   Risk function \code{rr}.
#'   
#' @param rr        \code{function} for Relative Risk which uses parameter 
#'   \code{theta}. The order of the parameters should be \code{rr(X, theta)}.
#'   
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
#' @return paf      Estimate of Population Attributable Fraction.
#'   
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'   
#' @note \code{"approximate"} method should be the last choice. In practice 
#'   \code{"empirical"} should be preferred as convergence is faster in 
#'   simulations than \code{"kernel"}. In addition, the scope
#'   of \code{"kernel"} is limited as it does not work with multivariate 
#'   exposures \code{X}.
#'   
#' @note \code{\link{paf}} is a wrapper for \code{\link{pif}} with 
#'   counterfactual of theoretical minimum risk exposure (\code{rr} = 1).
#'   
#' @note For more information on kernels see \code{\link[stats]{density}}.
#'   
#' @note Do not use the \code{$} operator when using \code{"approximate"}
#'   \code{method}.
#'   
#' @details The Relative Risk function \code{rr} and counterfactual \code{cft} 
#'   should consider \code{X} to be a \code{data.frame}. Each row of 
#'   \code{X} is an "individual" and each column a variable (exposure or covariate) 
#'   of said individual.
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
#' #Using the empirical method
#' paf(X, thetahat, rr)
#' 
#' #Same example with kernel method
#' paf(X, thetahat, rr, method = "kernel")
#' 
#' #Same example with approximate method
#' Xmean <- data.frame(Exposure = mean(X[,"Exposure"]))
#' Xvar  <- var(X[,"Exposure"])
#' paf(Xmean, thetahat, rr, method = "approximate", Xvar = Xvar)
#' 
#' #Additional options for approximate:
#' paf(Xmean, thetahat, rr, method = "approximate", Xvar = Xvar, 
#'    deriv.method = "Richardson",  deriv.method.args = list(eps=1e-3, d=0.1))
#' 
#' #Example 2: Linear Relative Risk with weighted sample
#' #--------------------------------------------
#' set.seed(18427)
#' X                   <- data.frame(Exposure = rbeta(100,3,1))
#' weights             <- runif(100)
#' normalized_weights  <- weights/sum(weights)
#' thetahat            <- 0.12
#' rr                  <- function(X, theta){theta*X^2 + 1}
#' paf(X, thetahat, rr, weights = normalized_weights)
#' 
#'    
#' #Additional options for kernel:
#' paf(X, thetahat, rr, weights = normalized_weights, 
#'      method = "kernel", ktype = "cosine", bw = "nrd0")
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
#' #When creating relative risks avoid using the $ operator
#' #as it doesn't work under approximate method
#' rr_not    <- function(X, theta){
#'                exp(theta[1]*X$Exposure + theta[2]*X$Covariate)
#'              }
#' rr_better <- function(X, theta){
#'                exp(theta[1]*X[,"Exposure"] + theta[2]*X[,"Covariate"])
#'              }
#'              
#' #For the empirical method it makes no difference:              
#' paf(X, thetahat, rr_better) 
#' paf(X, thetahat, rr_not) 
#' 
#' 
#' #But the approximate method crashes due to operator
#' Xmean <- data.frame(Exposure = mean(X[,"Exposure"]), 
#'                     Covariate = mean(X[,"Covariate"]))
#' Xvar  <- var(X)
#' 
#' paf(Xmean, thetahat, rr_better, method = "approximate", Xvar = Xvar)
#' \dontrun{
#' #Warning: $ operator in rr definitions don't work in approximate
#' paf(Xmean, thetahat, rr_not, method = "approximate", Xvar = Xvar)
#' }
#' 
#' 
#' \dontrun{
#' #Warning: Multivariate cases cannot be evaluated with kernel method
#' paf(X, thetahat, rr, method = "kernel") 
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
#' paf(X, thetahat, rr, check_rr = FALSE)
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
#' paf(BMI_adjusted, thetahat, rr, check_exposure = FALSE)
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
#' paf(X, theta, rr)
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
#' 
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
#' paf(X, thetahat, rr)
#' 
#' 
#' @seealso  \code{\link{paf.confidence}} for confidence interval estimation, 
#'   \code{\link{pif}} for Potential Impact Fraction estimation.
#'   
#'   See \code{\link{paf.exponential}} and \code{\link{paf.linear}} for 
#'   fractions with ready-to-use exponential and linear Relative Risks 
#'   respectively.
#'   
#'   Sensitivity analysis plots can be done with \code{\link{paf.plot}}, and 
#'   \code{\link{paf.sensitivity}}.
#'   
#'   
#' @references Vander Hoorn, S., Ezzati, M., Rodgers, A., Lopez, A. D., & 
#'   Murray, C. J. (2004). \emph{Estimating attributable burden of disease from 
#'   exposure and hazard data. Comparative quantification of health risks: 
#'   global and regional burden of disease attributable to selected major risk 
#'   factors}. Geneva: World Health Organization, 2129-40.
#'   
#' @importFrom stats var
#'   
#' @export

paf <- function(X, thetahat,   rr,         
                method  = "empirical",
                weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                Xvar    = var(X), 
                deriv.method.args = list(), 
                deriv.method      = "Richardson",
                adjust = 1, n = 512,
                ktype  = "gaussian", 
                bw     = "SJ",
                check_exposure = TRUE, check_rr = TRUE, 
                check_integrals = TRUE){
  
 
  #Estimation of PAF
  .paf <- pif(X = X, thetahat = thetahat, rr = rr, cft=NA,
              method = method, weights = weights, 
              Xvar = Xvar, deriv.method.args = deriv.method.args,
              deriv.method = deriv.method, adjust = adjust, n = n,
              ktype = ktype, bw = bw, 
              check_exposure = check_exposure, check_integrals = check_integrals,
              check_rr = check_rr,  is_paf = TRUE)
  
  return(.paf)

}