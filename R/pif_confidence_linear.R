#' @title  Confidence Intervals for the Potential Impact Fraction 
#' using Delta Method
#' 
#' @description Function that calculates approximate confidence intervals 
#' to the potential impact fraction.
#' 
#'@param X         Random sample (\code{data.frame}) which includes exposure and
#'  covariates.
#'  
#'@param thetahat  Estimator (\code{vector}) of \code{theta} for the Relative 
#'  Risk function.
#'  
#'@param rr        \code{function} for Relative Risk which uses parameter 
#'  \code{theta}. The order of the parameters shound be \code{rr(X, theta)}.
#' 
#' @param thetavar   Estimator of variance of \code{thetahat}
#' 
#' **Optional**
#' 
#'@param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'                 the Population Attributable Fraction \code{PAF} where 
#'                 counterfactual is 0 exposure
#' 
#' @param nsim      Number of simulations for estimation of variance
#' 
#' @param weights    Survey \code{weights} for the random sample \code{X}
#' 
#' @param confidence  Confidence level \% (default \code{95})
#' 
#' @param check_thetas Check that theta parameters are correctly inputed
#' 
#' @param check_integrals Check that counterfactual and relative risk's expected
#'   values are well defined for this scenario
#'      
#' @param check_exposure  Check that exposure \code{X} is positive and numeric
#'   
#' @param check_rr        Check that Relative Risk function \code{rr} equals 
#'   \code{1} when evaluated at \code{0}

#' @param is_paf    Boolean forcing evaluation of \code{\link{paf}}.
#' 
#' @author Rodrigo Zepeda Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho García Formentí \email{daliaf172@gmail.com}
#' 
#' @examples 
#' \dontrun{
#' #Example with risk given by HR (PAF)
#' set.seed(18427)
#' X <- rnorm(100,3,.5)
#' thetahat <- 0.12
#' thetavar <- 0.1
#' pif.confidence.linear(X, thetahat, function(X, theta){exp(theta*X)}, 
#'                       thetavar, nsim = 100)
#' 
#' #Example with linear counterfactual
#' cft      <- function(X){0.3*X}
#' pif.confidence.linear(X, thetahat, function(X, theta){exp(theta*X)}, 
#'                      thetavar, cft, nsim = 100)
#' 
#' #Example with theta and X multivariate
#' set.seed(18427)
#' X1 <- rnorm(100, 3,.5)
#' X2 <- rnorm(100,3,.5)
#' X  <- as.matrix(cbind(X1,X2))
#' thetahat <- c(0.1, 0.03)
#' thetavar <- matrix(c(0.1, 0, 0, 0.05), byrow = TRUE, nrow = 2)
#' rr       <- function(X, theta){
#'            .X <- as.matrix(X, ncol = 2)
#'            exp(theta[1]*.X[,1] + theta[2]*.X[,2])
#'            }
#' cft <- function(X){0.5*X}
#' pif.confidence.linear(X, thetahat, rr, thetavar, cft) 
#' }
#' 
#' @importFrom stats qnorm
#' @keywords internal
#' @export


pif.confidence.linear <- function(X, thetahat, rr, thetavar,
                                  cft = NA,
                                  weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                                  confidence = 95, nsim = 1000, check_thetas = TRUE,
                                  check_exposure = TRUE,
                                  check_rr = TRUE, check_integrals = TRUE,
                                  is_paf = FALSE){
  
  #Check confidence
  check.confidence(confidence)
  
  #Make thetavar matrix
  .thetavar <- as.matrix(thetavar)
  
  #Function for checking that thetas are correctly inputed
  if(check_thetas){ check.thetas(.thetavar, thetahat, NA, NA, "linear") }
  
  #Set the vector for the confidence intervals
  .ci <- c("Lower_CI" = NA, "Point_Estimate" = NA, 
           "Upper_CI" = NA, "Estimated_Variance" = NA)
  
  #Set Z for the confidence interval
  Z <- qnorm(1 - ((100-confidence)/200))
  
  #Get the point estimate and variance
  .ci["Point_Estimate"]     <- pif(X = X, thetahat = thetahat, rr = rr, cft = cft, 
                                   weights = weights, is_paf = is_paf,
                                   check_exposure = check_exposure, check_rr = check_rr,
                                   check_integrals = check_integrals)
  .ci["Estimated_Variance"] <- pif.variance.linear(X = X, thetahat = thetahat, rr = rr, 
                                                   thetavar = .thetavar, cft = cft, 
                                                   weights = weights, check_thetas = FALSE, 
                                                   nsim = nsim, is_paf = is_paf)
  .ci["Lower_CI"]           <- .ci["Point_Estimate"] - Z*sqrt(.ci["Estimated_Variance"])
  .ci["Upper_CI"]           <- .ci["Point_Estimate"] + Z*sqrt(.ci["Estimated_Variance"])
  
  # #Check ci < 1
  # if (.ci["Upper_CI"] >= 1){
  #   
  #   #Transform the problem to 0 <= 1 - pif to apply bounded CIs
  #   .transf_ci             <- 1 - .ci
  #   names(.transf_ci)      <- c("Upper_CI", "Point_Estimate", "Lower_CI", "Estimated_Variance")
  #   .transf_ci["Lower_CI"] <- (.transf_ci["Point_Estimate"]^2)/.transf_ci["Upper_CI"]           #Bound 1 - .ci below 
  #   .ci["Upper_CI"]        <- 1 - .transf_ci["Lower_CI"]                                        #Transform back
  #   
  # }
  
  .ci["Upper_CI"] <- min(.ci["Upper_CI"], 1)
  
  return(.ci)
  
}


