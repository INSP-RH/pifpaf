#' @title Potential Impact Fraction Sensitivity Analysis plot
#'   
#' @description Function that plots a sensitivity analysis for the potential
#'   impact fraction by checking how estimates vary when reducing the exposure
#'   sample \code{X}.
#'   
#' @param X         Random sample (vector or matrix) which includes exposure and
#'   covariates. or sample mean if approximate method is selected.
#'   
#' @param thetahat  Estimator (vector or matrix) of \code{theta} for the 
#'   Relative Risk function.
#'   
#' @param rr        Function for Relative Risk which uses parameter 
#'   \code{theta}. The order of the parameters shound be \code{rr(X, theta)}.
#'   
#'   
#'   **Optional**
#'   
#' @param weights   Normalized survey \code{weights} for the sample \code{X}.
#'   
#' @param method    Either \code{empirical} (default), or  \code{kernel}.
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
#' @param check_integrals Check that counterfactual and relative risk's expected
#'   values are well defined for this scenario
#'   
#' @param check_exposure  Check that exposure \code{X} is positive and numeric
#'   
#' @param check_rr        Check that Relative Risk function \code{rr} equals 
#'   \code{1} when evaluated at \code{0}
#'   
#' @param nsim      Integer with number of samples to include (for each removal)
#'   in order to conduct sensitivity analysis
#'   
#' @param mremove   Limit to number of measurements of \code{X} to remove
#'   
#' @param title     String with plot title
#'   
#' @param legendtitle   String title for the legend of plot
#'   
#' @param xlab          String label for the X-axis of the plot (corresponding
#'   to "a")
#'   
#' @param ylab          String label for the Y-axis of the plot (corresponding
#'   to "b")
#'   
#' @param colors        String vector with colors for plots
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí \email{daliaf172@gmail.com}
#'   
#' @import  ggplot2
#'   
#' @seealso \code{\link{paf}} for Population Attributable Fraction estimation,
#'   \code{\link{pif.heatmap}} for sensitivity analysis of the counterfactual,
#'   \code{\link{pif.plot}} for a plot of potential impact fraction as a
#'   function of theta
#'   
#' @examples 
#' 
#' #Example 1
#' #------------------------------------------------------------------
#' set.seed(3284)
#' X  <- rnorm(250,3)                        #Sample
#' rr <- function(X,theta){exp(X*theta)}     #Relative risk
#' theta <- 0.1                              #Estimate of theta
#' \donttest{
#' paf.sensitivity(X, thetahat = theta, rr = rr)
#' }
#' 
#' #Save file
#' #ggsave("My Potential Impact Fraction Sensitivity Analysis.pdf")
#' 
#' #Example 2
#' #--------------------------------------------------------------
#' set.seed(3284)
#' X     <- rbeta(1000, 1, 0.2)
#' theta <- c(0.12, 1)
#' rr    <- function(X, theta){X*theta[1] + theta[2]}
#' 
#' \donttest{
#' #Using empirical method
#' paf.sensitivity(X, thetahat = theta, rr = rr, 
#'                 mremove = 100, nsim = 50, 
#'                 title = "My Sensitivity Analysis for example 1")
#' }          
#' \donttest{
#' #Same example with kernel
#' paf.sensitivity(X, theta, rr = rr, 
#'                  mremove = 100, nsim = 50, method = "kernel", 
#'                  title = "Sensitivity Analysis for example 1 using kernel")
#' }                 
#' 
#' #Example 4: Plot counterfactual with categorical risks
#' #------------------------------------------------------------------
#' set.seed(18427)
#' X        <- sample(c("Normal","Overweight","Obese"), 1000, 
#'                    replace = TRUE, prob = c(0.4, 0.1, 0.5))
#' thetahat <- c(1, 1.7, 2)
#' 
#' #Categorical relative risk function
#' rr <- function(X, theta){
#' 
#'    #Create return vector with default risk of 1
#'    r_risk <- rep(1, length(X))
#'    
#'    #Assign categorical relative risk
#'    r_risk[which(X == "Normal")]      <- thetahat[1]
#'    r_risk[which(X == "Overweight")]  <- thetahat[2]
#'    r_risk[which(X == "Obese")]       <- thetahat[3]
#'    
#'    return(r_risk)
#' }
#' 
#' 
#' 
#' \dontrun{
#' pafplot <- paf.sensitivity(X, thetahat = thetahat, rr = rr, 
#'                            title = "Sensitivity analysis of PAF for excess-weight",
#'                            colors = rainbow(4), 
#'                            legendtitle = "Values", 
#'                            check_exposure = FALSE, check_rr = FALSE)              
#' pafplot              
#' 
#' #You can edit pafplot as it is a ggplot object
#' pafplot + theme_classic()
#' }
#' 
#' @import ggplot2
#' @export


paf.sensitivity <- function(X, thetahat, rr,         
                            weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                            method  = c("empirical", "kernel"),
                            adjust = 1, n = 512,
                            ktype  = c("gaussian", "epanechnikov", "rectangular", "triangular", 
                                       "biweight","cosine", "optcosine"), 
                            bw     = c("SJ", "nrd0", "nrd", "ucv", "bcv"),
                            nsim = 50, mremove = min(nrow(as.matrix(X))/2,100), ylab  = "PAF", 
                            xlab  = "Number of randomly deleted observations for X", 
                            legendtitle = "Sensitivity Analysis",
                            title = "Population Attributable Fraction (PAF) Sensitivity Analysis",
                            colors = c("red", "deepskyblue", "gray75", "gray25"),
                            check_exposure = TRUE, check_rr = TRUE, check_integrals = TRUE){
  
  pif.sensitivity(X = X, thetahat = thetahat, rr = rr, weights = weights,
                  method = method, adjust = adjust, n = n, ktype = ktype,
                  nsim = nsim, mremove = mremove, ylab = ylab, xlab = xlab,
                  legendtitle = legendtitle, title = title, colors = colors,
                  check_exposure = check_exposure, check_rr = check_rr,
                  check_integrals = check_integrals, is_paf = TRUE)
  
}