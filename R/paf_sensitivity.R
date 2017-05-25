#' @title Population Attributable Fraction Sensitivity Analysis Plot
#'   
#' @description Function that plots a sensitivity analysis for the Population 
#'   Attributable Fraction \code{\link{paf}} by checking how estimates vary when
#'   reducing the exposure's sample \code{X}.
#'   
#' @param X         Random sample (\code{data.frame}) which includes exposure 
#'   and covariates or sample \code{mean} if \code{"approximate"} method is 
#'   selected.
#'   
#' @param thetahat  Asymptotically consistent or Fisher consistent estimator
#'  (\code{vector}) of \code{theta} for the
#'   Relative Risk function.
#'   
#' @param rr        \code{function} for Relative Risk which uses parameter 
#'   \code{theta}. The order of the parameters should be \code{rr(X, theta)}.
#'   
#'  \strong{**Optional**}
#'   
#' @param weights   Normalized survey \code{weights} for the sample \code{X}.
#'   
#' @param method    Either \code{"empirical"} (default), \code{"kernel"} or 
#'   \code{"approximate"}. For details on estimation methods see 
#'   \code{\link{pif}}.
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
#' @param nsim      Integer with number of samples to include (for each removal)
#'   in order to conduct sensitivity analysis. See details for additional information.
#'   
#' @param mremove   Limit number of measurements of \code{X} to remove when
#'   resampling. See details for additional information.
#'   
#' @param title \code{string} Title of plot.
#'   
#' @param legendtitle   String title for the legend of plot.
#'   
#' @param xlab          \code{string} label for the X-axis of the plot.
#'   
#' @param ylab          \code{string} label for the Y-axis of the plot.
#'   
#' @param colors        String vector with colors for the plot.
#'   
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'   
#' @details \code{paf.sensitivity} conducts a sensitivity analysis of the 
#'   \code{\link{paf}} estimate by removing \code{mremove} elements \code{nsim}
#'   times and re-estimating \code{\link{paf}} with the reduced sample.
#'   
#' @return plotpaf      \code{\link[ggplot2]{ggplot}} object plotting a 
#'   sensitivity analysis of \code{\link{paf}}.
#'   
#' @seealso \code{\link{paf}} for Population Attributable Fraction estimation, 
#'   \code{\link{paf.plot}} for a plot of Population Attributable Fraction as a 
#'   function of the relative risk's parameter \code{theta}.
#'   
#' @examples 
#' \dontrun{
#' #Example 1
#' #------------------------------------------------------------------
#' set.seed(3284)
#' X  <- data.frame(rnorm(250,3))            #Sample
#' rr <- function(X,theta){exp(X*theta)}     #Relative risk
#' theta <- 0.1                              #Estimate of theta
#' paf.sensitivity(X, thetahat = theta, rr = rr)
#' 
#' 
#' #Save file
#' #require(ggplot2)
#' #ggsave("My Population Attributable Fraction Sensitivity Analysis.pdf")
#' 
#' #Example 2
#' #--------------------------------------------------------------
#' set.seed(3284)
#' X     <- data.frame(rbeta(1000, 1, 0.2))
#' theta <- c(0.12, 1)
#' rr    <- function(X, theta){X*theta[1] + theta[2]}
#' 
#' 
#' #Using empirical method
#' paf.sensitivity(X, thetahat = theta, rr = rr, 
#'                 mremove = 100, nsim = 50, 
#'                 title = "My Sensitivity Analysis for example 1")
#'                 
#' #Same example with kernel
#' paf.sensitivity(X, theta, rr = rr, 
#'                  mremove = 100, nsim = 50, method = "kernel", 
#'                  title = "Sensitivity Analysis for example 1 using kernel")
#'                  
#' 
#' #Example 3: Plot counterfactual with categorical risks
#' #------------------------------------------------------------------
#' set.seed(18427)
#' X        <- data.frame(Exposure = 
#'                sample(c("Normal","Overweight","Obese"), 1000, 
#'                       replace = TRUE, prob = c(0.4, 0.1, 0.5)))
#' thetahat <- c(1, 1.7, 2)
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
#' pafplot <- paf.sensitivity(X, thetahat = thetahat, rr = rr, 
#'                            title = "Sensitivity analysis of PAF for excess-weight",
#'                            colors = rainbow(4), 
#'                            legendtitle = "Values", 
#'                            check_exposure = FALSE, check_rr = FALSE)              
#' pafplot              
#' 
#' #You can edit pafplot as it is a ggplot object
#' #require(ggplot2)
#' #pafplot + theme_classic()
#' }
#' 
#' @import ggplot2
#' @export


paf.sensitivity <- function(X, thetahat, rr,   
                            method  = "empirical",
                            weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                            nsim    = 50, 
                            mremove = min(nrow(as.matrix(X))/2,100), 
                            adjust  = 1, 
                            n       = 512,
                            ktype   = "gaussian",
                            bw      = "SJ", 
                            ylab    = "PAF", 
                            xlab    = "Number of randomly deleted observations for X", 
                            legendtitle = "Sensitivity Analysis",
                            title   = paste0("Population Attributable Fraction (PAF) ", 
                                             "Sensitivity Analysis"),
                            colors  = c("red", "deepskyblue", "gray75", "gray25"),
                            check_exposure = TRUE, check_rr = TRUE, check_integrals = TRUE){
  
  pif.sensitivity(X = X, thetahat = thetahat, rr = rr, cft=NA, method = method,
                  weights = weights, nsim = nsim, mremove = mremove,
                   adjust = adjust, n = n, ktype = ktype, bw=bw,
                  ylab = ylab, xlab = xlab,
                  legendtitle = legendtitle, title = title, colors = colors,
                  check_exposure = check_exposure, check_rr = check_rr,
                  check_integrals = check_integrals, is_paf = TRUE)
  
}