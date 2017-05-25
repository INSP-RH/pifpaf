#' @title Potential Impact Fraction Sensitivity Analysis plot
#'   
#' @description Function that plots a sensitivity analysis for the Potential 
#'   Impact Fraction \code{\link{pif}} by checking how estimates vary when reducing 
#'   the exposure's sample \code{X}.
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
#'   
#'   \strong{**Optional**}
#'   
#' @param cft       \code{function} \code{cft(X)} for counterfactual. Leave empty for 
#'   the Population Attributable Fraction \code{\link{paf}} where 
#'   counterfactual is that of a theoretical minimum risk exposure 
#'   \code{X0} such that \code{rr(X0,theta) = 1}.
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
#' @param title \code{string} title of plot.
#'   
#' @param legendtitle   \code{string} title for the legend of plot.
#'   
#' @param xlab          \code{string} label for the X-axis of the plot.
#'   
#' @param ylab          \code{string} label for the Y-axis of the plot.
#'   
#' @param colors        \code{string} vector with colours for plots.
#'   
#' @param is_paf    Boolean forcing evaluation of \code{\link{paf}}. This forces
#'   the \code{\link{pif}} function ignore the inputed counterfactual and set it to the
#'   theoretical minimum risk value of \code{1}.
#'   
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'   
#' @return plotpif      \code{\link[ggplot2]{ggplot}} object plotting a
#'   sensitivity analysis of \code{\link{pif}}.
#'   
#' @seealso 
#' 
#' See \code{\link{pif}} for Potential Impact Fraction estimation, 
#'   \code{\link{pif.heatmap}} for sensitivity analysis of counterfactual, 
#'   \code{\link{pif.plot}} for a plot of Potential Impact Fraction as a
#'   function of the relative risk's parameter \code{theta}.
#'   
#' @details \code{pif.sensitivity} conducts a sensitivity analysis of the 
#'   \code{\link{pif}} estimate by removing \code{mremove} elements \code{nsim}
#'   times and re-estimating \code{\link{pif}} with the reduced sample.
#'    
#' @examples 
#' \dontrun{
#' #Example 1
#' #------------------------------------------------------------------
#' set.seed(3284)
#' X  <- data.frame(Exposure = rnorm(250,3)) #Sample
#' rr <- function(X,theta){exp(X*theta)}     #Relative risk
#' theta <- 0.1                              #Estimate of theta
#' 
#' pif.sensitivity(X, thetahat = theta, rr = rr)
#' 
#' 
#' #Save file
#' #require(ggplot2)
#' #ggsave("My Potential Impact Fraction Sensitivity Analysis.pdf")
#' 
#' #Example 2
#' #--------------------------------------------------------------
#' set.seed(3284)
#' X     <- data.frame(Exposure = rbeta(1000, 1, 0.2))
#' theta <- c(0.12, 1)
#' rr    <- function(X, theta){X*theta[1] + theta[2]}
#' cft   <- function(X){X/2}
#' 
#' 
#' #Using empirical method
#' pif.sensitivity(X, thetahat = theta, rr = rr, cft = cft,
#'                 mremove = 100, nsim = 50, 
#'                 title = "My Sensitivity Analysis for example 1")
#'                 
#' #Same example with kernel
#' pif.sensitivity(X, theta, rr = rr, cft = cft,
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
#' #Counterfactual of halving the percent of obesity and overweight cases
#' #to normality
#' cft <- function(X){
#' 
#'    #Find the overweight and obese individuals
#'    which_obese <- which(X == "Obese")
#'    which_over  <- which(X == "Overweight")
#'    
#'    #Reduce per.over % of overweight and per.obese % of obese
#'    X[sample(which_obese, floor(length(which_obese)*0.5)),1] <- "Normal"
#'    X[sample(which_over,  floor(length(which_over)*0.5)),1]  <- "Normal"
#'    
#'    return(X)
#' }
#' 
#' 
#' pifplot <- pif.sensitivity(X, thetahat = thetahat, rr = rr, cft = cft, 
#'                            title = "Sensitivity analysis of PIF for excess-weight",
#'                            colors = rainbow(4), 
#'                            legendtitle = "Values", 
#'                            check_exposure = FALSE, check_rr = FALSE)              
#' pifplot              
#' 
#' #You can edit pifplot as it is a ggplot object
#' #require(ggplot2)
#' #pifplot + theme_classic()
#' }
#' 
#' @import ggplot2
#' @export


pif.sensitivity <- function(X, thetahat, rr,         
                            cft = NA,
                            method  = "empirical",
                            weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                            nsim    = 50, 
                            mremove = min(nrow(as.matrix(X))/2,100), 
                            adjust = 1, 
                            n = 512,
                            ktype  = "gaussian", 
                            bw     = "SJ", 
                            ylab  = "PIF", 
                            xlab  = "Number of randomly deleted observations for X", 
                            legendtitle = "Sensitivity Analysis",
                            title = "Potential Impact Fraction (PIF) Sensitivity Analysis",
                            colors = c("red", "deepskyblue", "gray75", "gray25"),
                            check_exposure = TRUE, check_rr = TRUE, check_integrals = TRUE,
                            is_paf = FALSE){
  
  #Set X as matrix
  .X       <- as.data.frame(X)
  
  #Check m and n are correctly defined
  if(mremove <= 0){
    stop("m (maximum observations to remove) must be positive")
  } else if(mremove >= nrow(.X)){
    stop("m (maximum observations to remove) must be less than the amount of observations")
  }
  
  if(nsim <= 0){
    stop("n (number of samples) must be positive")
  } else if(nsim >= 500){
    warning("n (number of samples) is too big")
  }
  
  #Limit m values
  .m       <- min(ceiling(mremove), nrow(X))
  
  #Check that n is integer
  .n       <- min(ceiling(nsim), nrow(X))
  
  #Create matrix for saving values
  .pifdata <- matrix(data = NA, nrow = .n, ncol = .m)
  
  #Create matrix for saving values 
  .sumdata <- matrix(data = NA, nrow = .m, ncol = 7)
  
  #Get colnames
  colnames(.sumdata) <- c("Simulation", "Min", "First", "Median", "Mean", "Third", "Max")
  
  #Loop reducing i values from the sample and estimating the pif
  for( .j in 1:.m){
    
    #Loop through samples
    for (.i in 1:.n){
      
      #Update X
      if (.j == 1){
        
        .newX <- .X
        .newW <- weights
        
      } else {
        
        #Randomly sample elements to remove
        .todelete <- sample(1:nrow(.X), .j-1, replace = FALSE, prob = weights)
        
        #Remove from this sample 
        .newX           <- data.frame(.X[-.todelete,])
        colnames(.newX) <- colnames(.X)
        .newW <- weights[-.todelete]
        .newW <- .newW/sum(.newW)
      }
      
      
      #Estimate pif and save it matrix
      .pifdata[.i,.j] <- pif(X = .newX, thetahat = thetahat, rr = rr, cft = cft,
                             method = method, weights = .newW, Xvar = NA, adjust = adjust, n = n,
                             ktype = ktype, bw = bw,check_exposure = check_exposure, 
                             check_rr = check_rr, check_integrals = check_integrals,
                             is_paf = is_paf)
      
    }
    
    #Get summary of simulation for each removal .j
    .sumdata[.j,] <- c(.m - (.j-1), summary(.pifdata[,.j]))
    
  }
  
  #Create plot
  .plot <- ggplot(data.frame(.sumdata), aes(x = .m - .sumdata[,"Simulation"])) +
    geom_ribbon(aes(ymin = .sumdata[,"Min"], ymax = .sumdata[,"Max"], 
                    fill = "100% of cases"), alpha = 0.75) +
    geom_ribbon(aes(ymin = .sumdata[,"First"], ymax = .sumdata[,"Third"], 
                    fill = "75% of cases"), alpha = 0.75) +
    geom_line(aes(y  = .sumdata[,"Max"],    color = "100% of cases")) +
    geom_line(aes(y  = .sumdata[,"Min"],    color = "100% of cases")) +
    geom_line(aes(y  = .sumdata[,"First"],  color = "75% of cases")) +
    geom_line(aes(y  = .sumdata[,"Third"],  color = "75% of cases")) +
    geom_line(aes(y  = .sumdata[,"Median"], color = "Median")) +
    geom_point(aes(y = .sumdata[,"Mean"],   color = "Mean")) +
    theme_minimal() + ggtitle(title) + ylab(ylab) + xlab(xlab) + 
    scale_color_manual(name = legendtitle, 
                       values = c("Mean" = colors[1],"Median" = colors[2],
                                  "100% of cases" = colors[3], 
                                  "75% of cases" = colors[4]),
                       drop = FALSE) + 
    scale_fill_manual(values = c("100% of cases" = colors[3], 
                                 "75% of cases" = colors[4]), guide = FALSE) 
  
  return(.plot)
  
}