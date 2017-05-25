#' @title Create a plot of the distribution of exposure under counterfactual 
#'   scenario for continuous exposure
#'   
#' @description Function that creates a plot of the distribution of exposure
#'   \code{X} under counterfactual scenario when exposure is continuous.
#'   
#' @param X      Univariate \code{vector} continuous exposure levels.
#'   
#' @param cft    Counterfactual function of the exposure \code{cft(X)}.
#'   
#'\strong{**Optional**}
#'   
#' @param weights   Normalized survey \code{weights} for the sample \code{X}.
#'   
#' @param dnames    String vector indicating the names of the distributions for 
#'   the legend.
#'   
#' @param ktype    \code{kernel} type:  \code{"gaussian"}, 
#'   \code{"epanechnikov"}, \code{"rectangular"}, \code{"triangular"}, 
#'   \code{"biweight"}, \code{"cosine"}, \code{"optcosine"} (for \code{kernel} 
#'   method). Additional information on kernels in \code{\link[stats]{density}}.
#'   
#' @param bw        Smoothing bandwith parameter from density from
#'   \code{\link[stats]{density}}. Default \code{"SJ"}.
#'   
#' @param adjust    Adjust bandwith parameter from density from
#'   \code{\link[stats]{density}}.
#'   
#' @param n   Number of equally spaced points at which the density is to be
#'   estimated (see \code{\link[stats]{density}}).
#'   
#' @param check_exposure  Check that exposure \code{X} is positive and numeric.
#'   
#' @param title         String with plot title.
#'   
#' @param legendtitle   String title for the legend of plot.
#'   
#' @param xlab          String label for the X-axis of the plot (corresponding 
#'   to "a").
#'   
#' @param ylab          String label for the Y-axis of the plot (corresponding 
#'   to "b").
#'   
#' @param colors        String vector with colors for plots.
#'   
#' @param fill          Boolean that indicates whether there is interior colouring. Default \code{TRUE}.
#'   
#' @param fill_limits   Vector. Limits of subset of the exposure \code{X} such
#'   that only \code{fill_limits[1] < X < fill_limits[2]} are filled with color.
#'   
#' @return cft_plot   \code{\link[ggplot2]{ggplot}} object plotting the shift
#'   from actual to counterfactual distribution.
#'   
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'    
#' @import ggplot2
#' @importFrom stats density approx
#'   
#' @note This function reproduces the classic counterfactual plot from Figure 
#'   25.1 of Vander Hoorn as well as additional plots.
#'   
#' @references Vander Hoorn, S., Ezzati, M., Rodgers, A., Lopez, A. D., & 
#'   Murray, C. J. (2004). \emph{Estimating attributable burden of disease from 
#'   exposure and hazard data. Comparative quantification of health risks: 
#'   global and regional burden of disease attributable to selected major risk 
#'   factors}. Geneva: World Health Organization, 2129-40.
#'   
#' @seealso \code{\link{counterfactual.plot.discrete}} for plotting discrete counterfactuals, 
#'   \code{\link{pif}} for Potential Impact Fraction estimation, 
#'   \code{\link{pif.heatmap}} for sensitivity analysis of the counterfactual, 
#'   \code{\link{pif.plot}} for a plot of potential impact fraction as a 
#'   function of the relative risk's parameter \code{theta}.
#'   
#' @examples
#' 
#' #Example 1: Normal distribution and linear counterfactual
#' #--------------------------------------------------------
#' set.seed(2783569)
#' X   <- data.frame(rnorm(1000, 150, 15))
#' cft <- function(X){0.35*X + 75}  
#' counterfactual.plot.continuous(X, cft, xlab = "Usual SBP (mmHg)", 
#' ylab = "Relative risk of ischaemic heart disease",
#' dnames = c("Current distribution", "Theoretical Minimum Risk Distribution"),
#' title = paste0("Effect of a non-linear hazard function and choice", 
#'                "\nof baseline on total population risk", 
#'                "\n(Fig 25 from Vander Hoorn et al)"))
#'   
#' #Example 2: Counterfactual of BMI reduction only for those with excess-weight (BMI > 25)
#' #--------------------------------------------------------
#' set.seed(2783569)
#' X <- data.frame(Exposure = rlnorm(1000, 3, 0.2))
#' cft <- function(X){
#' 
#'      #Find individuals with excess weight
#'      excess_weight <- which(X[,"Exposure"] > 25)
#'      
#'      #Set those with excess weight to BMI of 25
#'      X[excess_weight, "Exposure"] <- 22.5
#'      
#'      return(X)
#' }     
#' 
#' counterfactual.plot.continuous(X, cft, ktype = "epanechnikov")   
#' 
#' #Change bandwidth method to reduce noice
#' counterfactual.plot.continuous(X, cft, ktype = "epanechnikov", bw = "nrd0")   
#' 
#' #Focus on what happens to the exposure > 23 
#' counterfactual.plot.continuous(X, cft, ktype = "epanechnikov", bw = "nrd0",
#' fill_limits = c(23, Inf)) 
#' 
#' #Delete fill
#' counterfactual.plot.continuous(X, cft, ktype = "epanechnikov", bw = "nrd0", fill = FALSE)   
#'   
#' @keywords internal
#' @export

counterfactual.plot.continuous <- function(X, cft,
                                  weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                                  adjust  = 1, n = 512,
                                  ktype   = c("gaussian", "epanechnikov", "rectangular", "triangular", 
                                             "biweight","cosine", "optcosine"), 
                                  bw      = c("SJ", "nrd0", "nrd", "ucv", "bcv"),
                                  title = "Exposure distribution under current and counterfactual scenarios",
                                  dnames = c("Current distribution", "Counterfactual distribution"),
                                  legendtitle = "Scenario",
                                  xlab = "Exposure", ylab = "Density",
                                  fill_limits = c(-Inf, Inf),
                                  fill = TRUE, 
                                  colors = c("deepskyblue", "tomato3"),
                                  check_exposure = TRUE){
  
  #Set X as data frame
  .X     <- data.frame(X)
  .cX    <- data.frame(cft(.X))
  
  #Get kernel parameters
  .ktype <- as.vector(ktype)[1]
  .bw    <- as.vector(bw)[1]
  
  #Check exposure values are greater than zero
  if(check_exposure){ check.exposure(.X) }
  
  #Create a kernel density plot
  .densX <- density(.X[,1],  kernel = .ktype, n = n, bw = .bw, adjust = adjust)
  .densY <- density(.cX[,1], kernel = .ktype, n = n, bw = .bw, adjust = adjust)
  
  #Densities to data frame
  .dX    <- as.data.frame(.densX[c("x","y")])
  .dY    <- as.data.frame(.densY[c("x","y")])
  
  #Create color vector
  names(colors) <- dnames
  
  #Generate density
  cft_plot <- ggplot() 
  
  #Check if fill subpopulation activated
  if(fill & fill_limits[1] < fill_limits[2] & fill_limits[1] < max(.dX[,1]) & fill_limits[2] > min(.dX[,1])){
    
    #Create subset between fill limits
    .sub_X    <- as.data.frame(subset(.dX$x, .dX[,1] >= fill_limits[1] & .dX[,1] <= fill_limits[2]))
    colnames(.sub_X) <- colnames(X)
    
    .sub_cft  <- cft(.sub_X)
    
    #Get density at new point
    .sub_dX  <- subset(.dX, .dX$x >= min(.sub_X)   & .dX$x <= max(.sub_X))
    .sub_dY  <- subset(.dY, .dY$x >= min(.sub_cft) & .dY$x <= max(.sub_cft))
    
    #Check that .sub_dy is not empty; if empty use approx
    if (nrow(.sub_dY) == 0){
      .sub_dY <- as.data.frame(approx(.dY$x, .dY$y, .sub_cft))
    }
    
    #Add to plot said population
    cft_plot <- cft_plot + 
      geom_ribbon(aes(x = .sub_dX[,1], ymax = .sub_dX[,2], ymin = 0,  
                      fill = dnames[1], color = dnames[1]), linetype = "dashed", alpha = 0.4, data = .sub_dX) +
      geom_ribbon(aes(x = .sub_dY[,1], ymax = .sub_dY[,2], ymin = 0, 
                      fill = dnames[2], color = dnames[2]), linetype = "dashed", alpha = 0.4, data = .sub_dY) +
      scale_fill_manual(name = legendtitle, values = colors)
  }
  
  #Finish plot
  cft_plot <- cft_plot + 
    geom_line(aes(x = .dX[,1], y = .dX[,2], color = dnames[1]), data = .dX) + 
    geom_line(aes(x = .dY[,1], y = .dY[,2], color = dnames[2]), data = .dY) +
    scale_color_manual(name = legendtitle, values = colors) + 
    ggtitle(title) + xlab(xlab) + ylab(ylab) + theme_classic()
  
  
  return(cft_plot)
            
}

