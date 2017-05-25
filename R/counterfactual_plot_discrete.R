#' @title Create a plot of the distribution of exposure under counterfactual
#'   scenario for discrete exposure.
#'   
#' @description Function that creates a plot of the distribution of exposure
#'   \code{X} under counterfactual scenario \code{cft} for discrete exposures.
#' 
#' @param X      One dimensional \code{vector} with exposure levels.
#'   
#' @param cft    Counterfactual function of the exposure \code{cft(X)}.
#'   
#'   \bold{**Optional**}
#'   
#' @param weights   Normalized survey \code{weights} for the sample \code{X}.
#'   
#' @param dnames    String vector indicating the names of the distributions for
#'   the legend.
#'   
#' @param title     String with plot title.
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
#' @param x_axis_order  Order of names in xaxis for plot.
#'   
#' @return cft_plot   \code{\link[ggplot2]{ggplot}} object plotting the shift from actual to
#'   counterfactual distribution.
#'   
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#'  
#'   
#' @import ggplot2
#'   
#' @references Vander Hoorn, S., Ezzati, M., Rodgers, A., Lopez, A. D., & 
#'   Murray, C. J. (2004). \emph{Estimating attributable burden of disease from 
#'   exposure and hazard data. Comparative quantification of health risks: 
#'   global and regional burden of disease attributable to selected major risk 
#'   factors}. Geneva: World Health Organization, 2129-40.
#'   
#' @seealso \code{\link{counterfactual.plot.continuous}} for plotting continuous counterfactuals, 
#'   \code{\link{pif}} for Potential Impact Fraction estimation, 
#'   \code{\link{pif.heatmap}} for sensitivity analysis of the counterfactual, 
#'   \code{\link{pif.plot}} for a plot of potential impact fraction as a 
#'   function of the relative risk's parameter \code{theta}.
#'   
#' @examples
#' 
#' #Example 1: Bivariate exposure
#' #--------------------------------------------------------
#' set.seed(2783569)
#' X   <- data.frame(Exposure = 
#'                    sample(c("Exposed","Unexposed"), 100, 
#'                    replace = TRUE, prob = c(0.3, 0.7)))
#' cft <- function(X){
#' 
#'      #Find which indivuals are exposed
#'      exposed      <- which(X[,"Exposure"] == "Exposed")
#'      
#'      #Change 1/3 of exposed to unexposed
#'      reduced                 <- sample(exposed, length(exposed)/3)
#'      X[reduced,"Exposure"]   <- "Unexposed"
#'      
#'      return(X)
#' }  
#' counterfactual.plot.discrete(X, cft)
#'   
#' #Example 2: Multivariate exposure
#' #--------------------------------------------------------
#' set.seed(2783569)
#' X   <- data.frame(
#'          Exposure = sample(c("Underweight","Normal","Overweight","Obese"), 
#'          1000, replace = TRUE, prob = c(0.05, 0.3, 0.25, 0.4)))
#'                
#' #Complex counterfactual of changing half of underweights to normal,
#' #1/2 of overweights to normal, 1/3 of obese to normal and 
#' #1/3 of obese to overweight
#' cft <- function(X){
#' 
#'      #Classify the individuals
#'      underweights    <- which(X[,"Exposure"] == "Underweight")
#'      overweights     <- which(X[,"Exposure"] == "Overweight")
#'      obese           <- which(X[,"Exposure"] == "Obese")
#'      
#'      #Sample 1/2 underweights and overweights and 2/3 of obese
#'      changed_under    <- sample(underweights, length(underweights)/2)
#'      changed_over     <- sample(overweights,  length(overweights)/2)
#'      changed_obese    <- sample(obese,        2*length(obese)/3)
#'      
#'      #Assign those obese that go to normal and those that go to overweight
#'      obese_to_normal  <- sample(changed_obese, length(changed_obese)/2)
#'      obese_to_over    <- which(!(changed_obese %in% obese_to_normal))
#'      
#'      #Change the individuals to normal and overweight
#'      X[changed_under,"Exposure"]   <- "Normal"
#'      X[changed_over,"Exposure"]    <- "Normal"
#'      X[obese_to_normal,"Exposure"] <- "Normal"
#'      X[obese_to_over,"Exposure"]   <- "Overweight"
#'      
#'      return(X)
#' }  
#' 
#' #Create plot of counterfactual distribution
#' counterfactual.plot.discrete(X, cft, x_axis_order = c("Underweight","Normal","Obese","Overweight")) 
#' 
#' @keywords internal 
#' @export

counterfactual.plot.discrete <- function(X, cft,
                                weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                                title = "Exposure distribution under current and counterfactual scenarios",
                                dnames = c("Current distribution", "Counterfactual distribution"),
                                legendtitle = "Scenario",
                                xlab = "Exposure", ylab = "Frequency",
                                colors = c("deepskyblue", "tomato3"),
                                x_axis_order = unique(X[,1])){
  
  #Set X as matrix
  .X  <- as.data.frame(X)
  colnames(.X) <- colnames(X)
  .cX <- data.frame(cft(.X), rep(dnames[2], nrow(.X)))
  .X  <- data.frame(.X, rep(dnames[1], nrow(.X)))
  
  colnames(.X) <- c(colnames(X)[1],"Distribution")
  colnames(.cX) <- c(colnames(X)[1],"Distribution")
  
  #Create a kernel density plot
  .dens_data         <- as.data.frame(rbind(.X,.cX))
  .dens_data$weights <- rep(weights, 2)
  
  #Create color vector
  names(colors) <- dnames
  
  ggplot(.dens_data) + 
    geom_bar(aes(x = .dens_data[,1], fill = .dens_data[,2], weight = .dens_data[,3]), 
             position="dodge", stat = "count") + 
    scale_fill_manual(name = legendtitle, values = colors) + 
    ggtitle(title) + xlab(xlab) + ylab(ylab) + theme_classic() +
    scale_x_discrete(limits = x_axis_order)
  
}

