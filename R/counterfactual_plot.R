#' @title Create a plot of the distribution of exposure under counterfactual 
#'   scenario
#'   
#' @description Function that creates a plot of the distribution of exposure 
#'   under counterfactual scenario.
#'   
#' @param X      Vector with exposure levels.
#'   
#' @param cft    Counterfactual function of the exposure \code{cft(X)}
#'   
#'   **Optional**
#'   
#' @param weights   Normalized survey \code{weights} for the sample \code{X}.
#'   
#' @param exposure.type Either \code{"continuous"} if distribution is continuous
#'   or \code{"discrete"} if distribution is discrete.
#'   
#' @param ktype    \code{kernel} type for \code{"continuous"} case: 
#'   \code{"gaussian"}, \code{"epanechnikov"}, \code{"rectangular"},
#'   \code{"triangular"}, \code{"biweight"}, \code{"cosine"}, \code{"optcosine"}
#'   (for \code{kernel} method). Additional information on kernels in
#'   \code{\link[stats]{density}}
#'   
#' @param bw        Smoothing bandwith parameter from density (for
#'   \code{"continuous"} case) from \code{\link[stats]{density}}. Default
#'   \code{"SJ"}.
#'   
#' @param adjust    Adjust bandwith parameter from density (for
#'   \code{"continuous"} case) from \code{\link[stats]{density}}.
#'   
#' @param n   Number of equally spaced points at which the density (for 
#'   \code{"continuous"} case) is to be estimated (see 
#'   \code{\link[stats]{density}}).
#'   
#' @param check_exposure  Check that exposure \code{X} is positive and numeric
#'   (if \code{"continuous"})
#'   
#' @param dnames    String vector indicating the names of the distributions for 
#'   the legend
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
#' @param x_axis_order  Order of names in xaxis for plot (\code{"discrete"}
#'   case)
#'   
#' @param fill          Colour the densities? Default \code{TRUE}
#'   
#' @param fill_limits   Vector. Limits of subset of the exposure \code{X} such
#'   that only \code{fill_limits[1] < X < fill_limits[2]} are filled with color.
#'   
#' @return cft_plot   ggplot object plotting the shift from actual to 
#'   counterfactual distribution
#'   
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí \email{daliaf172@gmail.com}
#'   
#' @import ggplot2
#'   
#' @references Vander Hoorn, S., Ezzati, M., Rodgers, A., Lopez, A. D., & 
#'   Murray, C. J. (2004). \emph{Estimating attributable burden of disease from 
#'   exposure and hazard data. Comparative quantification of health risks: 
#'   global and regional burden of disease attributable to selected major risk 
#'   factors}. Geneva: World Health Organization, 2129-40.
#'   
#' @details The function automatically tries to distinguish between
#'   \code{"continuous"} and \code{"discrete"} distribution inputs. By
#'   \code{"continuous"} we meen a vector of real numbers; by \code{"discrete"} 
#'   a vector of strings or factor variables.
#'   
#' @seealso \code{\link{pif}} for Potential Impact Fraction estimation, 
#'   \code{\link{pif.heatmap}} for sensitivity analysis of the counterfactual, 
#'   \code{\link{pif.plot}} for a plot of potential impact fraction as a 
#'   function of theta.
#'   
#' @note This function is a wrapper for
#'   \code{\link{counterfactual.plot.continuous}} and 
#'   \code{\link{counterfactual.plot.discrete}}
#'   
#' @examples
#' 
#' #Example 1: Bivariate exposure
#' #--------------------------------------------------------
#' set.seed(2783569)
#' X   <- sample(c("Exposed","Unexposed"), 100, replace = TRUE, prob = c(0.3, 0.7))
#' cft <- function(X){
#' 
#'      #Find which indivuals are exposed
#'      exposed    <- which(X == "Exposed")
#'      
#'      #Change 1/3 of exposed to unexposed
#'      reduced    <- sample(exposed, length(exposed)/3)
#'      X[reduced] <- "Unexposed"
#'      
#'      return(X)
#' }  
#' counterfactual.plot(X, cft)
#'   
#' #Example 2: Multivariate discrete
#' #--------------------------------------------------------
#' set.seed(2783569)
#' X   <- sample(c("Underweight","Normal","Overweight","Obese"), 1000, 
#'                replace = TRUE, prob = c(0.05, 0.3, 0.25, 0.4))
#'                
#' #Complex counterfactual of changing half of underweights to normal,
#' #1/2 of overweights to normal, 1/3 of obese to normal and 
#' #1/3 of obese to overweight
#' cft <- function(X){
#' 
#'      #Classify the individuals
#'      underweights    <- which(X == "Underweight")
#'      overweights     <- which(X == "Overweight")
#'      obese           <- which(X == "Obese")
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
#'      X[changed_under]   <- "Normal"
#'      X[changed_over]    <- "Normal"
#'      X[obese_to_normal] <- "Normal"
#'      X[obese_to_over]   <- "Overweight"
#'      
#'      return(X)
#' }  
#' 
#' #Create plot of counterfactual distribution
#' cftplot <- counterfactual.plot(X, cft, 
#'                x_axis_order = c("Underweight","Normal","Obese","Overweight")) 
#' cftplot 
#' 
#' #Objects returned are ggplot objects and you can play with them
#' require(ggplot2)
#' cftplot + coord_flip() + theme_grey()
#' 
#' #Example 3: Normal distribution and linear counterfactual
#' #--------------------------------------------------------
#' set.seed(2783569)
#' X   <- rnorm(1000, 150, 15)
#' cft <- function(X){0.35*X + 75}  
#' counterfactual.plot(X, cft, xlab = "Usual SBP (mmHg)", 
#' ylab = "Relative risk of ischaemic heart disease",
#' dnames = c("Current distribution", "Theoretical Minimum Risk Distribution"),
#' title = paste0("Effect of a non-linear hazard function and choice", 
#'                "\nof baseline on total population risk", 
#'                "\n(Fig 25 from Vander Hoorn et al)"))
#'   
#' #Example 4: Counterfactual of BMI reduction only for those with excess-weight (BMI > 25)
#' #--------------------------------------------------------
#' set.seed(2783569)
#' X <- rlnorm(1000, 3, 0.2)
#' cft <- function(X){
#' 
#'      #Find individuals with excess weight
#'      excess_weight <- which(X > 25)
#'      
#'      #Set those with excess weight to BMI of 25
#'      X[excess_weight] <- 25
#'      
#'      return(X)
#' }     
#' 
#' counterfactual.plot(X, cft, ktype = "epanechnikov")   
#' 
#' #Change bandwidth method to reduce noice
#' counterfactual.plot(X, cft, ktype = "epanechnikov", bw = "nrd0")   
#'   
#'   
#' @export

counterfactual.plot <- function(X, cft,
                                weights =  rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))), 
                                adjust  = 1, n = 512,
                                ktype   = c("gaussian", "epanechnikov", "rectangular", "triangular", 
                                            "biweight","cosine", "optcosine"), 
                                bw      = c("SJ", "nrd0", "nrd", "ucv", "bcv"),
                                title = "Exposure distribution under current and counterfactual scenarios",
                                dnames = c("Current distribution", "Counterfactual distribution"),
                                exposure.type = NA,
                                legendtitle = "Scenario",
                                xlab = "Exposure", ylab = "Density",
                                colors = c("deepskyblue", "tomato3"),
                                x_axis_order = unique(X),
                                fill_limits = c(-Inf, Inf),
                                fill = TRUE, 
                                check_exposure = TRUE){
  
  
  #Set X as matrix
  .X    <- as.matrix(X)
  
  #Get exposure type
  .type <- exposure.type
  
  #Identify if exposure is discrete or continuous
  if(is.na(.type)){
    if ( is.factor(.X) || is.character(.X)){
      .type <- "discrete"
    } else if (is.numeric(.X)){
      .type <- "continuous"
    } 
  }
  
  switch(.type,
         continuous = {
           .plot <- counterfactual.plot.continuous(X = .X, cft = cft, weights = weights, adjust = adjust, n = n,
                                                    ktype   = ktype, bw = bw, title = title, dnames = dnames, 
                                                    legendtitle = legendtitle, xlab = xlab, ylab = ylab, 
                                                    colors = colors, fill_limits = fill_limits, fill = fill, 
                                                    check_exposure = check_exposure)
         },
         discrete   = {
           .plot <- counterfactual.plot.discrete(X = .X, cft = cft, weights = weights, title = title, 
                                                  dnames = dnames, legendtitle = legendtitle, xlab = xlab, 
                                                  ylab = ylab, colors = colors, x_axis_order = x_axis_order)
         },
         {
           stop("Could not identify X type please select 'continuous' or 'discrete'")
         })
  
  return(.plot)
}

