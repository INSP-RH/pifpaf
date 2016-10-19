#' @title Potential Impact Fraction Sensitivity Analysis plot
#' 
#' @description Function that plots a sensitivity analysis for the potential impact fraction by checking how estimates vary when reducing the sample
#' 
#' @param X         Random sample (can be vector or matrix) which includes exposure and covariates.
#' 
#' @param thetahat  Estimative of \code{theta} for the Relative Risk function
#' 
#' 
#' @param rr        Function for relative risk
#' 
#' 
#' **Optional**
#' 
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for \code{PAF}
#' 
#' @param weights   Survey \code{weights} for the random sample \code{X}
#' 
#' @param n         Number of samples to include in order to conduct sensitivity analysis
#' 
#' @param m         Limit to number of variables to remove
#' 
#' @param filename  Name of file for saving plot
#' 
#' @param title  Plot title
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' 
#' @import matrixStats ggplot2
#' 
#' @examples 
#' 
#' #Example with risk given by HR
#' set.seed(18427)
#' X <- rgamma(1000, 1,1) #Careful as PIF might not exist for all Gammas
#' thetahat <- 0.12
#' sensitivity.pif(X, thetahat, rr = function(X, theta){exp(theta*X)}, 
#'                 m = 100,  title = "My Sensitivity Analysis")
#' 
#' @export
#' 


sensitivity.pif <- function(X, thetahat, rr, 
                            cft = function(Varx){matrix(0,ncol = ncol(as.matrix(Varx)), nrow = nrow(as.matrix(Varx)))}, 
                            weights = rep(1/nrow(as.matrix(X)),nrow(as.matrix(X))),
                            n = 100, m = 100, filename = NA, 
                            title = "Sensitivity Analysis for Potential Impact Fraction (PIF)"){
  
  #Set X as matrix
  .X       <- as.matrix(X)
  
  #Limit m values
  .m       <- min(ceiling(m), nrow(X))
  
  #Check that n is integer
  .n       <- min(ceiling(n), nrow(X))
  
  #Create matrix for saving values
  .pifdata <- matrix(data = NA, nrow = n, ncol = m)
  
  #Create matrix for saving values 
  .sumdata <- matrix(data = NA, nrow = m, ncol = 7)
  
  #Get colnames
  colnames(.sumdata) <- c("Simulation","Min","First","Median","Mean","Third","Max")
  
  #Loop reducing i values from the sample and estimating the pif
  for( .j in 1:.m){
    
    #Loop through samples
    for (.i in 1:.n){
      
      #Update X
      if (.j == 1){
        
        .newX <- .X
        .newW <- weights
        
      } else {
        
        #Sample
        .todelete <- sample(1:nrow(.X), .j-1, replace = FALSE)
        
        #Remove this sample 
        .newX <- .X[-.todelete,]  
        .newW <- weights[-.todelete]
      }
      
      
      #Estimate pif and save it in pifdata
      .pifdata[.i,.j] <- pif(.newX, thetahat, rr, cft, weights = .newW)
      
    }
    
    #Get minimum and maximum
    .sumdata[.j,] <- c(100 - (.j-1),summary(.pifdata[,.j]))
    
  }
  
  #Create plot
  .plot <- ggplot(as.data.frame(.sumdata),aes(x = .m - .sumdata[,"Simulation"])) +
    geom_errorbar(aes(ymin = .sumdata[,"First"], ymax = .sumdata[,"Third"], color = "75% cases")) +
    geom_line(aes(y =  .sumdata[,"Max"],    color = "Maximum"), linetype = 3) +
    geom_line(aes(y =  .sumdata[,"Min"],    color = "Minimum"), linetype = 3) +
    geom_point(aes(y = .sumdata[,"Mean"],   color = "Mean")) +
    geom_line(aes(y =  .sumdata[,"Median"], color = "Median")) +
    theme_minimal() +
    ggtitle(title) +
    ylab("PIF") +
    xlab("Number of randomly deleted observations for X") + 
    scale_color_manual("Analysis", 
                       values = c("Mean" = "red","Median" = "blue",
                                  "Maximum" = "gray75", "Minimum" = "gray75",
                                  "75% cases" = "gray25"))
  
  if(!is.na(filename)){
    
    ggsave(filename, .plot)
    
  }
  
  return(.plot)
  
}
  