#' @title INTERNAL: Check that expected value of counterfactual decreases exposure levels
#'   
#' @description Function that checks whether the counterfactual function of
#'   \code{\link{pif}} decreases the exposure values \code{X}.
#'   
#' @param cft       Function \code{cft(X)} for counterfactual. Leave empty for 
#'   the Population Attributable Fraction \code{\link{paf}} where counterfactual
#'   is 0 exposure.
#'   
#' @param X         Random sample (data.frame) which includes exposure and
#'   covariates or sample mean.
#'   
#' @author Rodrigo Zepeda Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho García Formentí \email{daliaf172@gmail.com}
#'   
#' @return TRUE if counterfactual \code{cft} is well defined.
#'   
#' @examples 
#' #Example 1 
#' cft <- function(X){0.5*X}
#' X   <- runif(100, 0,2)
#' check.cft(cft, X)
#' 
#' @export


check.cft <- function(cft, X){
  
  #Convert to matrix
  .X <- as.data.frame(X)
  
  #Check rows and columns
  n  <- nrow(.X)
  m  <- ncol(.X)
  
  #Loop checking counterfactual evaluation
  bool <- TRUE
  i    <- 1
  while(i <= n && bool){
    
    #Calculate counterfactual
    Cfti <- cft(.X[i,])
    
    #Loop checking that cft > X
    j <- 1
    while(j <= m && bool){
      if(Cfti[j] > .X[i,j]){
        warning(paste0("Counterfactual function increases some ", 
                       "exposure levels, verifiy it is correct."))
        bool <- FALSE
      }
      j <- j + 1
    }
    i <- i + 1
  }
  return(bool)
}
