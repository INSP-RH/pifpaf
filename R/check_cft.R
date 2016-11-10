#' @title Check counterfactual decreases exposure levels
#' 
#' @description Function that checks whether the counterfactual function decreases the exposure values
#' 
#' @param cft   counterfactual function
#' 
#' @param X     sample of exposure values / or mean exposure, when only point estimates are available
#' 
#' 
#' @examples 
#' #Example 1 
#' cft <- function(X){0.5*X}
#' X   <- runif(100, 0,2)
#' check.cft(cft, X)
#' 
#' 
#' @export

check.cft <- function(cft, X){
  .X <- as.matrix(X)
  n  <- dim(.X)[1]
  m  <- dim(.X)[2]
  bool <- TRUE
  i    <- 1
  while(i <= n && bool){
    Cfti <- cft(.X[i,])
    
    j <- 1
    while(j <= m && bool){
      if(Cfti[j] > .X[i,j]){
        warning("Counterfactual function increases exposure some exposure levels, verifiy it is correct.")
        bool <- FALSE
      }
      j <- j + 1
    }
    i <- i+1
  }
}
