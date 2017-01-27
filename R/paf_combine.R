#' @title Combine point estimates of PAF of same exposure for exclusive 
#'   subpopulations
#'   
#' @description Function for fast-computing an overall PAF from subpopulation 
#'   PAF
#'   
#' @param paf_vector Vector of the \code{\link{paf}} for each specific 
#'   subpopulation.
#'   
#' @param proportions Vector establishing the proportion of individuals in each 
#'   subpopulation of \code{paf_vector}.
#'   
#' @return An overall point-estimate of \code{\link{paf}} combining all 
#'   subpopulations.
#'   
#' @author Rodrigo Zepeda Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho García Formentí \email{daliaf172@gmail.com}
#'   
#' @seealso \code{\link{paf}} for Population Attributable Fraction estimation, 
#'   \code{\link{pif}} for Population Impact Fraction estimation, and 
#'   \code{\link{pif.combine}} for combining several PIF
#'   
#' @examples 
#' 
#' #Example 1
#' #-------------
#' 
#' #Estimate PAF for each subpopulation
#' pafmen   <- paf(X = data.frame(2.7), thetahat = 0.12, 
#'                 rr = function(X, theta){X*theta + 1},
#'                 Xvar = 0.11, method = "approximate")
#' pafwomen <- paf(X = data.frame(3.1), thetahat = 0.12, 
#'                 rr = function(X, theta){exp(X*theta/3)},
#'                 Xvar = 0.17, method = "approximate")
#' 
#' #Combine estimates
#' paf.combine(c(pafmen, pafwomen), c(0.45, 0.55))  
#' 
#' 
#' @export

paf.combine <- function(paf_vector, proportions){
  
  pif.combine(paf_vector, paf_vector, proportions)
  
}


