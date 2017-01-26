#' @title Combine point estimates of PIF of same exposure for exclusive
#'   subpopulations
#'   
#' @description Function for fast-computing an overall \code{\link{pif}} from
#'   subpopulation \code{\link{pif}}s and \code{\link{paf}}s
#'   
#' @param paf_vector Vector of the \code{\link{paf}} for each specific
#'   subpopulation.
#'   
#' @param pif_vector Vector of the \code{\link{pif}} for each specific
#'   subpopulation.
#'   
#' @param proportions Vector establishing the proportion of individuals in each
#'   subpopulation of \code{pif_vector}.
#'   
#' @note To combine \code{pif}s both \code{pif}s and \code{paf}s are required.
#'   
#' @return An overall point-estimate of \code{pif} combining all subpopulations.
#'   
#' @author Rodrigo Zepeda Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho García Formentí \email{daliaf172@gmail.com}
#'   
#' 
#' @seealso \code{\link{paf}} for Population Attributable Fraction estimation, 
#'   \code{\link{pif}} for Population Impact Fraction estimation, and 
#'   \code{\link{paf.combine}} for combining several PAF
#'   
#' @examples 
#' 
#' #Example 1
#' #-------------
#' 
#' #Estimate PAF for each subpopulation
#' pafmen   <- paf(X = 2.7, thetahat = 0.12, rr = function(X, theta){X*theta + 1},
#'                 Xvar = 0.11, method = "approximate")
#' pafwomen <- paf(X = 3.1, thetahat = 0.12, rr = function(X, theta){exp(X*theta/3)},
#'                 Xvar = 0.17, method = "approximate")
#' 
#' #Estimate PIF for each subpopulation
#' pifmen   <- pif(X = 2.7, thetahat = 0.12, rr = function(X, theta){X*theta + 1}, 
#'                 cft = function(X){X/2}, Xvar = 0.11, method = "approximate")
#' pifwomen <- pif(X = 3.1, thetahat = 0.12, rr = function(X, theta){exp(X*theta/3)},
#'                 cft = function(X){X/2}, Xvar = 0.17, method = "approximate")
#' 
#' #Combine estimates
#' pif.combine(c(pifmen, pifwomen), c(pafmen, pafwomen), c(0.45, 0.55))  
#' 
#' 
#' @export

pif.combine <- function(pif_vector, paf_vector, proportions){
  
  #Get length (# of pafs)
  .n <- length(pif_vector)
  
  #Check they have same length
  if (.n != length(proportions) || .n != length(paf_vector)){ 
    warning("pif_vector, paf_vector, and proportions have different length")
  }
  
  #Loop through each paf getting the Expected value of the RR
  .denominator <- 0
  .numerator   <- 0
  for (i in 1:.n){
    
    #Get denominator of pif
    .rr          <- 1/(1 - paf_vector[i])
    .denominator <- .denominator + .rr*proportions[i]
    
    #Get numerator of pif
    .cft         <- (1 - pif_vector[i])*.rr 
    .numerator   <- .numerator + .cft*proportions[i]
  }
  
  .pif <- 1 - .numerator/.denominator
  
  return(.pif)
  
}


