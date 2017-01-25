#' @title Combine point estimates of PIF of same exposure for exclusive subpopulations
#' 
#' @description Function for fast-computing an overall PIF from subpopulation PIFs and PAFs
#' 
#' @param paf_vector Vector of the \code{paf} for each specific subpopulation. 
#' 
#' @param pif_vector Vector of the \code{pif} for each specific subpopulation. 
#' 
#' @param proportions Vector establishing the proportion of individuals in each subpopulation of \code{paf_vector}.
#' 
#' @note To combine \code{pif}s both \code{pif}s and \code{paf}s are required.
#' 
#' @return An overall point-estimate of \code{pif} combining all subpopulations. 
#' 
#' @author Rodrigo Zepeda Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho García Formentí \email{daliaf172@gmail.com}
#' 
#' @examples 
#' 
#' #Example 1
#' #-------------
#' pif   <- c("Men" = 0.23, "Women" = 0.19)
#' paf   <- c("Men" = 0.47, "Women" = 0.51)
#' props <- c("Men" = 0.45, "Women" = 0.55)
#' pif.combine(pif, paf, props)  #total PAF
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


