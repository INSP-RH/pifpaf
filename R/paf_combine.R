#' @title Combine point estimates of PAF of same exposure for exclusive subpopulations
#' 
#' @description Function for fast-computing an overall PAF from subpopulation PAF
#' 
#' @param paf_vector Vector of the \code{paf} for each specific subpopulation. 
#' 
#' @param proportions Vector establishing the proportion of individuals in each subpopulation of \code{paf_vector}.
#' 
#' @return An overall point-estimate of \code{paf} combining all subpopulations. 
#' 
#' @author Rodrigo Zepeda Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho García Formentí \email{daliaf172@gmail.com}
#' 
#' @examples 
#' 
#' #Example 1
#' #-------------
#' pafmen   <- 0.23
#' pafwomen <- 0.17
#' paf.combine(c(pafmen, pafwomen), c(0.45, 0.55))  #total PAF
#' 
#' 
#' @export

paf.combine <- function(paf_vector, proportions){
  
  pif.combine(paf_vector, paf_vector, proportions)
  
}


