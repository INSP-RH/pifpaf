#' @title Sum of entry multiplication
#' 
#' @description Sum of entry multiplication.
#' 
#' @param X matrix
#' 
#' @param Y matrix
#' 
#' @author Rodrigo Zepeda Tello \email{rodrigo.zepeda@insp.mx}
#' @author Dalia Camacho García Formentí
#' 
#'@example 
#'
#' X   <- matrix(c(1,2,3,4,5,6), ncol = 3)
#' Y   <- matrix(2, ncol = 3, nrow = 2)
#'EntryMult(X,Y)
#'  
#' @export

EntryMult <- function(X,Y){
  if(dim(X)[1]==dim(Y)[1] && dim(X)[2]==dim(Y)[2]){
    m   <- dim(X)[1]
    n   <- dim(X)[2]
    Sum <- 0
    for(i in 1:m){
      for(j in 1:n){
        Sum <- Sum + X[i,j]*Y[i,j]
      }
    }
    return(Sum)
  }else{
    warning("Dimensions don't agree")
  }
}