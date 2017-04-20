#' @title Check csv file was correctly defined in R shiny input
#' 
#' @description Function that verifies the input is not empty, and that the input is in a csv format
#' 
#' @param file  Contains file information from "fileInput"
#' 
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#' 
#' @import tools
#' @export 
check.read.csv <- function(file){
  if(is.null(file)){
    sol <- "warn1"
  }else if(file_ext(file$name)!="csv"){
    sol <- "warn2"
  }else{
    sol <- as.numeric(as.matrix(read.csv(file$datapath,stringsAsFactors = FALSE)))
  }
  return(sol)
}

