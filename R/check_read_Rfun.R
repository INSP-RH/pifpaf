#' @title Check an R function has been given with a specific name
#' 
#' @description Function that verifies the input is not empty. Verifies the input is an R file. Verifies the R file contained a function with a specific name. If everything is correct returns TRUE, if not a warning label.
#' 
#' @param file  Contains file information from "fileInput"
#' @param fname Name of the function as a string
#' 
#' @author Rodrigo Zepeda-Tello \email{rzepeda17@gmail.com}
#' @author Dalia Camacho-García-Formentí \email{daliaf172@gmail.com}
#' 
#' @import tools
#' @export 
check.read.Rfun <- function(file, fname){
  sol <- TRUE
  if(is.null(file)){
    sol <- "warn1"
  }else if(file_ext(file$name)!="R"){
    sol <- "warn3"
  }else{
    source(file$datapath)
    if(exists(fname) == FALSE){
      sol <- "warn4"
    }else if(is.function(eval(parse(text = fname)))==FALSE){
      sol <- "warn5"
    }
  }
  return(sol)
}
