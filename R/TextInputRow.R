#' @title textInputRow
#' 
#' @description This function arranges the display of Input boxes,
#' Source: Alex Brown, answer on 719016, "shiny 4 small textInput boxes side-by-side", Stack Overflow, January 15 '14, 
#'         http://stackoverflow.com/a/21132918
#' @param inputId Input that will be used to acces the value 
#' @param label   Display label
#' @param value   Initial value
#' 
#' @export

#Text Input
textInputRow<-function (inputId, label, value = "")
{
  div(style="display:inline-block",
      tags$label(label, `for` = inputId),
      tags$input(id = inputId, type = "text", value = value,class="input-small"))
}
