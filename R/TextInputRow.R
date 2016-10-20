#' @title Text input row
#' 
#' @description This function arranges the display of Input boxes,
#'  code comes form answers in http://stackoverflow.com/questions/20637248/shiny-4-small-textinput-boxes-side-by-side .
#' 
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