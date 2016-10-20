#' @title Text input row
#' 
#' @description Reorganize the display of Input boxes, code comes form answers in http://stackoverflow.com/questions/20637248/shiny-4-small-textinput-boxes-side-by-side
#' 
#' @export

#Text Input
textInputRow<-function (inputId, label, value = "")
{
  div(style="display:inline-block",
      tags$label(label, `for` = inputId),
      tags$input(id = inputId, type = "text", value = value,class="input-small"))
}


