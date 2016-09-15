
shinyUI(fluidPage( theme = "prueba.css",
                       withMathJax(),
                         titlePanel(h1("PAF and PIF" )),
                         #Layout
                         fluidRow(
                           column(3, h3("Relative Risk Function and Counterfactual Function"),
                                  #Select which are going to be calculated
                                  selectInput(inputId =  "PIForPAF", label = h4("Calculate")
                                              , choices=list("PAF"= 1, "PIF"=2,
                                                             "PAF and PIF"=3),
                                              selected = 1),
                                  uiOutput("PIFselect"),
                                  uiOutput("Cftload"),
                                  uiOutput("aval"),
                                  uiOutput("bval"),
                                  #Select a function
                                  uiOutput("RR"),
                                  #Insert theta values
                                  uiOutput("theta"),
                                  uiOutput("theta0"),
                                  uiOutput("theta1"),
                                  uiOutput("theta2"),
                                  uiOutput("theta3"),

                                  #Other RR function

                                  uiOutput("RROther"),
                                  uiOutput("LoadThetahat")),
                           #Load maximum and minimum values for theta








                           ##########################################################################
                           column(4, offset = 1,
                                  tabsetPanel(
                                              tabPanel("Numeric results",uiOutput("PAFandPIFout")),
                                              tabPanel("Plots", 
                                                       h4("Sensitivity Plot"),
                                                       plotOutput("SenstivityPlot")
                                                       ))),
                           column(3, offset = 1,
                                  #Load X values
                                  h3("Exposure levels"),
                                  fileInput("X", label = h4("Load file with exposure levels in .csv format"), multiple = FALSE,
                                            accept =c("text/csv")),

                                  #Minimum exposure value
                                  radioButtons("xmin0", label = h5("Is the minimum exposure level equal to zero?"),
                                               choices = list("Yes"=1,"No"=0), selected = 1, inline = TRUE),
                                  uiOutput("Xmin"),

                                  #Particular weights
                                  radioButtons("w",label = h5("Are there particular weights for each X value?"),
                                               choices = list("Yes"=1,"No"=0),selected = 0, inline = TRUE),

                                  uiOutput("weights"))
                         )
                       )



)
