shinyUI(fluidPage( theme = "prueba2.css",
                   withMathJax(),
                   titlePanel(h1("Population Attributable Fraction and Population Impact Fraction" )),
                   #Layout
                   
                   sidebarLayout(
                     sidebarPanel( selectInput(inputId =  "PIForPAF", label = h3("Calculate")
                                               , choices=list("PAF"= 1, "PIF"=2,
                                                              "PAF and PIF"=3),
                                               selected = 1),
                                  #Select which are going to be calculated
                                  
                                 
                                  
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
                                  uiOutput("LoadThetahat"),
                                  
                                  uiOutput("PIFselect"),
                                  uiOutput("Cftload"),
                                  uiOutput("aval"),
                                  uiOutput("bval"),

                                  #Load X values
                                  h4("Exposure levels"),
                                  fileInput("X", label = h5("Load file with exposure levels in .csv format"), multiple = FALSE,
                                            accept =c("text/csv")),
                                  
                                  #Minimum exposure value
                                  radioButtons("xmin0", label = h5("Is the minimum exposure level equal to zero?"),
                                               choices = list("Yes"=1,"No"=0), selected = 1, inline = TRUE),
                                  uiOutput("Xmin"),
                                  
                                  #Particular weights
                                  radioButtons("w",label = h5("Are there particular weights for each X value?"),
                                               choices = list("Yes"=1,"No"=0),selected = 0, inline = TRUE),
                                  
                                  uiOutput("weights")),
                     #Load maximum and minimum values for theta
                     
                     
                     
                     
                     
                     
                     
                     
                     ##########################################################################
                     mainPanel(
                       tabsetPanel(
                         tabPanel("$$\\text{Numeric results}$$", uiOutput("PAFandPIFout")),
                         tabPanel("$$\\text{Sensitivity  Plot}$$", 
                                  h4("Sensitivity Plot"),
                                  plotOutput("SenstivityPlot")
                         ),
                         tabPanel("$$\\text{Plot for different } \\theta \\text{ values (Univariate)}$$",
                                  textInputRow("minTheta", label = h6("$$ \\text{Minimum } \\theta = $$  "),value = 0),
                                  textInputRow("maxTheta", label = h6("$$\\text{Maximum } \\theta =$$ "),value = 1),
                                  plotOutput("PIFPlot")),
                         tabPanel("$$\\text{Confidence Intervals}$$",
                                  textInputRow("Conf", label = h6("$$ \\text{Confidence Level } = $$  "),value = 95),
                                  uiOutput("sdtheta"),
                                  uiOutput("sdMat"),
                                  uiOutput("IntervalPAF"))
                       ))
                   )
                   
                   
)



)
