shinyServer(function(input,output){


  ################################################################
  #If PIF needs to be calculated
  ################################################################


  output$PIFselect <- renderUI({
    if(input$PIForPAF>=2){
      selectInput(inputId =  "Cft", label = h5("Select Counterfactual")
                  , choices=list("$$ aX + b $$"= 1, "Other"=2),
                  selected = 1)
    }

  })

  #Load counterfactual
  output$Cftload <- renderUI({
    if(length(input$Cft)>0 && input$PIForPAF>=2){
      if(input$Cft==2){
        fileInput("Cftfun", label = h6("Load Counterfactual function as an R global function"),
                  multiple = FALSE)
      }
    }
  })


  output$aval <- renderUI({withMathJax(
    if(input$PIForPAF>=2 && length(input$Cft)>0){
      if(input$Cft==1){
        textInputRow("a", label = h6("$$a = $$ "),value = 0)
      }
    }
  ) })

  output$bval <- renderUI({withMathJax(
    if(input$PIForPAF>=2 && length(input$Cft)>0){
      if(input$Cft!=2){
        textInputRow("b", label = h6("$$b = $$ "),value = 0)
      }
    }
  )})

  #RR function
  output$RR<-renderUI({(
    if(input$xmin0==1){
      selectInput(inputId =  "rr", label = h4("Choose a Relative Risk function")
                  , choices=list("$$e^{\\theta X}$$"= 1, "$$\\theta_1 X+ 1$$"=2,
                                 "$$\\theta_2 X^2+\\theta_1 X+1$$"=3,
                                 "$$\\theta_3X^3+\\theta_2 X^2+\\theta_1 X+1$$"=4,
                                 "Other"=5),
                  selected = 1)
    }else{
      selectInput(inputId =  "rr", label = h4("Choose a Relative Risk function")
                  , choices=list("$$e^{\\theta_1 X+\\theta_2}+\\theta_0$$"= 1, "$$\\theta_1 X+ \\theta_0$$"=2,
                                 "$$\\theta_2 X^2+\\theta_1 X+\\theta_0$$"=3,
                                 "$$\\theta_3X^3+\\theta_2 X^2+\\theta_1 X+\\theta_0$$"=4,
                                 "Other"=5),
                  selected = 1)
    }
  )})

  #Insert theta values manually


  output$theta <- renderUI({withMathJax(
    if(input$rr == 1 && input$xmin0 == 1){
      textInputRow("theta", label = h6("$$\\theta =$$ "),value = 1)
    })
  })

  output$theta0 <- renderUI({withMathJax(
    if(input$xmin0 == 0 & input$rr < 5){
      textInputRow("theta0", label = h6("$$\\theta_0 =$$ "),value = 1)
    }
  )})

  output$theta1 <- renderUI({ withMathJax(
    if((input$rr > 1 && input$rr < 5) || (input$xmin0 == 0 && input$rr < 5) ){
      textInputRow("theta1", label = h6("$$\\theta_1 =$$ "),value = 1)
    }
  )})

  output$theta2 <- renderUI({ withMathJax(
    if((input$rr >2 && input$rr < 5)||(input$xmin0 == 0 && input$rr == 1)){
      textInputRow("theta2", label = h6("$$\\theta_2 =$$ "),value = 1)
    }
  )
  })
  output$theta3 <- renderUI({
    if(input$rr ==4){
      textInputRow("theta3", label = h6("$$\\theta_3 =$$ "),value = 1)
    }
  })
  ################################################################
  #If Other was chosen
  #Load RR function
  output$RROther <- renderUI({
    if(input$rr==5){
      fileInput("RRfun", label = h6("Load Relative Risk function as an R global function"),
                multiple = FALSE)
    }
  })
  #Introduce theta values
  #Load values for thetahat
  output$LoadThetahat <- renderUI({
    if(input$rr==5){
      fileInput("thetahat", label = h6("Insert  file with the values for $${\\theta}$$ "), multiple = FALSE,
                accept =c("text/csv"))
    }
  })

  #Specify minimum exposure level
  output$Xmin <- renderUI({
    if(input$xmin0==0){
      numericInput("XMin", label = h6("Minimum exposure level"), value = 0, min = 0)
    }
  })

  #Load weight values
  output$weights<-renderUI({
    if(input$w==1){
      fileInput("weight", label = h6("Insert file with the weights for X as .csv"),
                multiple = FALSE, accept = c("text/csv"))
    }
  })


  #######################################
  #Calculate PAF
  ###############################

  #Shows the value for PAF
  output$PAFandPIFout <- renderText({
    rr    <- as.numeric(input$rr)
    if(is.null(input$X)){
      return("Please load file with exposure levels in .csv format")
    }else{
      if(file_ext(input$X$name)!="csv"){
        return("Warning: Incorrect file extension. File with exposure levels must be .csv format")
      }else{
        X   <- as.numeric(as.matrix(read.csv(input$X$datapath,stringsAsFactors = FALSE)))
        if (input$xmin0 ==1) {
          if (rr <5){
            theta <- switch(rr, c(input$theta),
                            c( input$theta1),
                            c( input$theta1, input$theta2),
                            c(input$theta1, input$theta2, input$theta3))
            theta <- as.numeric(theta)
            RR <- switch(rr, function(X,theta){exp(theta[1]*X)},
                         function(X,theta){theta[1]*X+1},
                         function(X,theta){theta[2]*X^2+theta[1]*X+1},
                         function(X,theta){theta[3]*X^3+ theta[2]*X^2+theta[1]*X+1})
            if(theta<0){
              return("$$\\theta$$ values should be greater or equal to zero")
            }else{
              if (input$w==0) {

                PAF <- pif(X=X, thetahat = theta, rr = RR, cft = function(X){0},eval.cvx = FALSE)
                if(input$PIForPAF == 1){
                  #Checar qué onda con multiconvex
                  if(ismulticonvex(f = RR, xmin = 0.99 * min(X), xmax = 1.01 * max(X), tmin = theta-1e-05 * abs(theta) -1e05 ,
                                   tmax = theta + 1e-05 * abs(theta) +1e-05, maxval = 1000)){
                    return(sprintf("PAF = %3f",PAF))
                  }else{
                    return(sprintf("PAF = %3f . Warning: Relative Risk Function not convex",PAF))
                  }
                }
                if(input$PIForPAF > 1){
                  if(input$Cft == 1){
                    Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                    PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft,eval.cvx = FALSE)
                    if(ismulticonvex(f = RR, xmin = 0.99 * min(X), xmax = 1.01 * max(X), tmin = theta-1e-05 * abs(theta) -1e05 ,
                                     tmax = theta + 1e-05 * abs(theta) +1e-05, maxval = 1000)){
                      if(input$PIForPAF == 2){
                        return(sprintf("PIF = %3f",PIF))
                      }else{
                        return(sprintf("PAF = %3f \n PIF = %3f",PAF,PIF))

                      }
                    }else{
                      if(input$PIForPAF == 2){
                        return(sprintf("PIF = %3f . Warning: Relative Risk Function not convex",PIF))
                      }else{
                        return(sprintf("PAF = %3f PIF = %3f . Warning: Relative Risk Function not convex",PAF,PIF))

                      }
                    }
                  }else{
                    if(is.null(input$Cftfun)){
                      return("Please load Counterfactual function in an R file with the function specified as: \n Counterfactual <- function(X){function goes here}. Example: Counterfactual <- function(X){2*X} ")
                    }else{
                      if(file_ext(input$Cftfun$name)!="R"){
                        return("Warning: Incorrect file extension. Counterfactual function must be written in an R file")
                      }else{
                        source(input$Cftfun$datapath)
                        if(exists("Counterfactual") == FALSE){
                          "Warning: There is no function named Counterfactual. The counterfactual function loaded must be named Counterfactual and it must be defined as Counterfactual <- function(X){function goes here} . Example: Counterfactual <- function(X){2*X} "
                        }else{
                          if(is.function(Counterfactual) == FALSE){
                            "Warning: Counterfactual must be a function and must be defined as Counterfactual <- function(X){function goes here} . Example: Counterfactual <- function(X){2*X} "
                          }else{
                            if(is.numeric(try(Counterfactual(0),silent = TRUE)) == FALSE || is.na(Counterfactual(0))){
                              return("Warning: Counterfactual function could not be evaluated.")
                            }else{
                              Cft <- Counterfactual
                              PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft,eval.cvx = FALSE)
                              if(ismulticonvex(f = RR, xmin = 0.99 * min(X), xmax = 1.01 * max(X), tmin = theta-1e-05 * abs(theta) -1e05 ,
                                               tmax = theta + 1e-05 * abs(theta) +1e-05, maxval = 1000)){
                                if(input$PIForPAF == 2){
                                  return(sprintf("PIF = %3f",PIF))
                                }else{
                                  return(sprintf("PAF = %3f \n PIF = %3f",PAF,PIF))

                                }
                              }else{
                                if(input$PIForPAF == 2){
                                  return(sprintf("PIF = %3f . Warning: Relative Risk Function not convex",PIF))
                                }else{
                                  return(sprintf("PAF = %3f PIF = %3f . Warning: Relative Risk Function not convex",PAF,PIF))

                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }

              else{
                #Para x xmin 0 con pesos distintos
                if(is.null(input$weight)){
                  return("Load weights in .csv format")
                }else{
                  if(file_ext(input$weight$name)!="csv"){
                    return("Incorrect format for weights, please load in .csv format")
                  }else{
                    Weight <- as.numeric(as.matrix(read.csv(input$weight$datapath,stringsAsFactors = FALSE)))

                    if(length(Weight)!=length(X)){
                      return(sprintf("Warning:The amount of exposure values and weight values are different, there are
                                    %d exposure values, %d weight values." , length(X), length(Weight)))
                    }else{
                      PAF <- pif(X=X, thetahat = theta, rr = RR, cft = function(X){0},eval.cvx = FALSE, weights = Weight)
                      if(input$PIForPAF==1){
                        #Checar qué onda con multiconvex
                        if(ismulticonvex(f = RR, xmin = 0.99 * min(X), xmax = 1.01 * max(X), tmin = theta-1e-05 * abs(theta) -1e05 ,
                                         tmax = theta + 1e-05 * abs(theta) +1e-05, maxval = 1000)){
                          return(sprintf("PAF = %3f",PAF))
                        }else{
                          return(sprintf("PAF = %3f . Warning: Relative Risk Function not convex",PAF))
                        }
                      }
                      if(input$PIForPAF > 1){
                        if(input$Cft == 1){
                          Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                          PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft, weights =  Weight, eval.cvx = FALSE)
                          if(ismulticonvex(f = RR, xmin = 0.99 * min(X), xmax = 1.01 * max(X), tmin = theta-1e-05 * abs(theta) -1e05 ,
                                           tmax = theta + 1e-05 * abs(theta) +1e-05, maxval = 1000)){
                            if(input$PIForPAF == 2){
                              return(sprintf("PIF = %3f",PIF))
                            }else{
                              return(sprintf("PAF = %3f \n PIF = %3f",PAF,PIF))

                            }
                          }else{
                            if(input$PIForPAF == 2){
                              return(sprintf("PIF = %3f . Warning: Relative Risk Function not convex",PIF))
                            }else{
                              return(sprintf("PAF = %3f PIF = %3f . Warning: Relative Risk Function not convex",PAF,PIF))

                            }
                          }
                        }else{
                          if(is.null(input$Cftfun)){
                            return("Please load Counterfactual function in an R file with the function specified as: \n Counterfactual <- function(X){function goes here}. Example: Counterfactual <- function(X){2*X} ")
                          }else{
                            if(file_ext(input$Cftfun$name)!="R"){
                              return("Warning: Incorrect file extension. Counterfactual function must be written in an R file")
                            }else{
                              source(input$Cftfun$datapath)
                              if(exists("Counterfactual") == FALSE){
                                "Warning: There is no function named Counterfactual. The counterfactual function loaded must be named Counterfactual and it must be defined as Counterfactual <- function(X){function goes here} . Example: Counterfactual <- function(X){2*X} "
                              }else{
                                if(is.function(Counterfactual) == FALSE){
                                  "Warning: Counterfactual must be a function and must be defined as Counterfactual <- function(X){function goes here} . Example: Counterfactual <- function(X){2*X} "
                                }else{
                                  if(is.numeric(try(Counterfactual(0),silent = TRUE)) == FALSE || is.na(Counterfactual(0))){
                                    return("Warning: Counterfactual function could not be evaluated.")
                                  }else{
                                    Cft <- Counterfactual
                                    PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft, weights = Weight, eval.cvx = FALSE)
                                    if(ismulticonvex(f = RR, xmin = 0.99 * min(X), xmax = 1.01 * max(X), tmin = theta-1e-05 * abs(theta) -1e05 ,
                                                     tmax = theta + 1e-05 * abs(theta) +1e-05, maxval = 1000)){
                                      if(input$PIForPAF == 2){
                                        return(sprintf("PIF = %3f",PIF))
                                      }else{
                                        return(sprintf("PAF = %3f \n PIF = %3f",PAF,PIF))

                                      }
                                    }else{
                                      if(input$PIForPAF == 2){
                                        return(sprintf("PIF = %3f . Warning: Relative Risk Function not convex",PIF))
                                      }else{
                                        return(sprintf("PAF = %3f PIF = %3f . Warning: Relative Risk Function not convex",PAF,PIF))

                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }

                    }
                  }
                }
              }
            }
          }else{#para xmin 0 con rr=5
            if(is.null(input$RRfun)){
              return("Please load Relative Risk function in an R file with the function specified as: \n RRfunction <- function(X, theta){function goes here}. Example: RRfunction <- function(X, theta){exp(theta*X)} ")
            }else{
              if(file_ext(input$RRfun$name)!="R"){
                return("Warning: Incorrect file extension. Relative Risk function must be written in an R file")
              }else{
                source(input$RRfun$datapath)
                if(exists("RRfunction") == FALSE){
                  "Warning: There is no function named RRfunction. The Relative Risk function loaded must be named RRfunction and it must be defined as RRfunction <- function(X, theta){function goes here}. Example: RRfunction <- function(X, theta){exp(theta*X)}  "
                }else{
                  if(is.function(RRfunction) == FALSE){
                    "Warning: RRfunction must be a function and must be defined as RRfunction <- function(X, theta){function goes here}. Example: RRfunction <- function(X, theta){exp(theta*X)} "
                  }else{
                    if(is.null(input$thetahat)){
                      return("Please load file with the values for $$\\theta$$")
                    }else{
                      if(file_ext(input$thetahat$name)!="csv"){
                        return("Warning: Incorrect file extension. $$\\theta$$ values must be in a .csv format")
                      }else{
                        theta   <- as.numeric(as.matrix(read.csv(input$thetahat$datapath,stringsAsFactors = FALSE)))
                        if(is.numeric(try(RRfunction(0,theta))) == FALSE || is.na(RRfunction(0,theta))){
                          return("Relative Risk Function could not be evaluated.")
                        }else{
                          if(RRfunction(0,theta) != 1){
                            return("Warning: Relative Risk function at minimum exposure level is not equal to one.")
                          }else{
                            RR <- RRfunction
                            if (input$w==0) {

                              PAF <- pif(X=X, thetahat = theta, rr = RR, cft = function(X){0},eval.cvx = FALSE)
                              if(input$PIForPAF == 1){
                                #Checar qué onda con multiconvex
                                if(ismulticonvex(f = RR, xmin = 0.99 * min(X), xmax = 1.01 * max(X), tmin = theta-1e-05 * abs(theta) -1e05 ,
                                                 tmax = theta + 1e-05 * abs(theta) +1e-05, maxval = 1000)){
                                  return(sprintf("PAF = %3f",PAF))
                                }else{
                                  return(sprintf("PAF = %3f . Warning: Relative Risk Function not convex",PAF))
                                }
                              }
                              if(input$PIForPAF > 1){
                                if(input$Cft == 1){
                                  Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                                  PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft,eval.cvx = FALSE)
                                  if(ismulticonvex(f = RR, xmin = 0.99 * min(X), xmax = 1.01 * max(X), tmin = theta-1e-05 * abs(theta) -1e05 ,
                                                   tmax = theta + 1e-05 * abs(theta) +1e-05, maxval = 1000)){
                                    if(input$PIForPAF == 2){
                                      return(sprintf("PIF = %3f",PIF))
                                    }else{
                                      return(sprintf("PAF = %3f \n PIF = %3f",PAF,PIF))

                                    }
                                  }else{
                                    if(input$PIForPAF == 2){
                                      return(sprintf("PIF = %3f . Warning: Relative Risk Function not convex",PIF))
                                    }else{
                                      return(sprintf("PAF = %3f PIF = %3f . Warning: Relative Risk Function not convex",PAF,PIF))

                                    }
                                  }
                                }else{
                                  if(is.null(input$Cftfun)){
                                    return("Please load Counterfactual function in an R file with the function specified as: \n Counterfactual <- function(X){function goes here}. Example: Counterfactual <- function(X){2*X} ")
                                  }else{
                                    if(file_ext(input$Cftfun$name)!="R"){
                                      return("Warning: Incorrect file extension. Counterfactual function must be written in an R file")
                                    }else{
                                      source(input$Cftfun$datapath)
                                      if(exists("Counterfactual") == FALSE){
                                        "Warning: There is no function named Counterfactual. The counterfactual function loaded must be named Counterfactual and it must be defined as Counterfactual <- function(X){function goes here} . Example: Counterfactual <- function(X){2*X} "
                                      }else{
                                        if(is.function(Counterfactual) == FALSE){
                                          "Warning: Counterfactual must be a function and must be defined as Counterfactual <- function(X){function goes here} . Example: Counterfactual <- function(X){2*X} "
                                        }else{
                                          if(is.numeric(try(Counterfactual(0),silent = TRUE)) == FALSE || is.na(Counterfactual(0))){
                                            return("Warning: Counterfactual function could not be evaluated.")
                                          }else{
                                            Cft <- Counterfactual
                                            PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft,eval.cvx = FALSE)
                                            if(ismulticonvex(f = RR, xmin = 0.99 * min(X), xmax = 1.01 * max(X), tmin = theta-1e-05 * abs(theta) -1e05 ,
                                                             tmax = theta + 1e-05 * abs(theta) +1e-05, maxval = 1000)){
                                              if(input$PIForPAF == 2){
                                                return(sprintf("PIF = %3f",PIF))
                                              }else{
                                                return(sprintf("PAF = %3f \n PIF = %3f",PAF,PIF))

                                              }
                                            }else{
                                              if(input$PIForPAF == 2){
                                                return(sprintf("PIF = %3f . Warning: Relative Risk Function not convex",PIF))
                                              }else{
                                                return(sprintf("PAF = %3f PIF = %3f . Warning: Relative Risk Function not convex",PAF,PIF))

                                              }
                                            }
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }

                            else{
                              #Para x xmin 0 con pesos distintos
                              if(is.null(input$weight)){
                                return("Load weights in .csv format")
                              }else{
                                if(file_ext(input$weight$name)!="csv"){
                                  return("Incorrect format for weights, please load in .csv format")
                                }else{
                                  Weight <- as.numeric(as.matrix(read.csv(input$weight$datapath,stringsAsFactors = FALSE)))

                                  if(length(Weight)!=length(X)){
                                    return(sprintf("Warning:The amount of exposure values and weight values are different, there are
                                                   %d exposure values, %d weight values." , length(X), length(Weight)))
                                  }else{
                                    PAF <- pif(X=X, thetahat = theta, rr = RR, cft = function(X){0},eval.cvx = FALSE, weights = Weight)
                                    if(input$PIForPAF==1){
                                      #Checar qué onda con multiconvex
                                      if(ismulticonvex(f = RR, xmin = 0.99 * min(X), xmax = 1.01 * max(X), tmin = theta-1e-05 * abs(theta) -1e05 ,
                                                       tmax = theta + 1e-05 * abs(theta) +1e-05, maxval = 1000)){
                                        return(sprintf("PAF = %3f",PAF))
                                      }else{
                                        return(sprintf("PAF = %3f . Warning: Relative Risk Function not convex",PAF))
                                      }
                                    }
                                    if(input$PIForPAF > 1){
                                      if(input$Cft == 1){
                                        Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                                        PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft, weights =  Weight, eval.cvx = FALSE)
                                        if(ismulticonvex(f = RR, xmin = 0.99 * min(X), xmax = 1.01 * max(X), tmin = theta-1e-05 * abs(theta) -1e05 ,
                                                         tmax = theta + 1e-05 * abs(theta) +1e-05, maxval = 1000)){
                                          if(input$PIForPAF == 2){
                                            return(sprintf("PIF = %3f",PIF))
                                          }else{
                                            return(sprintf("PAF = %3f \n PIF = %3f",PAF,PIF))

                                          }
                                        }else{
                                          if(input$PIForPAF == 2){
                                            return(sprintf("PIF = %3f . Warning: Relative Risk Function not convex",PIF))
                                          }else{
                                            return(sprintf("PAF = %3f PIF = %3f . Warning: Relative Risk Function not convex",PAF,PIF))

                                          }
                                        }
                                      }else{
                                        if(is.null(input$Cftfun)){
                                          return("Please load Counterfactual function in an R file with the function specified as: \n Counterfactual <- function(X){function goes here}. Example: Counterfactual <- function(X){2*X} ")
                                        }else{
                                          if(file_ext(input$Cftfun$name)!="R"){
                                            return("Warning: Incorrect file extension. Counterfactual function must be written in an R file")
                                          }else{
                                            source(input$Cftfun$datapath)
                                            if(exists("Counterfactual") == FALSE){
                                              "Warning: There is no function named Counterfactual. The counterfactual function loaded must be named Counterfactual and it must be defined as Counterfactual <- function(X){function goes here} . Example: Counterfactual <- function(X){2*X} "
                                            }else{
                                              if(is.function(Counterfactual) == FALSE){
                                                "Warning: Counterfactual must be a function and must be defined as Counterfactual <- function(X){function goes here} . Example: Counterfactual <- function(X){2*X} "
                                              }else{
                                                if(is.numeric(try(Counterfactual(0),silent = TRUE)) == FALSE || is.na(Counterfactual(0))){
                                                  return("Warning: Counterfactual function could not be evaluated.")
                                                }else{
                                                  Cft <- Counterfactual
                                                  PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft, weights = Weight, eval.cvx = FALSE)
                                                  if(ismulticonvex(f = RR, xmin = 0.99 * min(X), xmax = 1.01 * max(X), tmin = theta-1e-05 * abs(theta) -1e05 ,
                                                                   tmax = theta + 1e-05 * abs(theta) +1e-05, maxval = 1000)){
                                                    if(input$PIForPAF == 2){
                                                      return(sprintf("PIF = %3f",PIF))
                                                    }else{
                                                      return(sprintf("PAF = %3f \n PIF = %3f",PAF,PIF))

                                                    }
                                                  }else{
                                                    if(input$PIForPAF == 2){
                                                      return(sprintf("PIF = %3f . Warning: Relative Risk Function not convex",PIF))
                                                    }else{
                                                      return(sprintf("PAF = %3f PIF = %3f . Warning: Relative Risk Function not convex",PAF,PIF))

                                                    }
                                                  }
                                                }
                                              }
                                            }
                                          }
                                        }
                                      }
                                    }

                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }

          }
        }else{#para xmin distinto de 0

        }
      }
    }
  })

})
