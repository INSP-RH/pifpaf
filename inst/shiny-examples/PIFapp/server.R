
shinyServer(function(input,output){
  
  
  
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
      selectInput(inputId =  "rr", label = h5("Choose a Relative Risk function")
                  , choices=list("$$e^{\\theta X}$$"= 1, "$$\\theta_1 X+ 1$$"=2,
                                 "$$\\theta_2 X^2+\\theta_1 X+1$$"=3,
                                 "$$\\theta_3X^3+\\theta_2 X^2+\\theta_1 X+1$$"=4,
                                 "Other"=5),
                  selected = 1)
    }else{
      selectInput(inputId =  "rr", label = h5("Choose a Relative Risk function")
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
      fileInput("thetahat", label = h6("$$ \\text{Load  file with the values of } \\theta$$ "), multiple = FALSE,
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
  
  #Reactive sd levels
  
  output$sdtheta <- renderUI({withMathJax(
    if(input$rr <= 2 && input$xmin0 == 1){
      textInputRow("SDtheta", label = h6("$$\\text{Standard deviation of } \\theta =$$ "),value = 1)
    })
  })
  
  output$sdMat <- renderUI({withMathJax(
    if((input$rr > 2 && input$xmin0 == 1)){
      fileInput("SDMat", label = h6("$$\\text{Load .csv file with the square root of the variance and covariance matrix of  } \\theta $$ "), 
                multiple = FALSE, accept = c("text/csv"))
    })
  })
  
  #######################################
  #Calculate PAF
  ###############################
  
  #Shows the value for PAF
  output$PAFandPIFout <- renderText({
    rr    <- as.numeric(input$rr)
    if(is.null(input$X)){
      Error <- 1
      return(Warnings1(Error))
    }else{
      if(file_ext(input$X$name)!="csv"){
        Error <- 2
        return(Warnings1(Error))
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
              Error <- 3
              return(Warnings1(Error))
            }else{
              if (input$w==0) {
                PAF <- pif(X=X, thetahat = theta, rr = RR, eval.cvx = FALSE)
                
                if(input$PIForPAF == 1){
                  return(sprintf("PAF = %3f",PAF))
                }
                if(input$PIForPAF > 1){
                  if(input$Cft == 1){
                    Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                    PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft,eval.cvx = FALSE)
                    
                    if(input$PIForPAF == 2){
                      return(sprintf("PIF = %3f",PIF))
                    }else{
                      return(sprintf("PAF = %3f <br/> PIF = %3f",PAF,PIF))
                    }
                  }else{
                    if(is.null(input$Cftfun)){
                      Error <- 4
                      return(Warnings1(Error))
                    }else{
                      if(file_ext(input$Cftfun$name)!="R"){
                        Error <- 5
                        return(Warnings1(Error))
                      }else{
                        source(input$Cftfun$datapath)
                        if(exists("Counterfactual") == FALSE){
                          Error <- 6
                          return(Warnings1(Error))
                        }else{
                          if(is.function(Counterfactual) == FALSE){
                            Error <- 7
                            return(Warnings1(Error))
                          }else{
                            if(is.numeric(try(Counterfactual(0),silent = TRUE)) == FALSE || is.na(Counterfactual(0))){
                              Error <- 8
                              return(Warnings1(Error))
                            }else{
                              Cft <- Counterfactual
                              PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft,eval.cvx = FALSE)
                              if(input$PIForPAF == 2){
                                return(sprintf("PIF = %3f",PIF))
                              }else{
                                return(sprintf("PAF = %3f <br/> PIF = %3f",PAF,PIF))
                                
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
                  Error <- 9
                  return(Warnings1(Error))
                }else{
                  if(file_ext(input$weight$name)!="csv"){
                    Error <-10
                    return(Warnings1(Error))
                  }else{
                    Weight <- as.numeric(as.matrix(read.csv(input$weight$datapath,stringsAsFactors = FALSE)))
                    
                    if(length(Weight)!=length(X)){
                      return(sprintf("Warning: The amount of exposure values and weight values are different, there are %d exposure values, %d weight values." , length(X), length(Weight)))
                    }else{
                      PAF <- pif(X=X, thetahat = theta, rr = RR, eval.cvx = FALSE, weights = Weight)
                      if(input$PIForPAF==1){
                        return(sprintf("PAF = %3f",PAF))
                      }
                      if(input$PIForPAF > 1){
                        if(input$Cft == 1){
                          Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                          PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft, weights =  Weight, eval.cvx = FALSE)
                          
                          if(input$PIForPAF == 2){
                            return(sprintf("PIF = %3f",PIF))
                          }else{
                            return(sprintf("PAF = %3f <br/> PIF = %3f",PAF,PIF))
                            
                          }
                          
                        }else{
                          if(is.null(input$Cftfun)){
                            Error <- 4
                            return(Warnings1(Error))
                          }else{
                            if(file_ext(input$Cftfun$name)!="R"){
                              Error <- 5
                              return(Warnings1(Error))
                            }else{
                              source(input$Cftfun$datapath)
                              if(exists("Counterfactual") == FALSE){
                                Error <- 6
                                return(Warnings1(Error))     
                              }else{
                                if(is.function(Counterfactual) == FALSE){
                                  Error <- 7
                                  return(Warnings1(Error))                             
                                }else{
                                  if(is.numeric(try(Counterfactual(0),silent = TRUE)) == FALSE || is.na(Counterfactual(0))){
                                    Error <- 8
                                    return(Warnings1(Error)) 
                                  }else{
                                    Cft <- Counterfactual
                                    PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft, weights = Weight, eval.cvx = FALSE)
                                    if(input$PIForPAF == 2){
                                      return(sprintf("PIF = %3f",PIF))
                                    }else{
                                      return(sprintf("PAF = %3f <br/> PIF = %3f",PAF,PIF))
                                      
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
              Error <- 11
              return(Warnings1(Error))
            }else{
              if(file_ext(input$RRfun$name)!="R"){
                Error <- 12
                return(Warnings1(Error))
              }else{
                source(input$RRfun$datapath)
                if(exists("RRfunction") == FALSE){
                  Error <- 13
                  return(Warnings1(Error))
                }else{
                  if(is.function(RRfunction) == FALSE){
                    Error <- 14
                    return(Warnings1(Error))
                  }else{
                    if(is.null(input$thetahat)){
                      Error <- 15
                      return(Warnings1(Error))
                    }else{
                      if(file_ext(input$thetahat$name)!="csv"){
                        Error <- 16
                        return(Warnings1(Error))
                      }else{
                        theta <- as.numeric(as.matrix(read.csv(input$thetahat$datapath,stringsAsFactors = FALSE)))
                        if(is.numeric(try(RRfunction(0,theta))) == FALSE || is.na(RRfunction(0,theta))){
                          Error <- 17
                          return(Warnings1(Error))
                        }else{
                          if(RRfunction(0,theta) != 1){
                            Error <- 18
                            return(Warnings1(Error))
                          }else{
                            RR <- RRfunction
                            if (input$w==0) {
                              PAF <- pif(X=X, thetahat = theta, rr = RR, eval.cvx = FALSE)
                              if(input$PIForPAF == 1){
                                return(sprintf("PAF = %3f",PAF))
                              }
                              if(input$PIForPAF > 1){
                                if(input$Cft == 1){
                                  Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                                  PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft,eval.cvx = FALSE)
                                  if(input$PIForPAF == 2){
                                    return(sprintf("PIF = %3f",PIF))
                                  }else{
                                    return(sprintf("PAF = %3f <br/> PIF = %3f",PAF,PIF))
                                  }
                                  
                                }else{
                                  if(is.null(input$Cftfun)){
                                    Error <- 4
                                    return(Warnings1(Error))
                                  }else{
                                    if(file_ext(input$Cftfun$name)!="R"){
                                      Error <- 5
                                      return(Warnings1(Error))
                                    }else{
                                      source(input$Cftfun$datapath)
                                      if(exists("Counterfactual") == FALSE){
                                        Error <- 6
                                        return(Warnings1(Error))     
                                      }else{
                                        if(is.function(Counterfactual) == FALSE){
                                          Error <- 7
                                          return(Warnings1(Error))    
                                        }else{
                                          if(is.numeric(try(Counterfactual(0),silent = TRUE)) == FALSE || is.na(Counterfactual(0))){
                                            Error <- 8
                                            return(Warnings1(Error))
                                          }else{
                                            Cft <- Counterfactual
                                            PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft,eval.cvx = FALSE)
                                            if(input$PIForPAF == 2){
                                              return(sprintf("PIF = %3f",PIF))
                                            }else{
                                              return(sprintf("PAF = %3f <br/> PIF = %3f",PAF,PIF))
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
                                Error <- 9
                                return(Warnings1(Error))
                              }else{
                                if(file_ext(input$weight$name)!="csv"){
                                  Error <- 10
                                  return(Warnings1(Error))
                                }else{
                                  Weight <- as.numeric(as.matrix(read.csv(input$weight$datapath,stringsAsFactors = FALSE)))
                                  
                                  if(length(Weight)!=length(X)){
                                    return(sprintf("Warning: The amount of exposure values and weight values are different, there are %d exposure values, %d weight values." , length(X), length(Weight)))
                                  }else{
                                    PAF <- pif(X=X, thetahat = theta, rr = RR, eval.cvx = FALSE, weights = Weight)
                                    if(input$PIForPAF==1){
                                      return(sprintf("PAF = %3f",PAF))
                                      
                                    }
                                    if(input$PIForPAF > 1){
                                      if(input$Cft == 1){
                                        Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                                        PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft, weights =  Weight, eval.cvx = FALSE)
                                        if(input$PIForPAF == 2){
                                          return(sprintf("PIF = %3f",PIF))
                                        }else{
                                          return(sprintf("PAF = %3f <br/> PIF = %3f",PAF,PIF))
                                          
                                        }
                                        
                                      }else{
                                        if(is.null(input$Cftfun)){
                                          Error <- 4
                                          return(Warnings1(Error))
                                        }else{
                                          if(file_ext(input$Cftfun$name)!="R"){
                                            Error <- 5
                                            return(Warnings1(Error))
                                          }else{
                                            source(input$Cftfun$datapath)
                                            if(exists("Counterfactual") == FALSE){
                                              Error <- 6
                                              return(Warnings1(Error))
                                            }else{
                                              if(is.function(Counterfactual) == FALSE){
                                                Error <- 7
                                                return(Warnings1(Error))
                                              }else{
                                                if(is.numeric(try(Counterfactual(0),silent = TRUE)) == FALSE || is.na(Counterfactual(0))){
                                                  Error <- 8
                                                  return(Warnings1(Error))
                                                }else{
                                                  Cft <- Counterfactual
                                                  PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft, weights = Weight, eval.cvx = FALSE)
                                                  if(input$PIForPAF == 2){
                                                    return(sprintf("PIF = %3f",PIF))
                                                  }else{
                                                    return(sprintf("PAF = %3f <br/> PIF = %3f",PAF,PIF))
                                                    
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
          xmin <- input$XMin
          if(xmin<0){
            Error <- 19
            return(Warnings1(Error))
          }else{
            if(rr < 5){
              theta <- switch(rr, c(input$theta0, input$theta1, input$theta2),
                              c( input$theta0, input$theta1),
                              c( input$theta0, input$theta1, input$theta2),
                              c(input$theta0, input$theta1, input$theta2, input$theta3))
              theta <- as.numeric(theta)
              RR <- switch(rr, function(X,theta){exp(theta[2]*X+theta[3])+theta[1]},
                           function(X,theta){theta[2]*X+theta[1]},
                           function(X,theta){theta[3]*X^2+theta[2]*X+theta[1]},
                           function(X,theta){theta[4]*X^3+ theta[3]*X^2+theta[2]*X+theta[1]})
              if(RR(xmin,theta)!=1){
                Error <- 20
                return(Warnings1(Error))
              }else{
                if (input$w==0) {
                  
                  PAF <- pif(X=X, thetahat = theta, rr = RR, cft = function(X){sapply(X,function(x){xmin})},eval.cvx = FALSE)
                  if(input$PIForPAF == 1){
                    return(sprintf("PAF = %3f",PAF))
                  }
                  if(input$PIForPAF > 1){
                    if(input$Cft == 1){
                      Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                      PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft,eval.cvx = FALSE)
                      if(input$PIForPAF == 2){
                        return(sprintf("PIF = %3f",PIF))
                      }else{
                        return(sprintf("PAF = %3f <br/> PIF = %3f",PAF,PIF))
                      }
                      
                    }else{
                      if(is.null(input$Cftfun)){
                        Error <- 4
                        return(Warnings1(Error))      
                      }else{
                        if(file_ext(input$Cftfun$name)!="R"){
                          Error <- 5
                          return(Warnings1(Error))
                        }else{
                          source(input$Cftfun$datapath)
                          if(exists("Counterfactual") == FALSE){
                            Error <- 6
                            return(Warnings1(Error))
                          }else{
                            if(is.function(Counterfactual) == FALSE){
                              Error <- 7
                              return(Warnings1(Error))
                            }else{
                              if(is.numeric(try(Counterfactual(xmin),silent = TRUE)) == FALSE || is.na(Counterfactual(xmin))){
                                Error <- 8
                                return(Warnings1(Error))
                              }else{
                                Cft <- Counterfactual
                                PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft,eval.cvx = FALSE)
                                if(input$PIForPAF == 2){
                                  return(sprintf("PIF = %3f",PIF))
                                }else{
                                  return(sprintf("PAF = %3f <br/> PIF = %3f",PAF,PIF))
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }else{
                  if(is.null(input$weight)){
                    Error <- 9
                    return(Warnings1(Error))
                  }else{
                    if(file_ext(input$weight$name)!="csv"){
                      Error <- 10
                      return(Warnings1(Error))
                    }else{
                      Weight <- as.numeric(as.matrix(read.csv(input$weight$datapath,stringsAsFactors = FALSE)))
                      
                      if(length(Weight)!=length(X)){
                        return(sprintf("Warning: The amount of exposure values and weight values are different, there are %d exposure values, %d weight values." , length(X), length(Weight)))
                      }else{
                        PAF <- pif(X=X, thetahat = theta, rr = RR, cft = function(X){sapply(X,function(x){xmin})},eval.cvx = FALSE, weights = Weight)
                        if(input$PIForPAF==1){
                          return(sprintf("PAF = %3f",PAF))
                        }
                        if(input$PIForPAF > 1){
                          if(input$Cft == 1){
                            Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                            PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft, weights =  Weight, eval.cvx = FALSE)
                            
                            if(input$PIForPAF == 2){
                              return(sprintf("PIF = %3f",PIF))
                            }else{
                              return(sprintf("PAF = %3f <br/> PIF = %3f",PAF,PIF))
                            }
                            
                          }else{
                            if(is.null(input$Cftfun)){
                              Error <- 4
                              return(Warnings1(Error))
                            }else{
                              if(file_ext(input$Cftfun$name)!="R"){
                                Error <- 5
                                return(Warnings1(Error))
                              }else{
                                source(input$Cftfun$datapath)
                                if(exists("Counterfactual") == FALSE){
                                  Error <- 6
                                  return(Warnings1(Error))
                                }else{
                                  if(is.function(Counterfactual) == FALSE){
                                    Error <- 7
                                    return(Warnings1(Error))
                                  }else{
                                    if(is.numeric(try(Counterfactual(xmin),silent = TRUE)) == FALSE || is.na(Counterfactual(xmin))){
                                      Error <- 8
                                      return(Warnings1(Error))
                                    }else{
                                      Cft <- Counterfactual
                                      PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft, weights = Weight, eval.cvx = FALSE)
                                      if(input$PIForPAF == 2){
                                        return(sprintf("PIF = %3f",PIF))
                                      }else{
                                        return(sprintf("PAF = %3f <br/> PIF = %3f",PAF,PIF))
                                        
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
            }else{
              if(is.null(input$RRfun)){
                Error <- 11
                return(Warnings1(Error))
              }else{
                if(file_ext(input$RRfun$name)!="R"){
                  Error <- 12
                  return(Warnings1(Error))
                }else{
                  source(input$RRfun$datapath)
                  if(exists("RRfunction") == FALSE){
                    Error <- 13
                    return(Warnings1(Error))
                  }else{
                    if(is.function(RRfunction) == FALSE){
                      Error <- 14
                      return(Warnings1(Error))
                    }else{
                      if(is.null(input$thetahat)){
                        Error <- 15
                        return(Warnings1(Error))
                      }else{
                        if(file_ext(input$thetahat$name)!="csv"){
                          Error <- 16
                          return(Warnings1(Error))
                        }else{
                          theta   <- as.numeric(as.matrix(read.csv(input$thetahat$datapath,stringsAsFactors = FALSE)))
                          if(is.numeric(try(RRfunction(xmin,theta))) == FALSE || is.na(RRfunction(xmin,theta))){
                            Error <- 17
                            return(Warnings1(Error))
                          }else{
                            if(RRfunction(xmin,theta) != 1){
                              Error <- 18
                              return(Warnings1(Error))
                            }else{
                              RR <- RRfunction
                              if (input$w==0) {
                                
                                PAF <- pif(X=X, thetahat = theta, rr = RR, cft = function(X){sapply(X,function(x){xmin})},eval.cvx = FALSE)
                                if(input$PIForPAF == 1){
                                  return(sprintf("PAF = %3f",PAF))
                                  
                                }
                                if(input$PIForPAF > 1){
                                  if(input$Cft == 1){
                                    Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                                    PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft,eval.cvx = FALSE)
                                    
                                    if(input$PIForPAF == 2){
                                      return(sprintf("PIF = %3f",PIF))
                                    }else{
                                      return(sprintf("PAF = %3f <br/> PIF = %3f",PAF,PIF))
                                      
                                    }
                                  }else{
                                    if(is.null(input$Cftfun)){
                                      Error <- 4
                                      return(Warnings1(Error))
                                    }else{
                                      if(file_ext(input$Cftfun$name)!="R"){
                                        Error <- 5
                                        return(Warnings1(Error))
                                      }else{
                                        source(input$Cftfun$datapath)
                                        if(exists("Counterfactual") == FALSE){
                                          Error <- 6
                                          return(Warnings1(Error))
                                        }else{
                                          if(is.function(Counterfactual) == FALSE){
                                            Error <- 7
                                            return(Warnings1(Error))
                                          }else{
                                            if(is.numeric(try(Counterfactual(xmin),silent = TRUE)) == FALSE || is.na(Counterfactual(xmin))){
                                              Error <- 8
                                              return(Warnings1(Error))
                                            }else{
                                              Cft <- Counterfactual
                                              PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft,eval.cvx = FALSE)
                                              if(input$PIForPAF == 2){
                                                return(sprintf("PIF = %3f",PIF))
                                              }else{
                                                return(sprintf("PAF = %3f <br/> PIF = %3f",PAF,PIF))
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
                                  Error <- 9
                                  return(Warnings1(Error))
                                }else{
                                  if(file_ext(input$weight$name)!="csv"){
                                    Error <- 10
                                    return(Warnings1(Error))
                                  }else{
                                    Weight <- as.numeric(as.matrix(read.csv(input$weight$datapath,stringsAsFactors = FALSE)))
                                    
                                    if(length(Weight)!=length(X)){
                                      return(sprintf("Warning: The amount of exposure values and weight values are different, there are %d exposure values, %d weight values." , length(X), length(Weight)))
                                    }else{
                                      PAF <- pif(X=X, thetahat = theta, rr = RR, cft = function(X){sapply(X,function(x){xmin})},eval.cvx = FALSE, weights = Weight)
                                      if(input$PIForPAF==1){
                                        return(sprintf("PAF = %3f",PAF))
                                      }
                                      if(input$PIForPAF > 1){
                                        if(input$Cft == 1){
                                          Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                                          PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft, weights =  Weight, eval.cvx = FALSE)
                                          if(input$PIForPAF == 2){
                                            return(sprintf("PIF = %3f",PIF))
                                          }else{
                                            return(sprintf("PAF = %3f <br/> PIF = %3f",PAF,PIF))
                                            
                                          }
                                        }else{
                                          if(is.null(input$Cftfun)){
                                            Error <- 4
                                            return(Warnings1(Error))
                                          }else{
                                            if(file_ext(input$Cftfun$name)!="R"){
                                              Error <- 5
                                              return(Warnings1(Error))
                                            }else{
                                              source(input$Cftfun$datapath)
                                              if(exists("Counterfactual") == FALSE){
                                                Error <- 6
                                                return(Warnings1(Error))
                                              }else{
                                                if(is.function(Counterfactual) == FALSE){
                                                  Error <- 7
                                                  return(Warnings1(Error))
                                                }else{
                                                  if(is.numeric(try(Counterfactual(xmin),silent = TRUE)) == FALSE || is.na(Counterfactual(xmin))){
                                                    Error <- 8
                                                    return(Warnings1(Error))
                                                  }else{
                                                    Cft <- Counterfactual
                                                    PIF <- pif(X=X, thetahat = theta, rr = RR, cft = Cft, weights = Weight, eval.cvx = FALSE)
                                                    if(input$PIForPAF == 2){
                                                      return(sprintf("PIF = %3f",PIF))
                                                    }else{
                                                      return(sprintf("PAF = %3f <br/> PIF = %3f",PAF,PIF))
                                                      
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
        }
      }
    }
  })
  
  
  ################################
  #Sensitivity Plot
  output$SenstivityPlot <- renderPlot({
    rr    <- as.numeric(input$rr)
    if(is.null(input$X)){
      Error <- 1; return(Warnings2(Error));
    }else{
      if(file_ext(input$X$name)!="csv"){
        Error <- 2; return(Warnings2(Error));
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
              Error <- 3; return(Warnings2(Error));
            }else{
              if (input$w==0) {
                
                PAF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR,  title = "Sensitivity Analysis for Population Attributable Fraction \n (PAF)")
                if(input$PIForPAF == 1){
                  return(PAF.Plot)
                }
                if(input$PIForPAF > 1){
                  if(input$Cft == 1){
                    Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                    PIF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR, cft = Cft)
                    
                    if(input$PIForPAF == 2){
                      return(PIF.Plot)
                    }else{
                      return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
                    }
                    
                  }else{
                    if(is.null(input$Cftfun)){
                      Error <- 4; return(Warnings2(Error));
                    }else{
                      if(file_ext(input$Cftfun$name)!="R"){
                        Error <- 5; return(Warnings2(Error));
                      }else{
                        source(input$Cftfun$datapath)
                        if(exists("Counterfactual") == FALSE){
                          Error <- 6; return(Warnings2(Error));
                        }else{
                          if(is.function(Counterfactual) == FALSE){
                            Error <- 7; return(Warnings2(Error));
                          }else{
                            if(is.numeric(try(Counterfactual(0),silent = TRUE)) == FALSE || is.na(Counterfactual(0))){
                              Error <- 8; return(Warnings2(Error));
                            }else{
                              Cft <- Counterfactual
                              PIF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR, cft = Cft)
                              if(input$PIForPAF == 2){
                                return(PIF.Plot)
                              }else{
                                return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
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
                  Error <- 9; return(Warnings2(Error));
                }else{
                  if(file_ext(input$weight$name)!="csv"){
                    Error <- 10; return(Warnings2(Error));
                  }else{
                    Weight <- as.numeric(as.matrix(read.csv(input$weight$datapath,stringsAsFactors = FALSE)))
                    
                    if(length(Weight)!=length(X)){
                      return(ggplot()+ggtitle(sprintf("Warning: The amount of exposure values and weight values \n  are different, there are \n %d exposure values, %d weight values." , length(X), length(Weight))))
                    }else{
                      PAF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR,  weights = Weight, title = "Sensitivity Analysis for Population Attributable Fraction \n (PAF)")
                      if(input$PIForPAF==1){
                        if(input$PIForPAF == 1){
                          return(PAF.Plot)
                        }
                      }
                      if(input$PIForPAF > 1){
                        if(input$Cft == 1){
                          Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                          PIF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR, cft = Cft, weights =  Weight)
                          if(input$PIForPAF == 2){
                            return(PIF.Plot)
                          }else{
                            return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
                          }
                        }else{
                          if(is.null(input$Cftfun)){
                            Error <- 4; return(Warnings2(Error));
                          }else{
                            if(file_ext(input$Cftfun$name)!="R"){
                              Error <- 5; return(Warnings2(Error));
                            }else{
                              source(input$Cftfun$datapath)
                              if(exists("Counterfactual") == FALSE){
                                Error <- 6; return(Warnings2(Error));
                              }else{
                                if(is.function(Counterfactual) == FALSE){
                                  Error <- 7; return(Warnings2(Error));
                                }else{
                                  if(is.numeric(try(Counterfactual(0),silent = TRUE)) == FALSE || is.na(Counterfactual(0))){
                                    Error <- 8; return(Warnings2(Error));
                                  }else{
                                    Cft <- Counterfactual
                                    PIF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR, cft = Cft, weights = Weight)
                                    if(input$PIForPAF == 2){
                                      return(PIF.Plot)
                                    }else{
                                      return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
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
              Error <- 11; return(Warnings2(Error));
            }else{
              if(file_ext(input$RRfun$name)!="R"){
                Error <- 12; return(Warnings2(Error));
              }else{
                source(input$RRfun$datapath)
                if(exists("RRfunction") == FALSE){
                  Error <- 13; return(Warnings2(Error));
                }else{
                  if(is.function(RRfunction) == FALSE){
                    Error <- 14; return(Warnings2(Error));
                  }else{
                    if(is.null(input$thetahat)){
                      Error <- 15; return(Warnings2(Error));
                    }else{
                      if(file_ext(input$thetahat$name)!="csv"){
                        Error <- 16; return(Warnings2(Error));
                      }else{
                        theta   <- as.numeric(as.matrix(read.csv(input$thetahat$datapath,stringsAsFactors = FALSE)))
                        if(is.numeric(try(RRfunction(0,theta))) == FALSE || is.na(RRfunction(0,theta))){
                          Error <- 17; return(Warnings2(Error));
                        }else{
                          if(RRfunction(0,theta) != 1){
                            Error <- 18; return(Warnings2(Error));
                          }else{
                            RR <- RRfunction
                            if (input$w==0) {
                              
                              PAF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR,  title = "Sensitivity Analysis for Population Attributable Fraction \n (PAF)")
                              if(input$PIForPAF == 1){
                                if(input$PIForPAF == 1){
                                  return(PAF.Plot)
                                }
                              }
                              if(input$PIForPAF > 1){
                                if(input$Cft == 1){
                                  Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                                  PIF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR, cft = Cft)
                                  if(input$PIForPAF == 2){
                                    return(PIF.Plot)
                                  }else{
                                    return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
                                  }
                                }else{
                                  if(is.null(input$Cftfun)){
                                    Error <- 4; return(Warnings2(Error));
                                  }else{
                                    if(file_ext(input$Cftfun$name)!="R"){
                                      Error <- 5; return(Warnings2(Error));
                                    }else{
                                      source(input$Cftfun$datapath)
                                      if(exists("Counterfactual") == FALSE){
                                        Error <- 6; return(Warnings2(Error));
                                      }else{
                                        if(is.function(Counterfactual) == FALSE){
                                          Error <- 7; return(Warnings2(Error));                                        }else{
                                            if(is.numeric(try(Counterfactual(0),silent = TRUE)) == FALSE || is.na(Counterfactual(0))){
                                              Error <- 8; return(Warnings2(Error));
                                            }else{
                                              Cft <- Counterfactual
                                              PIF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR, cft = Cft)
                                              if(input$PIForPAF == 2){
                                                return(PIF.Plot)
                                              }else{
                                                return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
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
                                Error <- 9; return(Warnings2(Error));
                              }else{
                                if(file_ext(input$weight$name)!="csv"){
                                  Error <- 10; return(Warnings2(Error));
                                }else{
                                  Weight <- as.numeric(as.matrix(read.csv(input$weight$datapath,stringsAsFactors = FALSE)))
                                  
                                  if(length(Weight)!=length(X)){
                                    return(ggplot()+ggtitle(sprintf("Warning: The amount of exposure values and weight values \n are different, there are \n  %d exposure values, %d weight values." , length(X), length(Weight))))
                                  }else{
                                    PAF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR,  weights = Weight, title = "Sensitivity Analysis for Population Attributable Fraction \n (PAF)")
                                    if(input$PIForPAF==1){
                                      if(input$PIForPAF == 1){
                                        return(PAF.Plot)
                                      }
                                    }
                                    if(input$PIForPAF > 1){
                                      if(input$Cft == 1){
                                        Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                                        PIF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR, cft = Cft, weights =  Weight)
                                        if(input$PIForPAF == 2){
                                          return(PIF.Plot)
                                        }else{
                                          return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
                                        }
                                      }else{
                                        if(is.null(input$Cftfun)){
                                          Error <- 4; return(Warnings2(Error));
                                        }else{
                                          if(file_ext(input$Cftfun$name)!="R"){
                                            Error <- 5; return(Warnings2(Error));
                                          }else{
                                            source(input$Cftfun$datapath)
                                            if(exists("Counterfactual") == FALSE){
                                              Error <- 6; return(Warnings2(Error));
                                            }else{
                                              if(is.function(Counterfactual) == FALSE){
                                                Error <- 7; return(Warnings2(Error));
                                              }else{
                                                if(is.numeric(try(Counterfactual(0),silent = TRUE)) == FALSE || is.na(Counterfactual(0))){
                                                  Error <- 8; return(Warnings2(Error));
                                                }else{
                                                  Cft <- Counterfactual
                                                  PIF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR, cft = Cft, weights = Weight)
                                                  if(input$PIForPAF == 2){
                                                    return(PIF.Plot)
                                                  }else{
                                                    return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
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
          xmin <- input$XMin
          if(xmin<0){
            Error <- 19; return(Warnings2(Error));
          }else{
            if(rr < 5){
              theta <- switch(rr, c(input$theta0, input$theta1, input$theta2),
                              c( input$theta0, input$theta1),
                              c( input$theta0, input$theta1, input$theta2),
                              c(input$theta0, input$theta1, input$theta2, input$theta3))
              theta <- as.numeric(theta)
              RR <- switch(rr, function(X,theta){exp(theta[2]*X+theta[3])+theta[1]},
                           function(X,theta){theta[2]*X+theta[1]},
                           function(X,theta){theta[3]*X^2+theta[2]*X+theta[1]},
                           function(X,theta){theta[4]*X^3+ theta[3]*X^2+theta[2]*X+theta[1]})
              if(RR(xmin,theta)!=1){
                Error <- 18; return(Warnings2(Error));
              }else{
                if (input$w==0) {
                  
                  PAF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR, cft = function(X){sapply(X,function(x){xmin})}, title = "Sensitivity Analysis for Population Attributable Fraction \n (PAF)")
                  if(input$PIForPAF == 1){
                    if(input$PIForPAF == 1){
                      return(PAF.Plot)
                    }
                  }
                  if(input$PIForPAF > 1){
                    if(input$Cft == 1){
                      Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                      PIF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR, cft = Cft)
                      if(input$PIForPAF == 2){
                        return(PIF.Plot)
                      }else{
                        return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
                      }
                    }else{
                      if(is.null(input$Cftfun)){
                        Error <- 4; return(Warnings2(Error));
                      }else{
                        if(file_ext(input$Cftfun$name)!="R"){
                          Error <- 5; return(Warnings2(Error));
                        }else{
                          source(input$Cftfun$datapath)
                          if(exists("Counterfactual") == FALSE){
                            Error <- 6; return(Warnings2(Error));
                          }else{
                            if(is.function(Counterfactual) == FALSE){
                              Error <- 7; return(Warnings2(Error));
                            }else{
                              if(is.numeric(try(Counterfactual(xmin),silent = TRUE)) == FALSE || is.na(Counterfactual(xmin))){
                                Error <- 8; return(Warnings2(Error));
                              }else{
                                Cft <- Counterfactual
                                PIF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR, cft = Cft)
                                if(input$PIForPAF == 2){
                                  return(PIF.Plot)
                                }else{
                                  return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }else{
                  if(is.null(input$weight)){
                    Error <- 9; return(Warnings2(Error));
                  }else{
                    if(file_ext(input$weight$name)!="csv"){
                      Error <- 10; return(Warnings2(Error));
                    }else{
                      Weight <- as.numeric(as.matrix(read.csv(input$weight$datapath,stringsAsFactors = FALSE)))
                      
                      if(length(Weight)!=length(X)){
                        return(ggplot()+ggtitle(sprintf("Warning: The amount of exposure values and weight values \n are different, there are \n  %d exposure values, %d weight values." , length(X), length(Weight))))
                      }else{
                        PAF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR, cft = function(X){sapply(X,function(x){xmin})}, weights = Weight, title = "Sensitivity Analysis for Population Attributable Fraction \n (PAF)")
                        if(input$PIForPAF==1){
                          if(input$PIForPAF == 1){
                            return(PAF.Plot)
                          }
                        }
                        if(input$PIForPAF > 1){
                          if(input$Cft == 1){
                            Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                            PIF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR, cft = Cft, weights =  Weight)
                            if(input$PIForPAF == 2){
                              return(PIF.Plot)
                            }else{
                              return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
                            }
                          }else{
                            if(is.null(input$Cftfun)){
                              Error <- 4; return(Warnings2(Error));
                            }else{
                              if(file_ext(input$Cftfun$name)!="R"){
                                Error <- 5; return(Warnings2(Error));
                              }else{
                                source(input$Cftfun$datapath)
                                if(exists("Counterfactual") == FALSE){
                                  Error <- 6; return(Warnings2(Error));
                                }else{
                                  if(is.function(Counterfactual) == FALSE){
                                    Error <- 7; return(Warnings2(Error));
                                  }else{
                                    if(is.numeric(try(Counterfactual(xmin),silent = TRUE)) == FALSE || is.na(Counterfactual(xmin))){
                                      Error <- 8; return(Warnings2(Error));
                                    }else{
                                      Cft <- Counterfactual
                                      PIF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR, cft = Cft, weights = Weight)
                                      if(input$PIForPAF == 2){
                                        return(PIF.Plot)
                                      }else{
                                        return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
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
            }else{
              if(is.null(input$RRfun)){
                Error <- 11; return(Warnings2(Error));
              }else{
                if(file_ext(input$RRfun$name)!="R"){
                  Error <- 12; return(Warnings2(Error));
                }else{
                  source(input$RRfun$datapath)
                  if(exists("RRfunction") == FALSE){
                    Error <- 13; return(Warnings2(Error));
                  }else{
                    if(is.function(RRfunction) == FALSE){
                      Error <- 14; return(Warnings2(Error));
                    }else{
                      if(is.null(input$thetahat)){
                        Error <- 15; return(Warnings2(Error));
                      }else{
                        if(file_ext(input$thetahat$name)!="csv"){
                          Error <- 16; return(Warnings2(Error));
                        }else{
                          theta   <- as.numeric(as.matrix(read.csv(input$thetahat$datapath,stringsAsFactors = FALSE)))
                          if(is.numeric(try(RRfunction(xmin,theta))) == FALSE || is.na(RRfunction(xmin,theta))){
                            Error <- 17; return(Warnings2(Error));
                          }else{
                            if(RRfunction(xmin,theta) != 1){
                              Error <- 18; return(Warnings2(Error));
                            }else{
                              RR <- RRfunction
                              if (input$w==0) {
                                
                                PAF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR, cft = function(X){sapply(X,function(x){xmin})}, title = "Sensitivity Analysis for Population Attributable Fraction \n (PAF)")
                                if(input$PIForPAF == 1){
                                  if(input$PIForPAF == 1){
                                    return(PAF.Plot)
                                  }
                                }
                                if(input$PIForPAF > 1){
                                  if(input$Cft == 1){
                                    Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                                    PIF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR, cft = Cft)
                                    if(input$PIForPAF == 2){
                                      return(PIF.Plot)
                                    }else{
                                      return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
                                    }
                                  }else{
                                    if(is.null(input$Cftfun)){
                                      Error <- 4; return(Warnings2(Error));
                                    }else{
                                      if(file_ext(input$Cftfun$name)!="R"){
                                        Error <- 5; return(Warnings2(Error));
                                      }else{
                                        source(input$Cftfun$datapath)
                                        if(exists("Counterfactual") == FALSE){
                                          Error <- 6; return(Warnings2(Error));
                                        }else{
                                          if(is.function(Counterfactual) == FALSE){
                                            Error <- 7; return(Warnings2(Error));
                                          }else{
                                            if(is.numeric(try(Counterfactual(xmin),silent = TRUE)) == FALSE || is.na(Counterfactual(xmin))){
                                              Error <- 8; return(Warnings2(Error));
                                            }else{
                                              Cft <- Counterfactual
                                              PIF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR, cft = Cft)
                                              if(input$PIForPAF == 2){
                                                return(PIF.Plot)
                                              }else{
                                                return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
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
                                  Error <- 9; return(Warnings2(Error));
                                }else{
                                  if(file_ext(input$weight$name)!="csv"){
                                    Error <- 10; return(Warnings2(Error));
                                  }else{
                                    Weight <- as.numeric(as.matrix(read.csv(input$weight$datapath,stringsAsFactors = FALSE)))
                                    
                                    if(length(Weight)!=length(X)){
                                      return(ggplot()+ggtitle(sprintf("Warning: The amount of exposure values and weight values \n are different, there are \n  %d exposure values, %d weight values." , length(X), length(Weight))))
                                    }else{
                                      PAF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR, cft = function(X){sapply(X,function(x){xmin})}, weights = Weight, title = "Sensitivity Analysis for Population Attributable Fraction \n (PAF)")
                                      if(input$PIForPAF==1){
                                        if(input$PIForPAF == 1){
                                          return(PAF.Plot)
                                        }
                                      }
                                      if(input$PIForPAF > 1){
                                        if(input$Cft == 1){
                                          Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                                          PIF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR, cft = Cft, weights =  Weight)
                                          if(input$PIForPAF == 2){
                                            return(PIF.Plot)
                                          }else{
                                            return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
                                          }
                                        }else{
                                          if(is.null(input$Cftfun)){
                                            Error <- 4; return(Warnings2(Error));
                                          }else{
                                            if(file_ext(input$Cftfun$name)!="R"){
                                              Error <- 5; return(Warnings2(Error));
                                            }else{
                                              source(input$Cftfun$datapath)
                                              if(exists("Counterfactual") == FALSE){
                                                Error <- 6; return(Warnings2(Error));
                                              }else{
                                                if(is.function(Counterfactual) == FALSE){
                                                  Error <- 7; return(Warnings2(Error));
                                                }else{
                                                  if(is.numeric(try(Counterfactual(xmin),silent = TRUE)) == FALSE || is.na(Counterfactual(xmin))){
                                                    Error <- 8; return(Warnings2(Error));
                                                  }else{
                                                    Cft <- Counterfactual
                                                    PIF.Plot <- sensitivity.pif(X=X, thetahat = theta, rr = RR, cft = Cft, weights = Weight)
                                                    if(input$PIForPAF == 2){
                                                      return(PIF.Plot)
                                                    }else{
                                                      return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
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
        }
      }
    }
  })
  
  ################################################
  #Plot for different theta values
  ################################################
  
  
  output$PIFPlot <- renderPlot({
    thetamin <- as.numeric(input$minTheta)
    thetamax <- as.numeric(input$maxTheta)
    
    if(thetamin>=thetamax){
      Error <- 20; return(Warnings2(Error));
    }else{
      rr <- as.numeric(input$rr)
      if(input$xmin0 == 1){
        if(input$rr >2 && input$rr != 5){
          Error <- 21; return(Warnings2(Error));
        }else{
          
          if(is.null(input$X)){
            Error <- 1; return(Warnings2(Error));
          }else{
            if(file_ext(input$X$name)!="csv"){
              Error <- 2; return(Warnings2(Error));
            }else{
              X   <- as.numeric(as.matrix(read.csv(input$X$datapath,stringsAsFactors = FALSE)))
              if(rr < 5){
                RR <- switch(rr, function(X,theta){exp(theta[1]*X)},
                             function(X,theta){theta[1]*X+1})
                if(thetamin<0){
                  Error <- 3; return(Warnings2(Error));
                }else{
                  if (input$w==0) {
                    
                    PAF.Plot <- plotpif(X=X, thetamin = thetamin, thetamax = thetamax, rr = RR, )+ggtitle("Population Attributable Fraction (PAF) \n under different values of theta")
                    if(input$PIForPAF == 1){
                      return(PAF.Plot)
                    }
                    if(input$PIForPAF > 1){
                      if(input$Cft == 1){
                        Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                        PIF.Plot <- plotpif(X=X, thetamin = thetamin, thetamax = thetamax, rr = RR, cft = Cft)
                        
                        if(input$PIForPAF == 2){
                          return(PIF.Plot)
                        }else{
                          return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
                        }
                        
                      }else{
                        if(is.null(input$Cftfun)){
                          Error <- 4; return(Warnings2(Error));
                        }else{
                          if(file_ext(input$Cftfun$name)!="R"){
                            Error <- 5; return(Warnings2(Error));
                          }else{
                            source(input$Cftfun$datapath)
                            if(exists("Counterfactual") == FALSE){
                              Error <- 6; return(Warnings2(Error));
                            }else{
                              if(is.function(Counterfactual) == FALSE){
                                Error <- 7; return(Warnings2(Error));
                              }else{
                                if(is.numeric(try(Counterfactual(0),silent = TRUE)) == FALSE || is.na(Counterfactual(0))){
                                  Error <- 8; return(Warnings2(Error));
                                }else{
                                  Cft <- Counterfactual
                                  PIF.Plot <- plotpif(X=X, thetamin = thetamin, thetamax = thetamax, rr = RR, cft = Cft)
                                  if(input$PIForPAF == 2){
                                    return(PIF.Plot)
                                  }else{
                                    return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
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
                      Error <- 9; return(Warnings2(Error));
                    }else{
                      if(file_ext(input$weight$name)!="csv"){
                        Error <- 10; return(Warnings2(Error));
                      }else{
                        Weight <- as.numeric(as.matrix(read.csv(input$weight$datapath,stringsAsFactors = FALSE)))
                        
                        if(length(Weight)!=length(X)){
                          return(ggplot()+ggtitle(sprintf("Warning: The amount of exposure values and weight values \n  are different, there are \n %d exposure values, %d weight values." , length(X), length(Weight))))
                        }else{
                          PAF.Plot <- plotpif(X=X, thetamin = thetamin, thetamax = thetamax, rr = RR,  weights = Weight) + ggtitle("Population Attributable Fraction (PAF) \n under different values of theta")
                          if(input$PIForPAF==1){
                            if(input$PIForPAF == 1){
                              return(PAF.Plot)
                            }
                          }
                          if(input$PIForPAF > 1){
                            if(input$Cft == 1){
                              Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                              PIF.Plot <- plotpif(X=X, thetamin = thetamin, thetamax = thetamax, rr = RR, cft = Cft, weights =  Weight)
                              if(input$PIForPAF == 2){
                                return(PIF.Plot)
                              }else{
                                return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
                              }
                            }else{
                              if(is.null(input$Cftfun)){
                                Error <- 4; return(Warnings2(Error));
                              }else{
                                if(file_ext(input$Cftfun$name)!="R"){
                                  Error <- 5; return(Warnings2(Error));
                                }else{
                                  source(input$Cftfun$datapath)
                                  if(exists("Counterfactual") == FALSE){
                                    Error <- 6; return(Warnings2(Error));
                                  }else{
                                    if(is.function(Counterfactual) == FALSE){
                                      Error <- 7; return(Warnings2(Error));
                                    }else{
                                      if(is.numeric(try(Counterfactual(0),silent = TRUE)) == FALSE || is.na(Counterfactual(0))){
                                        Error <- 8; return(Warnings2(Error));
                                      }else{
                                        Cft <- Counterfactual
                                        PIF.Plot <- plotpif(X=X, thetamin = thetamin, thetamax = thetamax, rr = RR, cft = Cft, weights = Weight)
                                        if(input$PIForPAF == 2){
                                          return(PIF.Plot)
                                        }else{
                                          return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
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
              }else{#rr==5
                if(is.null(input$RRfun)){
                  Error <- 11; return(Warnings2(Error));
                }else{
                  if(file_ext(input$RRfun$name)!="R"){
                    Error <- 12; return(Warnings2(Error));
                  }else{
                    source(input$RRfun$datapath)
                    if(exists("RRfunction") == FALSE){
                      Error <- 13; return(Warnings2(Error));
                    }else{
                      if(is.function(RRfunction) == FALSE){
                        Error <- 14; return(Warnings2(Error));
                      }else{
                        
                        if(is.numeric(try(RRfunction(0,thetamin))) == FALSE || is.na(RRfunction(0,thetamin))){
                          Error <- 17; return(Warnings2(Error));
                        }else{
                          if(RRfunction(0,thetamin) != 1 && RRfunction(0,thetamax) != 1){
                            Error <- 22; return(Warnings2(Error));
                          }else{
                            RR <- RRfunction
                            if (input$w==0) {
                              
                              PAF.Plot <- plotpif(X=X, thetamin = thetamin, thetamax = thetamax, rr = RR, ) + ggtitle("Population Attributable Fraction (PAF) \n under different values of theta")
                              if(input$PIForPAF == 1){
                                if(input$PIForPAF == 1){
                                  return(PAF.Plot)
                                }
                              }
                              if(input$PIForPAF > 1){
                                if(input$Cft == 1){
                                  Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                                  PIF.Plot <- plotpif(X=X, thetamin = thetamin, thetamax = thetamax, rr = RR, cft = Cft)
                                  if(input$PIForPAF == 2){
                                    return(PIF.Plot)
                                  }else{
                                    return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
                                  }
                                }else{
                                  if(is.null(input$Cftfun)){
                                    Error <- 4; return(Warnings2(Error));
                                  }else{
                                    if(file_ext(input$Cftfun$name)!="R"){
                                      Error <- 5; return(Warnings2(Error));
                                    }else{
                                      source(input$Cftfun$datapath)
                                      if(exists("Counterfactual") == FALSE){
                                        Error <- 6; return(Warnings2(Error));
                                      }else{
                                        if(is.function(Counterfactual) == FALSE){
                                          Error <- 7; return(Warnings2(Error));
                                        }else{
                                          if(is.numeric(try(Counterfactual(0),silent = TRUE)) == FALSE || is.na(Counterfactual(0))){
                                            Error <- 8; return(Warnings2(Error));
                                          }else{
                                            Cft <- Counterfactual
                                            PIF.Plot <- plotpif(X=X, thetamin = thetamin, thetamax = thetamax, rr = RR, cft = Cft)
                                            if(input$PIForPAF == 2){
                                              return(PIF.Plot)
                                            }else{
                                              return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
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
                                Error <- 9; return(Warnings2(Error));
                              }else{
                                if(file_ext(input$weight$name)!="csv"){
                                  Error <- 10; return(Warnings2(Error));
                                }else{
                                  Weight <- as.numeric(as.matrix(read.csv(input$weight$datapath,stringsAsFactors = FALSE)))
                                  
                                  if(length(Weight)!=length(X)){
                                    return(ggplot()+ggtitle(sprintf("Warning: The amount of exposure values and weight values \n are different, there are \n  %d exposure values, %d weight values." , length(X), length(Weight))))
                                  }else{
                                    PAF.Plot <- plotpif(X=X, thetamin = thetamin, thetamax = thetamax, rr = RR,  weights = Weight) + ggtitle("Population Attributable Fraction (PAF) \n under different values of theta")
                                    if(input$PIForPAF==1){
                                      if(input$PIForPAF == 1){
                                        return(PAF.Plot)
                                      }
                                    }
                                    if(input$PIForPAF > 1){
                                      if(input$Cft == 1){
                                        Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                                        PIF.Plot <- plotpif(X=X, thetamin = thetamin, thetamax = thetamax, rr = RR, cft = Cft, weights =  Weight)
                                        if(input$PIForPAF == 2){
                                          return(PIF.Plot)
                                        }else{
                                          return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
                                        }
                                      }else{
                                        if(is.null(input$Cftfun)){
                                          Error <- 4; return(Warnings2(Error));
                                        }else{
                                          if(file_ext(input$Cftfun$name)!="R"){
                                            Error <- 5; return(Warnings2(Error));
                                          }else{
                                            source(input$Cftfun$datapath)
                                            if(exists("Counterfactual") == FALSE){
                                              Error <- 6; return(Warnings2(Error));
                                            }else{
                                              if(is.function(Counterfactual) == FALSE){
                                                Error <- 7; return(Warnings2(Error));
                                              }else{
                                                if(is.numeric(try(Counterfactual(0),silent = TRUE)) == FALSE || is.na(Counterfactual(0))){
                                                  Error <- 8; return(Warnings2(Error));
                                                }else{
                                                  Cft <- Counterfactual
                                                  PIF.Plot <- plotpif(X=X, thetamin = thetamin, thetamax = thetamax, rr = RR, cft = Cft, weights = Weight)
                                                  if(input$PIForPAF == 2){
                                                    return(PIF.Plot)
                                                  }else{
                                                    return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
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
      }else{#si xmin==0
        xmin <- input$XMin
        if(xmin<0){
          Error <- 19; return(Warnings2(Error));
        }else{
          if(rr != 5){
            Error <- 21; return(Warnings2(Error));
          }else{
            if(is.null(input$X)){
              Error <- 1; return(Warnings2(Error));
            }else{
              if(file_ext(input$X$name)!="csv"){
                Error <- 2; return(Warnings2(Error));
              }else{
                X   <- as.numeric(as.matrix(read.csv(input$X$datapath,stringsAsFactors = FALSE)))
                
                if(is.null(input$RRfun)){
                  Error <- 11; return(Warnings2(Error));
                }else{
                  if(file_ext(input$RRfun$name)!="R"){
                    Error <- 12; return(Warnings2(Error));
                  }else{
                    source(input$RRfun$datapath)
                    if(exists("RRfunction") == FALSE){
                      Error <- 13; return(Warnings2(Error));
                    }else{
                      if(is.function(RRfunction) == FALSE){
                        Error <- 14; return(Warnings2(Error));
                      }else{
                        
                        if(is.numeric(try(RRfunction(xmin,thetamin))) == FALSE || is.na(RRfunction(xmin,thetamin))){
                          Error <- 17; return(Warnings2(Error));
                        }else{
                          if(RRfunction(xmin,thetamin) != 1 && RRfunction(xmin,thetamax) != 1){
                            Error <- 22; return(Warnings2(Error));
                          }else{
                            RR <- RRfunction
                            if (input$w==0) {
                              
                              PAF.Plot <- plotpif(X=X, thetamin = thetamin, thetamax = thetamax, rr = RR, cft = function(X){sapply(X,function(x){xmin})}) + ggtitle("Population Attributable Fraction (PAF) \n under different values of theta")
                              if(input$PIForPAF == 1){
                                if(input$PIForPAF == 1){
                                  return(PAF.Plot)
                                }
                              }
                              if(input$PIForPAF > 1){
                                if(input$Cft == 1){
                                  Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                                  PIF.Plot <- plotpif(X=X, thetamin = thetamin, thetamax = thetamax, rr = RR, cft = Cft)
                                  if(input$PIForPAF == 2){
                                    return(PIF.Plot)
                                  }else{
                                    return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
                                  }
                                }else{
                                  if(is.null(input$Cftfun)){
                                    Error <- 4; return(Warnings2(Error));
                                  }else{
                                    if(file_ext(input$Cftfun$name)!="R"){
                                      Error <- 5; return(Warnings2(Error));
                                    }else{
                                      source(input$Cftfun$datapath)
                                      if(exists("Counterfactual") == FALSE){
                                        Error <- 6; return(Warnings2(Error));
                                      }else{
                                        if(is.function(Counterfactual) == FALSE){
                                          Error <- 7; return(Warnings2(Error));
                                        }else{
                                          if(is.numeric(try(Counterfactual(xmin),silent = TRUE)) == FALSE || is.na(Counterfactual(xmin))){
                                            Error <- 8; return(Warnings2(Error));
                                          }else{
                                            Cft <- Counterfactual
                                            PIF.Plot <- plotpif(X=X, thetamin = thetamin, thetamax = thetamax, rr = RR, cft = Cft)
                                            if(input$PIForPAF == 2){
                                              return(PIF.Plot)
                                            }else{
                                              return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
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
                              #
                              if(is.null(input$weight)){
                                Error <- 9; return(Warnings2(Error));
                              }else{
                                if(file_ext(input$weight$name)!="csv"){
                                  Error <- 10; return(Warnings2(Error));
                                }else{
                                  Weight <- as.numeric(as.matrix(read.csv(input$weight$datapath,stringsAsFactors = FALSE)))
                                  
                                  if(length(Weight)!=length(X)){
                                    return(ggplot()+ggtitle(sprintf("Warning: The amount of exposure values and weight values \n are different, there are \n  %d exposure values, %d weight values." , length(X), length(Weight))))
                                  }else{
                                    PAF.Plot <- plotpif(X=X, thetamin = thetamin, thetamax = thetamax, rr = RR, cft = function(X){sapply(X,function(x){xmin})}, weights = Weight) + ggtitle("Population Attributable Fraction (PAF) \n under different values of theta")
                                    if(input$PIForPAF==1){
                                      if(input$PIForPAF == 1){
                                        return(PAF.Plot)
                                      }
                                    }
                                    if(input$PIForPAF > 1){
                                      if(input$Cft == 1){
                                        Cft <- function(X){as.numeric(input$a)*X+as.numeric(input$b)}
                                        PIF.Plot <- plotpif(X=X, thetamin = thetamin, thetamax = thetamax, rr = RR, cft = Cft, weights =  Weight)
                                        if(input$PIForPAF == 2){
                                          return(PIF.Plot)
                                        }else{
                                          return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
                                        }
                                      }else{
                                        if(is.null(input$Cftfun)){
                                          Error <- 4; return(Warnings2(Error));
                                        }else{
                                          if(file_ext(input$Cftfun$name)!="R"){
                                            Error <- 5; return(Warnings2(Error));
                                          }else{
                                            source(input$Cftfun$datapath)
                                            if(exists("Counterfactual") == FALSE){
                                              Error <- 6; return(Warnings2(Error));
                                            }else{
                                              if(is.function(Counterfactual) == FALSE){
                                                Error <- 7; return(Warnings2(Error));
                                              }else{
                                                if(is.numeric(try(Counterfactual(xmin),silent = TRUE)) == FALSE || is.na(Counterfactual(xmin))){
                                                  Error <- 8; return(Warnings2(Error));
                                                }else{
                                                  Cft <- Counterfactual
                                                  PIF.Plot <- plotpif(X=X, thetamin = thetamin, thetamax = thetamax, rr = RR, cft = Cft, weights = Weight)
                                                  if(input$PIForPAF == 2){
                                                    return(PIF.Plot)
                                                  }else{
                                                    return(grid.arrange(PAF.Plot,PIF.Plot,ncol = 1))
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
      }
    }
  })
  
  #############################################
  #Confidence intervals
  
  output$IntervalPAF <- renderText({
    conf <- as.numeric(input$Conf)
    rr   <- as.numeric(input$rr)
    if(input$PIForPAF>1 || input$xmin0 != 1){
      return("Confidence intervals only for PAF with minimum exposure level equal to zero.")
    }else{
      if(conf <= 0 || conf >= 100){
        return("Warning: Confidence level must be greater than zero and less than a hundred.")
      }else{
        if(rr <= 2 && input$xmin0 == 1){
          SD <- as.matrix(as.numeric(input$SDtheta))
          if(SD < 0){
            return("Standard deviation must be greater or equal to zero.")
          }else{
            if(is.null(input$X)){
              Error <- 1; return(Warnings1(Error));
            }else{
              if(file_ext(input$X$name)!="csv"){
                Error <- 2; return(Warnings1(Error));
              }else{
                X     <- as.numeric(as.matrix(read.csv(input$X$datapath,stringsAsFactors = FALSE)))
                theta <- switch(rr, c(input$theta),
                                c( input$theta1))
                theta <- as.numeric(theta)
                RR    <- switch(rr, function(X,theta){exp(theta[1]*X)},
                                function(X,theta){theta[1]*X+1})
                if(theta<0){
                  Error <- 3; return(Warnings1(Error));
                }else{
                  if (input$w==0) {
                    
                    PAF.interval <- paf.confidence.linear(X=X, thetahat = theta, thetasd = SD, rr = RR, confidence = conf)
                    
                    return(sprintf("Confidence interval PAF :  (%f, %f) <br/> Point Estimate :  %f <br/> Estimated Variance :  %f", PAF.interval[1], PAF.interval[3], PAF.interval[2], PAF.interval[4]))
                    
                  }
                }
              }
            }
          }
        }else{
          if(is.null(input$SDMat)){
            Error <- 21; return(Warnings1(Error));
          }else{
            if(file_ext(input$SDMat$name)!="csv"){
              Error <- 22; return(Warnings1(Error));
            }else{
              SD  <- as.matrix(read.csv(input$SDMat$datapath,stringsAsFactors = FALSE))
              if(is.matrix(SD) == FALSE){
                Error <- 23; return(Warnings1(Error));
              }else{
                if(min(SD)<0){
                  Error <- 24; return(Warnings1(Error));
                }else{
                  if(is.positive.definite(SD) == FALSE){
                    Error <- 25; return(Warnings1(Error));
                  }else{
                    if(is.null(input$X)){
                      Error <- 1; return(Warnings1(Error));
                    }else{
                      if(file_ext(input$X$name)!="csv"){
                        Error <- 2; return(Warnings1(Error));
                      }else{
                        X   <- as.numeric(as.matrix(read.csv(input$X$datapath,stringsAsFactors = FALSE)))
                        
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
                            Error <- 3; return(Warnings1(Error));
                          }else{
                            if (input$w==0) {
                              
                              PAF.interval <- paf.confidence.linear(X=X, thetahat = theta, thetasd = SD, rr = RR, confidence = conf)
                              return(sprintf("Confidence interval PAF :  (%f, %f) <br/> Point Estimate :  %f <br/> Estimated Variance :  %f", PAF.interval[1], PAF.interval[3], PAF.interval[2], PAF.interval[4]))
                            }else{
                              #Para x xmin 0 con pesos distintos
                              if(is.null(input$weight)){
                               Error <- 9; return(Warnings1(Error));
                              }else{
                                if(file_ext(input$weight$name)!="csv"){
                                  Error <- 10; return(Warnings1(Error));
                                }else{
                                  Weight <- as.numeric(as.matrix(read.csv(input$weight$datapath,stringsAsFactors = FALSE)))
                                  
                                  if(length(Weight)!=length(X)){
                                    return(sprintf("Warning: The amount of exposure values and weight values are different, there are %d exposure values, %d weight values." , length(X), length(Weight)))
                                  }else{
                                    PAF.interval <- paf.confidence.linear(X=X, thetahat = theta, thetasd = SD, rr = RR, weights = Weight)
                                    return(sprintf("Confidence interval PAF :  (%f, %f) <br/> Point Estimate :  %f <br/> Estimated Variance :  %f", PAF.interval[1], PAF.interval[3], PAF.interval[2], PAF.interval[4]))
                                  }
                                }
                              }
                            }
                          }
                        }else{#para xmin 0 con rr=5
                          if(is.null(input$RRfun)){
                            Error <- 11; return(Warnings1(Error));
                          }else{
                            if(file_ext(input$RRfun$name)!="R"){
                              Error <- 12; return(Warnings1(Error));         
                              }else{
                              source(input$RRfun$datapath)
                              if(exists("RRfunction") == FALSE){
                                Error <- 13; return(Warnings1(Error));     
                                }else{
                                if(is.function(RRfunction) == FALSE){
                                  Error <- 14; return(Warnings1(Error));      
                                  }else{
                                  if(is.null(input$thetahat)){
                                    Error <- 15; return(Warnings1(Error));      
                                    }else{
                                    if(file_ext(input$thetahat$name)!="csv"){
                                      Error <- 16; return(Warnings1(Error));   
                                      }else{
                                      theta   <- as.numeric(as.matrix(read.csv(input$thetahat$datapath,stringsAsFactors = FALSE)))
                                      if(is.numeric(try(RRfunction(0,theta))) == FALSE || is.na(RRfunction(0,theta))){
                                        Error <- 17; return(Warnings1(Error));
                                        }else{
                                        if(RRfunction(0,theta) != 1){
                                          Error <- 18; return(Warnings1(Error));
                                          }else{
                                          RR <- RRfunction
                                          if (input$w==0) {
                                            
                                            PAF.interval <- paf.confidence.linear(X=X, thetahat = theta, thetasd = SD, rr = RR, confidence = conf)
                                            return(sprintf("Confidence interval PAF :  (%f, %f) <br/> Point Estimate :  %f <br/> Estimated Variance :  %f", PAF.interval[1], PAF.interval[3], PAF.interval[2], PAF.interval[4]))
                                          } else{
                                            #Para x xmin 0 con pesos distintos
                                            if(is.null(input$weight)){
                                             Error <- 9; return(Warnings1(Error));
                                            }else{
                                              if(file_ext(input$weight$name)!="csv"){
                                                Error <- 10; return(Warnings1(Error));
                                              }else{
                                                Weight <- as.numeric(as.matrix(read.csv(input$weight$datapath,stringsAsFactors = FALSE)))
                                                
                                                if(length(Weight)!=length(X)){
                                                  Error <- 11; return(Warnings1(Error));
                                                  }else{
                                                  PAF.interval <- paf.confidence.linear(X=X, thetahat = theta, thetasd = SD, rr = RR, weights = Weight)
                                                  
                                                  return(sprintf("Confidence interval PAF :  (%f, %f) <br/> Point Estimate :  %f <br/> Estimated Variance :  %f", PAF.interval[1], PAF.interval[3], PAF.interval[2], PAF.interval[4]))
                                                  
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
    }
  })
  
})
