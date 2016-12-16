context("Confidence intervals for paf approximate loglinear")

test_that("Checking paf_confidence_approximate_loglinear",{
  
  #Check that if variance of X is zero and X is zero paf is zero and so the upper and lower bounds of its CI
  expect_equal({
    Xmean    <- 0
    Xvar     <- 0
    thetahat <- 0.2
    thetavar <- 0.003
    rr       <- function(X,theta){exp(X*theta)}
    paf.confidence.approximate.loglinear(Xmean, Xvar, thetahat, thetavar, rr)
  },
  c("Lower" = 0, "Point_Estimate" =  0, "Upper" = 0) 
  )
  
  #Check that if there's no variance in theta the result is the same
  expect_equal({
    Xmean   <- 1
    Xvar    <- .3
    rr      <- function(X,theta){exp(X*theta)}
    theta   <- 0.2
    thetasd <- 0
    paf.confidence.approximate.loglinear(Xmean, Xvar, theta, thetasd, rr)
  },
  {Xmean   <- 1
  Xvar    <- .3
  rr      <- function(X,theta){exp(X*theta)}
  theta   <- 0.2
  thetasd <- 0
  paf.confidence.approximate.loglinear(Xmean, Xvar, theta, thetasd, rr)
  })
  
  #Check errors 
  #Covariance matrix of theta not positive semi-definite
  expect_error({
    Xmean   <- 1
    Xvar    <- .3
    rr      <- function(X,theta){exp(X*theta)}
    theta   <- 0.2
    thetasd <- - 0.02
    paf.confidence.approximate.loglinear(Xmean, Xvar, theta, thetasd, rr)
  })
  
  #Exposure mean less than zero
  expect_warning({
    Xmean   <- -1
    Xvar    <- .3
    rr      <- function(X,theta){exp(X*theta)}
    theta   <- 0.2
    thetasd <- 0.02
    paf.confidence.approximate.loglinear(Xmean, Xvar, theta, thetasd, rr)
  })
  
  #Covariance matrix of exposure values not positive semi-definite
  expect_error({
    Xmean   <- 1
    Xvar    <- - .3
    rr      <- function(X,theta){exp(X*theta)}
    theta   <- 0.2
    thetasd <- 0.02
    paf.confidence.approximate.loglinear(Xmean, Xvar, theta, thetasd, rr)
  })
  
  #Exposure diminishes risk
  expect_warning({
    Xmean   <- 1
    Xvar    <- .3
    rr      <- function(X,theta){exp(X*theta)}
    theta   <- -0.2
    thetasd <- 0.02
    paf.confidence.approximate.loglinear(Xmean, Xvar, theta, thetasd, rr)
  })
  
  
})
