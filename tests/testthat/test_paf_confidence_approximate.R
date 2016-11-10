context("Approximate Population Attributable Fraction Confidence Intervals")

test_that("Checking paf.confidence.approximate",{
  
  #Check that it gives correct results when exposure values are all zero
  expect_equal({
    
    Xmean    <- 0
    Xvar     <- 0
    thetahat <- 1.4
    thetavar <- 0.02
    rr       <- function(X,theta){exp(X*theta)}
    paf.confidence.approximate(Xmean, Xvar, thetahat, thetavar, rr)
    
  },
  c("Lower_CI" = 0, "Point_Estimate" = 0, "Upper_CI" = 0, "Estimated_Variance" = 0))
  
  #Check that it throws an error when exposure values are negative
  expect_warning({
    Xmean    <- -1
    Xvar     <- 0.2
    thetahat <- 1.4
    thetavar <- 0.02
    rr       <- function(X,theta){exp(X*theta)}
    paf.confidence.approximate(Xmean, Xvar, thetahat, thetavar, rr)
  })
  
  #Check that there's an error when variance is not positive definite
  expect_error({
    Xmean    <- 1
    Xvar     <- -0.2
    thetahat <- 1.4
    thetavar <- 0.02
    rr       <- function(X,theta){exp(X*theta)}
    paf.confidence.approximate(Xmean, Xvar, thetahat, thetavar, rr)
  })
  
  #Check that there's an error when variance is not positive definite
  expect_error({
    Xmean    <- 1
    Xvar     <- 0.2
    thetahat <- 1.4
    thetavar <- -0.02
    rr       <- function(X,theta){exp(X*theta)}
    paf.confidence.approximate(Xmean, Xvar, thetahat, thetavar, rr)
  })
  
  #Check there's a warning when risk decreases as exposure increases
  expect_warning({
    Xmean    <- 1
    Xvar     <- 0.2
    thetahat <- -1.4
    thetavar <- 0.02
    rr       <- function(X,theta){exp(X*theta)}
    paf.confidence.approximate(Xmean, Xvar, thetahat, thetavar, rr)
  })
  
  expect_error({
    Xmean    <- c(1,2)
    Xvar     <- c(0.2,0.1)
    thetahat <- 1.4
    thetavar <- 0.02
    rr       <- function(X,theta){exp(theta*(X[1]+X[2]))}
    paf.confidence.approximate(Xmean, Xvar, thetahat, thetavar, rr)
  })
})