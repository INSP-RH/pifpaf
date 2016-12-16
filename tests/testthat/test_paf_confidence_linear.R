context("CI of PAF using linear method")

test_that("Checking paf.confidence.linear",{
  
  #Expect error  when confidence levels are incorrectly specified
  expect_error({
    set.seed(347618)
    X        <- rnorm(100,3,1)
    thetahat <- 1.4
    thetavar <- 0.05
    rr       <- function(X,theta){exp(X*theta)}
    paf.confidence.linear(X, thetahat, thetavar, rr, confidence = -1)
  })
  
  #Expect error  when  variance is not positive semi-definite
  
  expect_error({
    set.seed(347618)
    X        <- rnorm(100,3,1)
    thetahat <- 1.4
    thetavar <- -0.05
    rr       <- function(X,theta){exp(X*theta)}
    paf.confidence.linear(X, thetahat, thetavar, rr, confidence = 95)
    })
  
  expect_error({
    set.seed(347618)
    X        <- rnorm(100,3,1)
    thetahat <- c(1.4, 1)
    thetavar <- matrix(c(.2,.9,.9, .1), ncol = 2, byrow = TRUE)
    rr       <- function(X,theta){X*theta[1]+theta[2]}
    paf.confidence.linear(X, thetahat, thetavar, rr, confidence = 95)
  })
  
  #Check dimensions are correctly specified
  expect_error({
    set.seed(347618)
    X        <- rnorm(100,3,1)
    thetahat <- .02
    thetavar <- matrix(c(.2,.01,.01, .1), ncol = 2, byrow = TRUE)
    rr       <- function(X,theta){X*theta[1]+theta[2]}
    paf.confidence.linear(X, thetahat, thetavar, rr, confidence = 95)
  })
  
  #All exposure levels are zero expect paf equal to zero
  expect_equal({
    set.seed(347618)
    X        <- rep(0,100)
    thetahat <- 0.2
    thetavar <- 0.02
    rr       <- function(X,theta){exp(X*theta)}
    paf.confidence.linear(X, thetahat, thetavar, rr, confidence = 95)
  },
  
  c("Lower_CI" = 0, "Point_Estimate" = 0, "Upper_CI" = 0, "Estimated_Variance" = 0))
  
  
})
