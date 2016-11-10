context("General PAF confidence intervals")

test_that("Checking paf.confidence",{
  
  #Check results are the same for corresponding methods when variance of thetas is equal to zero
  #Inverse Empirical
  expect_equal({
    set.seed(76787)
    X        <- rnorm(100,3,1)
    thetahat <- 0.2
    thetavar <- 0
    rr       <- function(X,theta){exp(X*theta)}
    paf.confidence(X, thetahat, thetavar, rr, method = "inverse", est.method = "empirical")
  },
  {
    paf.confidence.inverse(X, thetahat, thetavar, rr, method = "empirical")
  })
  
  #Inverse Approximate
  expect_equal({
    set.seed(76787)
    X        <- mean(rnorm(100,3,1))
    Xvar     <- 0.2
    thetahat <- 0.2
    thetavar <- 0
    rr       <- function(X,theta){exp(X*theta)}
    paf.confidence(X, thetahat, thetavar, rr, method = "inverse", est.method = "approximate", Xvar = Xvar)
  },
  {
    paf.confidence.inverse(X, thetahat, thetavar, rr, method = "approximate", Xvar = Xvar)
  })
  
  #One to one Empirical
  expect_equal({
    set.seed(76787)
    X        <- rnorm(100,3,1)
    thetahat <- 0.2
    thetalow <- 0.1
    thetaup  <- 0.3
    rr       <- function(X,theta){exp(X*theta)}
    paf.confidence(X, thetahat, thetamin = thetalow, thetamax = thetaup, rr = rr, method = "one2one", est.method = "empirical")
  },
  {
    paf.confidence.one2one(X, thetahat, thetalow, thetaup, rr, method = "empirical")
  })
  
  #One to one Approximate
  expect_equal({
    set.seed(76787)
    X        <- mean(rnorm(100,3,1))
    Xvar     <- 1
    thetahat <- 0.2
    thetalow <- 0.1
    thetaup  <- 0.3
    rr       <- function(X,theta){exp(X*theta)}
    paf.confidence(X, thetahat, thetamin = thetalow, thetamax = thetaup, rr = rr, method = "one2one", est.method = "approximate", Xvar = Xvar)
  },
  {
    paf.confidence.one2one(X, thetahat, thetalow, thetaup, rr, method = "approximate", Xvar = Xvar)
  })
  
  #Linear Empirical
  expect_equal({
    set.seed(76787)
    X        <- rnorm(100,3,1)
    thetahat <- 0.2
    thetavar <- 0
    rr       <- function(X,theta){exp(X*theta)}
    paf.confidence(X, thetahat, thetavar, rr, method = "linear", est.method = "empirical")
  },
  {
    paf.confidence.linear(X, thetahat, thetavar, rr)
  })
  
  #Linear Approximate
  
  expect_equal({
    set.seed(76787)
    X        <- mean(rnorm(100,3,1))
    Xvar     <- 1
    thetahat <- 0.2
    thetavar <- 0
    rr       <- function(X,theta){exp(X*theta)}
    paf.confidence(X, thetahat, thetavar = thetavar, rr = rr, method = "linear", est.method = "approximate", Xvar = Xvar)
  },
  {
    paf.confidence.approximate(X, Xvar, thetahat, thetavar, rr)
  })
  
  #Loglinear Empirical
  expect_equal({
    set.seed(76787)
    X        <- rnorm(100,3,1)
    thetahat <- 0.2
    thetavar <- 0
    rr       <- function(X,theta){exp(X*theta)}
    paf.confidence(X, thetahat, thetavar, rr, method = "log", est.method = "empirical")
  },
  {
    paf.confidence.loglinear(X, thetahat, thetavar, rr)
  })
  
  #Loglinear Approximate
  
  expect_equal({
    set.seed(76787)
    X        <- mean(rnorm(100,3,1))
    Xvar     <- 1
    thetahat <- 0.2
    thetavar <- 0
    rr       <- function(X,theta){exp(X*theta)}
    paf.confidence(X, thetahat, thetavar = thetavar, rr = rr, method = "log", est.method = "approximate", Xvar = Xvar)
  },
  {
    paf.confidence.approximate.loglinear(X, Xvar, thetahat, thetavar, rr)
  })
  
})