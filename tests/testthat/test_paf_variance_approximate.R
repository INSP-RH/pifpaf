context("Check PAF variance for approximate method")

test_that("Checking paf_variance_approximate",{
  
  expect_equal({
    Xmean    <- 2
    Xvar     <- 0.5
    thetahat <- 0.02
    thetasd  <- 0
    rr       <- function(X,theta){exp(X*theta)}
    paf.variance.approximate(Xmean, Xvar, thetahat, thetasd, rr)
  },
  {(thetahat/(exp(Xmean*thetahat)))^2*Xvar})
  
  expect_equal({
    Xmean    <- 2
    Xvar     <- 0.5
    thetahat <- 0.02
    thetasd  <- 0
    rr       <- function(X,theta){X*theta+1}
    paf.variance.approximate(Xmean, Xvar, thetahat, thetasd, rr)
  },
  {(thetahat/((thetahat*Xmean+1)^2))^2*Xvar})
  
})
