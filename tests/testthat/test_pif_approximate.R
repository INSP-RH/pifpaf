context("Approximate Population Impact Fraction point-estimate")



test_that("Checking pif.approximate function errors",{
  
  
  #Check that mean exposure values are greater than zero
  expect_warning({
    X     <- runif(1000, -1,0) 
    Xmean <- mean(X)
    Xvar  <- var(X)
    thetahat <- 1.4
    pif.approximate(Xmean, Xvar, thetahat, function(X, theta){theta*X + 1})
    
  })
  
  #Check that relative risk > 0
  expect_error({
    X     <- runif(1000,0,2)
    Xmean <- mean(X)
    Xvar  <- var(X)
    thetahat <- 1.4
    pif.approximate(Xmean, Xvar, thetahat, function(X, theta){-theta*X + 1})
    
  })
  

  
  #Check that counterfactual relative risk > 0
  expect_error({
    
    X     <- rnorm(100, 2, .5)
    Xmean <- mean(X)
    Xvar  <- var(X)
    thetahat <- 1.4
    pif.approximate(Xmean, Xvar, thetahat, rr = function(X, theta){X*theta + 1}, 
                  cft = function(X){-100*X})
    
  })
  
  #Check that variance and covariance matrix is positive semi-definite
  expect_error({
    
    X     <- runif(1000,0,3)
    Xmean <- mean(X)
    Xvar  <- -var(X)
    thetahat <- 1.4
    pif.approximate(Xmean, Xvar, thetahat, function(X, theta){X*theta + 1}, cft = function(X){0.5*X})
    
  })
  
  #Check that variance and covariance matrix dimensions go according to the dimension of X
  expect_error({
    
    X     <- runif(1000,0,3)
    Xmean <- mean(X)
    Xvar  <- cov(matrix(X, ncol = 2))
    thetahat <- 1.4
    pif.approximate(Xmean, Xvar, thetahat, function(X, theta){X*theta + 1}, cft = function(X){0.5*X})
    
  })
  
})

test_that("Checking pif.approximate function warnings",{
  
  #Check that RR(0, theta) = 1
  expect_warning({
    
    X        <- cbind(rnorm(100, 3,.8), rbeta(100, 0.5, 0.5))
    Xmean    <- matrix(colMeans(X), ncol = 2)
    Xvar     <- matrix(cov(X),  ncol = 2)
    thetahat <- 1.4
    pif.approximate(Xmean, Xvar, thetahat, function(X, theta){X[,1]*X[,2]*theta + 2})
    
  })
  
  #Check that RR under counterfactual is smaller than RR under normal circumstances
  expect_warning({
    
    X     <- 1:20
    Xmean <- mean(X)
    Xvar  <- var(X)
    thetahat <- 1.4
    pif.approximate(Xmean, Xvar, thetahat, function(X, theta){X*theta + 1}, cft = function(X){100*X})
    
  })
  
  
  
})

test_that("Checking pif.approximate convergence",{
  
  #Check that approximate PAF works
  expect_equal(
    pif.approximate(Xmean = 3, Xvar = 1, thetahat = 1, rr = function(X, theta){exp(theta*X)}),
    1 - 1/(1.5*exp(3))
  )
  
  #Check that approximate PIF works when counterfactual is identity
  expect_equal(
    pif.approximate(Xmean = 3, Xvar = 1, thetahat = 1, rr = function(X, theta){exp(theta*X)}, cft = function(X){X}),
    0
  )
  
  #Check that approximate PAF works when RR is constant 1
 
  expect_equal(
    pif.approximate(Xmean = 3, Xvar = 1, thetahat = 1, rr = function(X, theta){1}),
    0
  )
  
  #Check that approximate PIF works when RR is constant 1
  expect_equal(
    pif.approximate(Xmean = 3, Xvar = 1, thetahat = 1, rr = function(X, theta){1}, cft = function(X){(0.4*X)^2}),
    0
  )
  
  #Check that approximate PIF works when counterfactual is constant 0 is equal to empiric
  expect_equal(
    pif.approximate(Xmean = 3, Xvar = 1, thetahat = 1, rr = function(X, theta){1}, cft = function(X){0}),
    pif.approximate(Xmean = 3, Xvar = 1, thetahat = 1, rr = function(X, theta){1})
  )
  
})
