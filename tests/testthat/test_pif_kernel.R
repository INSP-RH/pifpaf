context("Kernel-based Population Impact Fraction point-estimate")

test_that("Checking pif.kernel function errors",{
  
  #Check that relative risk > 0
  expect_error({
    X <- rnorm(100, 4,1)
    thetahat <- 1.4
    pif.kernel(X, thetahat, function(X, theta){-theta*X + 1})
    
  })
  
  #Check that exposure levels are positive
  expect_warning({
    X <- rnorm(100, -1)
    thetahat <- 1.4
    pif.kernel(X, thetahat, function(X, theta){exp(theta*X)})
    
  })
  
  #Check that counterfactual relative risk > 0
  expect_error({
    X        <- rnorm(100, 4,1)
    thetahat <- 1.4
    pif.kernel(X, thetahat, rr = function(X, theta){X*theta + 1}, 
                  cft = function(X){-100*X})
    
  })
  
  #Check that X has only one column
  expect_error({
    
    X <- cbind(rnorm(100, 2, 0.5), rnorm(100, 2, 0.5))
    thetahat <- 1.4
    pif.kernel(X, thetahat, rr = function(X, theta){X[,1]*theta + 1}, 
               cft = function(X){-100*X})
    
  })
})

test_that("Checking pif.kernel function warnings",{
  
  #Check that RR(0, theta) = 1
  expect_warning({
    
    X <- rnorm(100,4,1)
    thetahat <- 1.4
    pif.kernel(X, thetahat, function(X, theta){X*theta + 2})
    
  })
  
  #Check that RR under counterfactual is smaller than RR under normal circumstances
  expect_warning({
    X <- c(1,2,3)
    thetahat <- 1.4
    pif.kernel(X, thetahat, function(X, theta){X*theta + 1}, cft = function(X){100*X})
    
  })
  
})

test_that("Checiking pif.kernel convergence",{
  
  #Check that empirical PAF works
  X <- rnorm(10000, 1, sd = 0.001)
  expect_equal(
    pif.kernel(X, 0.1, rr = function(X, theta){exp(theta*X)}),
    pif.empirical(X, 0.1, rr = function(X, theta){exp(theta*X)})
  )
  
  #Check that empirical PIF works when counterfactual is identity
  expect_equal(
    pif.kernel(c(1,2,3), 1, rr = function(X, theta){exp(theta*X)}, cft = function(X){X}),
    0
  )
  
  #Check that empirical PAF works when RR is constant 1
  expect_equal(
    pif.kernel(c(1,2,3), 1, rr = function(X, theta){1}),
    0
  )
  
  #Check that empirical PIF works when RR is constant 1
  expect_equal(
    pif.kernel(c(1,2,3), 1, rr = function(X, theta){1}, cft = function(X){X^2}),
    0
  )
  
})

