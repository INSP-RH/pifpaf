context("Empirical Population Impact Fraction point-estimate")

test_that("Checking pif.empirical function errors",{
  
  #Check that relative risk > 0
  expect_error({
    
    X <- rnorm(100, 1)
    thetahat <- 1.4
    pif.empirical(X, thetahat, function(X, theta){-theta*X + 1})
    
  })
  
  #Check that relative risk > 0 (second check)
  expect_error({
    
    X <- rnorm(100, -1)
    thetahat <- 1.4
    pif.empirical(X, thetahat, function(X, theta){theta*X + 1})
    
  })
  
  #Check that counterfactual relative risk > 0
  expect_error({
    
    X <- rnorm(100, 1)
    thetahat <- 1.4
    pif.empirical(X, thetahat, rr = function(X, theta){X*theta + 1}, 
        cft = function(X){-100*X})
    
  })
})

test_that("Checking pif.empirical function warnings",{
  
  #Check that RR(0, theta) = 1
  expect_warning({
    
    X <- cbind(rnorm(100), rbeta(100, 0.5, 0.5))
    thetahat <- 1.4
    pif.empirical(X, thetahat, function(X, theta){X[,1]*X[,2]*theta + 2})
    
  })
  
})

test_that("Checiking pif.empirical convergence",{
  
    #Check that empirical PAF works
    expect_equal(
      pif.empirical(c(1,2,3), 1, rr = function(X, theta){exp(theta*X)}),
      1 - 3/sum(exp(c(1,2,3)))
    )
    
    #Check that empirical PIF works when counterfactual is identity
    expect_equal(
      pif.empirical(c(1,2,3), 1, rr = function(X, theta){exp(theta*X)}, cft = function(X){X}),
      0
    )
    
    #Check that empirical PAF works when RR is constant 1
    expect_equal(
      pif.empirical(c(1,2,3), 1, rr = function(X, theta){1}),
      0
    )
    
    #Check that empirical PIF works when RR is constant 1
    expect_equal(
      pif.empirical(c(1,2,3), 1, rr = function(X, theta){1}, cft = function(X){X^2}),
      0
    )
    
    #Check that empirical PIF works when counterfactual is constant 0
    expect_equal(
      pif.empirical(c(1,2,3), 1, rr = function(X, theta){exp(theta*X)}, cft = function(X){0}),
      1 - 3/sum(exp(c(1,2,3)))
    )
  
})
