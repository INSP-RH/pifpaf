context("Check_exposure")

test_that("Checking check exposure",{
  
  #Check that if X has a negative (exposure) value the function throws an error
  expect_warning({
    X <- runif(100,-.2, 1)
    check.exposure(X)
  })
  
  expect_warning({
    X <- -.1
    check.exposure(X)
  })
  
  #Check that if all values are positive returns true
  expect_null({
    X <- runif(100)
    check.exposure(X)
  })
  
  expect_null({
    X <- rlnorm(100)
    check.exposure(X)
  })
})