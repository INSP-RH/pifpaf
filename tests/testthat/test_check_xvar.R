context("Check Xvar")

test_that("Checking check.xvar",{
  
  #Check that for positive values of meancft and meanrr gives TRUE
  expect_warning({
    check.xvar(NA)
  })
  
  expect_equal(
    {n <- runif(1)
    check.xvar(n)},
    n
  )
  
})