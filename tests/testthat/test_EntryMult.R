context("Entry by entry matrix multiplications")

test_that("Checking EntryMult",{
  
  #Check it gives sum (AijBij)
  expect_equal(
    EntryMult(2, 7),
    14
  )
  
  expect_equal(
    EntryMult(c(1,2,3), c(3,2,1)),
    10
  )
  
  expect_equal(
    EntryMult(matrix(c(1,2,3,4), ncol = 2), matrix(rep(2, 4), ncol = 2)),
    20
  )
  
  expect_equal(
    EntryMult(matrix(c(1,2,3,4, 1,2,3,4, 1,1,1,1, 1,1,1,1), ncol = 4), matrix(rep(1, 16), ncol = 4)),
    28
  )
  
  #Check errors when dimensions are different
  
  expect_error({
    X <- 3
    Y <- c(2,7)
    EntryMult(X,Y)
  })
  
  
})