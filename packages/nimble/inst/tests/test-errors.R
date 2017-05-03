source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of error handling.")

test_that("Testing of error handling in SEXP_2_NimArr", {
  # First construct a nimbleFunction.
  nimFun <- nimbleFunction(
    run = function (x = logical(1),
                    y = integer(1),
                    z = double(1)) {
      return(x + y + z)
      returnType(double(1))
    }
  )
  
  # Make sure the nimble function works on well formed inputs.
  x <- c(FALSE, TRUE)
  y <- c(2, 4)
  z <- c(8.0, 16.0)
  expectedResult <- c(10.0, 21.0)
  expect_equal(nimFun(x, y, z), expectedResult)
  
  # Next compile the nimbleFunction.
  compiledFun = compileNimble(nimFun)
  expect_equal(compiledFun(x, y, z), expectedResult)
  
  # Finally check that compiledFun behaves as expected on invalid inputs.
  wrong_type = c("x", "y")
  wrong_shape1 = TRUE
  wrong_shape2 = matrix(c(1, 2, 3, 4), 2, 2)
  
  # These DO NOT trigger errors.
  expect_failure(expect_error(compiledFun(wrong_shape1, y, z)))
  expect_failure(expect_error(compiledFun(x, wrong_shape1, z)))
  expect_failure(expect_error(compiledFun(x, y, wrong_shape1)))
  expect_failure(expect_error(compiledFun(x, wrong_shape2, z)))
  expect_failure(expect_error(compiledFun(x, y, wrong_shape2)))
  
  # These DO trigger errors.
  expect_error(compiledFun(wrong_type, y, z))
  expect_error(compiledFun(x, wrong_type, z))
  expect_error(compiledFun(x, y, wrong_type))
  expect_error(compiledFun(wrong_shape2, y, z))
})
