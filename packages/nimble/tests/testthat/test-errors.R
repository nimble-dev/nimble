source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

test_that("Testing of stopping on run-time size errors", {
  # First construct a nimbleFunction.
  current_option <- getNimbleOption("stopOnSizeErrors")
  nimbleOptions(stopOnSizeErrors = TRUE)
  on.exit( nimbleOptions(stopOnSizeErrors = current_option ) )
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
  
  # January 2024, branch stop-on-size-errors: We are making these nimStop by default
  expect_error(compiledFun(wrong_shape1, y, z))
  expect_error(compiledFun(x, wrong_shape1, z))
  expect_error(compiledFun(x, y, wrong_shape1))
  expect_error(compiledFun(x, wrong_shape2, z))
  expect_error(compiledFun(x, y, wrong_shape2))
  
  # These DO trigger errors.
  expect_error(compiledFun(wrong_type, y, z))
  expect_error(compiledFun(x, wrong_type, z))
  expect_error(compiledFun(x, y, wrong_type))
  expect_error(compiledFun(wrong_shape2, y, z))
})

test_that("Testing of warning on run-time size errors", {
  # First construct a nimbleFunction.
  current_option <- getNimbleOption("stopOnSizeErrors")
  nimbleOptions(stopOnSizeErrors = FALSE)
  on.exit( nimbleOptions(stopOnSizeErrors = current_option ) )
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

  # January 2024, branch stop-on-size-errors: With the nimbleOption shown, these can use nimPrint instead of nimStop
  # Hence we expect an error but we expect that to fail because execution does not stop but just emits a message.
  # Print this message so anyone running tests knows this is ok.
  cat("\nSeven run-time size error (really warning) messages known to occur here:\n")
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


test_that("Testing of error handling in indexing of undefined array", {
    mynf <- nimbleFunction(
        run = function() {
            x[1] <- 5
        })
    
    expect_error(cnf <- compileNimble(mynf), info = "has 'x' been created")
    mynf <- nimbleFunction(
        run = function() {
            x[1,2] <- 5
        })
    expect_error(cnf <- compileNimble(mynf), info = "has 'x' been created")
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
