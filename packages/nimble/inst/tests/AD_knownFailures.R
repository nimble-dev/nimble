## Just an idea for now. We'd use this list in test_AD rather than having the
## knownFailure information embedded in the test parameterization.

knownFailures <- list( ## not complete
  ##
  ## operator name
  ##
  exp = list(
    ##
    ## arg type when failure occured (copied from test param name
    ##
    "arg1 = double(2, c(3, 4))" = list(
      ##
      ## first part of output that led to failure (value, jacobian, or hessian)
      ##
      hessian = 'method1', ## give name of first method that led to failure
      notes = 'numerical precision error'
    )
  ),
  inprod = list(
  "arg1 = double(1, 4) arg2 = double(1, 4)" = list(
    compilation = TRUE ## specify a compilation failure
  )
)
