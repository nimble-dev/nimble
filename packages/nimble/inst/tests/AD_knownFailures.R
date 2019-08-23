## try to keep in alphabetical order
AD_knownFailures <- list( ## not complete:
  ##
  ## name the elements by operator name
  ##
  "%%" = list(
    '*' = list(
      compilation = TRUE
    )
  ),
  abs = list(
    "*" = list(
      input = list(
        cpp = 'negative values',
        R = NULL ## if the R version had errors too, we could say so here
      )
    )
  ),
  ceil = list(
    "*" = list(
      compilation = TRUE
    )
  ),
  exp = list(
    ##
    ## arg type when failure occured (copied from test param name)
    ##
    "arg1 = double(2, c(3, 4))" = list(
      ##
      ## first part of output that led to failure (value, jacobian, or hessian)
      ##
      hessian = 'method1', ## give name of first method that led to failure
      notes = 'numerical precision error'
    )
  ),
  floor = list(
    "*" = list(
      compilation = TRUE
    )
  ),
  ftrunc = list(
    "*" = list(
      compilation = TRUE
    )
  ),
  ## 'R implements lgamma as "the natural logarithm of _the absolute value of_
  ## the gamma function", and since we use exp(lgamma(x)) for gammafn(x) we get
  ## the sign of the output flipped when gamma(x) < 0.'
  factorial = list(
    "*" = list(
      input = list( ## currently doesn't do anything, just a reminder
        cpp = 'any input x where gamma(x) < 0, e.g. in (-1, 0)'
      )
    )
  ),
  gammafn = list(
    "*" = list(
      notes = 'See entry for factorial.'
    )
  ),
  inprod = list(
    "arg1 = double(1, 4) arg2 = double(1, 4)" = list(
      compilation = TRUE ## specify a compilation failure
    )
  ),
  nimRound = list(
    "*" = list(
      compilation = TRUE
    )
  ),
  var = list(
    'arg1 = double(2, c(3, 4))' = list(
      compilation = TRUE,
      notes = "nimble doesn't support var of matrices, see sizeUnaryReduction"
    )
  )
)
