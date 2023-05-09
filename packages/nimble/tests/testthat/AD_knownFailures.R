## try to keep in alphabetical order
AD_knownFailures <- list( ## not complete:
  ## name the elements by operator name
  ##
  "%%" = list(
    '*' = list( ## * means all inputs fail
      compilation = TRUE
    )
  )
  ## var = list(
  ##   'arg1 = double(2))' = list(
  ##     ## nimble doesn't support var of matrices, see sizeUnaryReduction
  ##     ## This has been removed from unaryReductionOpTests2
  ##     compilation = TRUE
  ##   )
  ## ),
  )
