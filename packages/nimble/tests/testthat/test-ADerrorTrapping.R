source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)

test_that("Warning message works for call not supported for derivs.", {
  expect_message(
    ADfun1 <- nimbleFunction(
      setup = function(){},
      run = function(y = double(1)) {
        outList <- derivs(testMethod(y), wrt = c('x'))
        returnType(ADNimbleList())
        return(outList)
      },
      methods = list(
        testMethod = function(x = double(1, 2)) {
          out <- pnorm(x[1],0,1) ## Not supported
          a <- nimStep(x[1]) ## not supported
          returnType(double())
          return(out)
        }
      ), buildDerivs = c('testMethod')
    )
  )

  output <- capture_messages(
    ADfun1 <- nimbleFunction(
      setup = function(){},
      run = function(y = double(1)) {
        outList <- derivs(testMethod(y), wrt = c('x'))
        returnType(ADNimbleList())
        return(outList)
      },
      methods = list(
        testMethod = function(x = double(1, 2)) {
          out <- dnorm(x[1],0,1) ## supported
          returnType(double())
          return(out)
        }
      ), buildDerivs = c('testMethod')
    )
  )

  expect_identical(output, character())
})

test_that("Warning messages work for checking if a user-defined distribution supports derivs when it must.", {
  # Case where both user-defined dist and model have buildDerivs TRUE, so no warning is emitted.
  dmyexp <- nimbleFunction(
    run = function(x = double(0), rate = double(0, default = 1), log = integer(0, default = 0)) {
      returnType(double(0))
      logProb <- log(rate) - x*rate
      if(log) {
        return(logProb)
      } else {
        return(exp(logProb))
      }
    }, buildDerivs = TRUE)
  temporarilyAssignInGlobalEnv(dmyexp)
  
  rmyexp <- nimbleFunction(
    run = function(n = integer(0), rate = double(0, default = 1)) {
      returnType(double(0))
      if(n != 1) nimPrint("rmyexp only allows n = 1; using n = 1.")
      dev <- runif(1, 0, 1)
      return(-log(1-dev) / rate)
    }
  )
  temporarilyAssignInGlobalEnv(rmyexp, replace = TRUE)
  
  code1 <- nimbleCode({
    for(i in 1:3) {
      y1[i] ~ dmyexp(rate = r1)
    }
    r1 <- 1 / s1
    s1 ~ dunif(0, 100)
  })
  
  msgs <- capture_messages(
    m <- nimbleModel(code1, buildDerivs = TRUE)
  )

  expect_false(any(grepl("deriv", msgs)))

  # Case where model has buildDerivs FALSE, so no warning is emitted.
  ## deregisterDistributions("dmyexp")
  msgs <- capture_messages(
    m <- nimbleModel(code1)
  )
  expect_false(any(grepl("deriv", msgs)))

  # Case where user-defined dist has buildDerivs FALSE (by omission) but model needs derivs, so a warning is emitted.
  ## deregisterDistributions("dmyexp")
  dmyexp <- nimbleFunction(
    run = function(x = double(0), rate = double(0, default = 1), log = integer(0, default = 0)) {
      returnType(double(0))
      logProb <- log(rate) - x*rate
      if(log) {
        return(logProb)
      } else {
        return(exp(logProb))
      }
    })
  temporarilyAssignInGlobalEnv(rmyexp, replace = TRUE)

  msgs <- capture_messages(
    m <- nimbleModel(code1, buildDerivs = TRUE)
  )
  expect_true(any(grepl("Distribution dmyexp", msgs) &
                    grepl("derivatives", msgs)))



# Follow-on case where model does not need derivatives so no message is needed.
  msgs <- capture_messages(
    m <- nimbleModel(code1, buildDerivs = FALSE)
  )
  expect_false(any(grepl("Distribution dmyexp", msgs) &
                     grepl("derivatives", msgs)))

  # Case where user-defined dist has buildDerivs FALSE (by literal) but model needs derivs, so a warning is emitted.
  ## deregisterDistributions("dmyexp")
  dmyexp <- nimbleFunction(
    run = function(x = double(0), rate = double(0, default = 1), log = integer(0, default = 0)) {
      returnType(double(0))
      logProb <- log(rate) - x*rate
      if(log) {
        return(logProb)
      } else {
        return(exp(logProb))
      }
    }, buildDerivs = FALSE)
  temporarilyAssignInGlobalEnv(rmyexp, replace = TRUE)
  msgs <- capture_messages(
    m <- nimbleModel(code1, buildDerivs = TRUE)
  )
  expect_true(any(grepl("Distribution dmyexp", msgs) &
                    grepl("derivatives", msgs)))
  # Follow-on case where model does not need derivatives so no message is needed.
  msgs <- capture_messages(
    m <- nimbleModel(code1, buildDerivs = FALSE)
  )
  expect_false(any(grepl("Distribution dmyexp", msgs) &
                     grepl("derivatives", msgs)))

  # Case where user-defined dist has buildDerivs FALSE (by empty list) but model needs derivs, so a warning is emitted.
  ## deregisterDistributions("dmyexp")
  dmyexp <- nimbleFunction(
    run = function(x = double(0), rate = double(0, default = 1), log = integer(0, default = 0)) {
      returnType(double(0))
      logProb <- log(rate) - x*rate
      if(log) {
        return(logProb)
      } else {
        return(exp(logProb))
      }
    }, buildDerivs = list())
  temporarilyAssignInGlobalEnv(rmyexp, replace = TRUE)
  msgs <- capture_messages(
    m <- nimbleModel(code1, buildDerivs = TRUE)
  )
  expect_true(any(grepl("Distribution dmyexp", msgs) &
                    grepl("derivatives", msgs)))
  # Follow-on case where model does not need derivatives so no message is needed.
  msgs <- capture_messages(
    m <- nimbleModel(code1, buildDerivs = FALSE)
  )
  expect_false(any(grepl("Distribution dmyexp", msgs) &
                     grepl("derivatives", msgs)))

})

test_that("Warning message works for buildDerivs not set for methods being differentiated.", {

    expect_silent(
        derivs_nf <- nimbleFunction(
            setup = function(model, with_respect_to_nodes, calc_nodes) {},
            run = function(order = integer(1),
                           reset = logical(0, default=FALSE)) {
                ans <- nimDerivs(model$calculate(calc_nodes), wrt = with_respect_to_nodes,
                                 order = order, reset = reset)
                return(ans)
                returnType(ADNimbleList())
            }
        )
    )
    
    expect_silent(
        ADfun1 <- nimbleFunction(
            setup = function(){},
            run = function(y = double(1)) {
                outList <- derivs(testMethod(y), wrt = c('x'))
                returnType(ADNimbleList())
                return(outList)
            },
            methods = list(
                testMethod = function(x = double(1, 2)) {
                    out <- dnorm(x[1],0,1)
                    returnType(double())
                    return(out)
                }
            ), buildDerivs = 'testMethod'
        )
    )

    expect_silent(
        ADfun1 <- nimbleFunction(
            setup = function(){},
            run = function(y = double(1)) {
                outList <- derivs(testMethod(y), wrt = c('x'))
                returnType(ADNimbleList())
                return(outList)
            },
            methods = list(
                testMethod = function(x = double(1, 2)) {
                    out <- dnorm(x[1],0,1)
                    returnType(double())
                    return(out)
                }
            ), buildDerivs = list(testMethod = list(ignore = 'foo'))
        )
    )

    expect_message(
        ADfun1 <- nimbleFunction(
            setup = function(){},
            run = function(y = double(1)) {
                outList <- derivs(testMethod(y), wrt = c('x'))
                returnType(ADNimbleList())
                return(outList)
            },
            methods = list(
                testMethod = function(x = double(1, 2)) {
                    out <- dnorm(x[1],0,1)
                    returnType(double())
                    return(out)
                }
            )
        ),
        "Detected use of `nimDerivs`")

    expect_message(
        ADfun1 <- nimbleFunction(
            setup = function(){},
            run = function(y = double(1)) {
                outList <- derivs(testMethod(y), wrt = c('x'))
                returnType(ADNimbleList())
                return(outList)
            },
            methods = list(
                testMethod = function(x = double(1, 2)) {
                    out <- dnorm(x[1],0,1)
                    returnType(double())
                    return(out)
                }
            ), buildDerivs = 'run'
        ),
        "Detected use of `nimDerivs`")
    
    expect_message(
        ADfun1 <- nimbleFunction(
            setup = function(){},
            run = function(y = double(1)) {
                outList <- derivs(testMethod(y), wrt = c('x'))
                returnType(ADNimbleList())
                return(outList)
            },
            methods = list(
                testMethod = function(x = double(1, 2)) {
                    out <- dnorm(x[1],0,1)
                    returnType(double())
                    return(out)
                }
            ), buildDerivs = list(run = list(ignore = 'foo'))
        ),
        "Detected use of `nimDerivs`")
})


test_that("Incorrect use of buildDerivs=TRUE in nimbleFunction with setup.", {
    expect_error(
        nf_sqrt <- nimbleFunction(
            setup = function() {},
            run = function(x = double(1)) {
                return(sqrt(x))
                returnType(double(1))
            },
            buildDerivs = TRUE
        ), "'buildDerivs' cannot be 'TRUE' when a setup function is provided"
    )
})

