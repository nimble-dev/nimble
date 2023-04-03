source(system.file(file.path('tests', 'testthat', 'AD_test_utils.R'), package = 'nimble'))
source(system.file(file.path('tests', 'testthat', 'AD_math_test_lists.R'), package = 'nimble'))
source(system.file(file.path('tests', 'testthat', 'AD_distribution_test_lists.R'), package = 'nimble'))
source(system.file(file.path('tests', 'testthat', 'AD_knownFailures.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)
nimbleOptions(allowDynamicIndexing = FALSE)

context("Testing of derivatives for nimbleFunctions.")

test_that('Derivatives of dnorm function correctly.',
  {
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
      ), buildDerivs = c('testMethod')
    )
    ADfunInst <- ADfun1()
    xRec <- matrix(c(1, -1))
    x <- matrix(c(2, -2))
    Rderivs <- ADfunInst$run(x)
    temporarilyAssignInGlobalEnv(ADfunInst)  
    cADfunInst <- compileNimble(ADfunInst)
    cADfunInst$run(xRec)
    cderivs <- cADfunInst$run(x)
    expect_equal(cderivs$value, Rderivs$value)
    expect_equal(cderivs$jacobian, Rderivs$jacobian)
    expect_equal(cderivs$hessian, Rderivs$hessian)
  }
) 

test_that('Derivatives of x^2 function correctly.',
          {
            listOfADLists <- nimbleList(list1 = ADNimbleList(),
                                        list2 = ADNimbleList())
            temporarilyAssignInGlobalEnv(listOfADLists)  
            ADfun2 <- nimbleFunction(
              setup = function(){},
              run = function(y = double(1, c(2))) {
                ADlistList <- listOfADLists$new()
                ADlistList$list1 <- derivs(testMethod(y, y))
                ADlistList$list2 <- derivs(testMethod2(y, y))
                returnType(listOfADLists())
                return(ADlistList)
              },
              methods = list(
                testMethod = function(x = double(1, c(2)), y = double(1, c(2))) {
                  returnType(double(1, c(2)))
                  return(x[]^2)
                },
                testMethod2 = function(x = double(1, c(2)), y = double(1, c(2))) {
                  returnType(double(0))
                  return(2^x[1])
                }
                
              ), buildDerivs = c('testMethod', 'testMethod2')
            )
            ADfunInst <- ADfun2()
            xRec <- c(2, 2)
            x <- c(1, 1)
            Rderivs <- ADfunInst$run(x)
            temporarilyAssignInGlobalEnv(ADfunInst)  
            suppressWarnings(cADfunInst <- compileNimble(ADfunInst))
            cADfunInst$run(xRec)
            cderivs <- cADfunInst$run(x)
            expect_equal(cderivs$list1$value, Rderivs$list1$value)
            expect_equal(cderivs$list1$jacobian, Rderivs$list1$jacobian, tolerance = 0.01)
            expect_equal(cderivs$list1$hessian, Rderivs$list1$hessian, tolerance = 0.01)
            expect_equal(cderivs$list2$value, Rderivs$list2$value)
            expect_equal(cderivs$list2$jacobian, Rderivs$list2$jacobian, tolerance = 0.01)
            expect_equal(cderivs$list2$hessian, Rderivs$list2$hessian, tolerance = 0.01)
          }
) 

test_that('Derivatives of sum(log(x)) function correctly.',
          {
            ADfun3 <- nimbleFunction(
              setup = function(){},
              run = function(x = double(1, c(2))) {
                outVals <- 2*derivs(testMethod(x))$jacobian
                returnType(double(2))
                return(outVals)
              },
              methods = list(
                testMethod = function(x = double(1, c(2))) {
                  returnType(double(0))
                  return(sum(log(x[])))
                }
              ), buildDerivs = c('testMethod')
            )
            ADfunInst <- ADfun3()
            xRec <- c(1, 1)
            x <- c(1, 1)
            Rderivs <- ADfunInst$run(x)
            temporarilyAssignInGlobalEnv(ADfunInst)  
            cADfunInst <- compileNimble(ADfunInst)
            cADfunInst$run(xRec)
            cderivs <- cADfunInst$run(x)
            expect_equal(cderivs, Rderivs, tolerance = 0.01)
          }
) 

test_that('Derivatives of model$calculate() work in expressions.', 
          {
            testModelCode <- nimbleCode({
              x[1:2] ~ dmnorm(zeroVec[1:2], diagMat[1:2, 1:2])
            })
            testModel <- nimbleModel(code = testModelCode,
                                     inits = list(x = c(0.5, 0.6)),
                                     constants = list(zeroVec = c(.3, .4),
                                                      diagMat = matrix(c(1.2, .2, .2, 1.6), nrow = 2)))
            ADfun4 <- nimbleFunction(
              setup = function(model){},
              run = function() {
                outVals <- 2*sum(derivs(model$calculate('x'), 
                                        wrt = 'x')$jacobian)
                returnType(double())
                return(outVals)
              }
            )
            ADfunInst <- ADfun4(testModel)
            Rderivs <- ADfunInst$run()
            temporarilyAssignInGlobalEnv(testModel)  
            temporarilyAssignInGlobalEnv(ADfunInst)  
            compileNimble(testModel)
            cADfunInst <- compileNimble(ADfunInst)
            cderivs <- cADfunInst$run()
            expect_equal(cderivs, Rderivs, tolerance = 0.01)
          }
) 

test_that('Derivatives of matrix multiplication function correctly.',
          {
            ADfun5 <- nimbleFunction(
              setup = function(){},
              run = function(x = double(2, c(2, 2)),
                             y = double(2, c(2, 4))) {
                outVals <- derivs(testMethod(x, y), wrt = 'x', order = 1)$jacobian
                returnType(double(2))
                return(outVals)
              },
              methods = list(
                testMethod = function(x = double(2, c(2, 2)),
                                      y = double(2, c(2, 4))) {
                  returnType(double(2, c(2, 4)))
                  return(x%*%y)
                }
              ), buildDerivs = c('testMethod')
            )
            ADfunInst <- ADfun5()
            x <- diag(2)
            set.seed(1)
            y <- matrix(rnorm(8), ncol = 4, nrow = 2)
            Rderivs <- ADfunInst$run(x, y)
            temporarilyAssignInGlobalEnv(ADfunInst)  
            cADfunInst <- compileNimble(ADfunInst)
            cderivs <- cADfunInst$run(x, y)
            expect_equal(cderivs, Rderivs, tolerance = 0.01)
          }
)

test_that('Derivatives of matrix component-wise exponentiation function correctly.',
          {
            ADfun6 <- nimbleFunction(
              setup = function(){},
              run = function(x = double(2, c(2, 2))) {
                outVals <- derivs(testMethod(x), wrt = 'x', order = 1)$jacobian
                returnType(double(2))
                return(outVals)
              },
              methods = list(
                testMethod = function(x = double(2, c(2, 2))) {
                  returnType(double(2, c(2, 2)))
                  return(exp(x))
                }
              ), buildDerivs = c('testMethod')
            )
            ADfunInst <- ADfun6()
            x <- diag(2)
            Rderivs <- ADfunInst$run(x)
            temporarilyAssignInGlobalEnv(ADfunInst)  
            cADfunInst <- compileNimble(ADfunInst)
            cderivs <- cADfunInst$run(x)
            expect_equal(cderivs, Rderivs, tolerance = 0.01)
          }
) 

test_that('Derivatives with rep() work correctly.',
{
  ADfun7 <- nimbleFunction(
    setup = TRUE,
    run = function(x = double(), n = double()) {
      y <- rep(x, n) # Regardless of n or n2, the second argument is lifted to a TYPE_
      ans <- sum(y^2)
      return(ans)
      returnType(double())
    },
    methods = list(
      derivsRun = function(x = double(), n = double()) {
        ans <- derivs(run(x, n), wrt = c(1), order = 0:2)
        return(ans)
        returnType(ADNimbleList())
      }),
    buildDerivs = list(run = list(ignore = c("n")))
  )

  ADfunInst <- ADfun7()
  temporarilyAssignInGlobalEnv(ADfunInst)
  cADfunInst <- compileNimble(ADfunInst)
  cderivs <- cADfunInst$derivsRun(3,4)
  derivs <- ADfunInst$derivsRun(3,4)
  expect_equivalent(cderivs$value, derivs$value)
  expect_equivalent(cderivs$jacobian, derivs$jacobian)
  expect_equivalent(cderivs$hessian, derivs$hessian)
}
) 

test_that('Derivatives with c() work correctly.',
{
  ADfun8 <- nimbleFunction(
    setup = TRUE,
    run = function(x = double(), y = double()) {
      z <- c(x, y)
      ans <- sum(z^2)
      return(ans)
      returnType(double())
    },
    methods = list(
      derivsRun = function(x = double(), y = double()) {
        ans <- derivs(run(x, y), wrt = 1:2, order = 0:2)
        return(ans)
        returnType(ADNimbleList())
      }),
    buildDerivs = list(run = list())
  )

  ADfunInst <- ADfun8()
  temporarilyAssignInGlobalEnv(ADfunInst)
  cADfunInst <- compileNimble(ADfunInst)
  cderivs <- cADfunInst$derivsRun(3,5)
  derivs <- ADfunInst$derivsRun(3,5)
  expect_equivalent(cderivs$value, derivs$value)
  expect_equivalent(cderivs$jacobian, derivs$jacobian)
  expect_equivalent(cderivs$hessian, derivs$hessian)
}
) 

## Dirichlet with log argument
dirch_test_log <- make_AD_test2(
  op = list(
    name = "ddirch manual using additive log-ratio transformation",
    opParam = list(name = "ddirch manual"),
    # X = Xtrans = log(Xorig_i / Xorig_n), i = 1... n-1
    # Xorig[1:n-1] = exp(Xtrans)
    # Xorig[n] = 1/(1+sum(Xorig[1:n-1]))
    expr = quote({
      Xorig_1_nm1_over_Xn <- exp(x)
      Xorig_n <- 1/(1 + sum(Xorig_1_nm1_over_Xn))
      Xorig <- c(Xorig_1_nm1_over_Xn * Xorig_n, Xorig_n)
      out <- ddirch(x = Xorig, alpha=alpha, log = log)
    }),
    args = list(
      x = quote(double(1)),
      alpha = quote(double(1)),
      log = quote(double())
    ),
    outputType = quote(double())
  ),
  argTypes = c(x='double(1)', alpha='double(1)', log='double()'),
  wrt = c('x', 'alpha'),
  inputs = list(record = list(x = c(log(.2/.5), log(.3/.5)), alpha = c(2, 4, 4), log = 1),
                test   = list(x = c(log(.4/.45), log(.15/.45)), alpha = c(3, 2, 5), log = 0)))

dirch_test_out <- test_AD2(dirch_test_log) 

## Dirichlet without log argument
dirch_test_fixedlog <- make_AD_test2(
  op = list(
    name = "ddirch manual using additive log-ratio transformation",
    opParam = list(name = "ddirch manual"),
    # X = Xtrans = log(Xorig_i / Xorig_n), i = 1... n-1
    # Xorig[1:n-1] = exp(Xtrans)
    # Xorig[n] = 1/(1+sum(Xorig[1:n-1]))
    expr = quote({
      Xorig_1_nm1_over_Xn <- exp(x)
      Xorig_n <- 1/(1 + sum(Xorig_1_nm1_over_Xn))
      Xorig <- c(Xorig_1_nm1_over_Xn * Xorig_n, Xorig_n)
      out <- ddirch(x = Xorig, alpha=alpha, log = 1)
    }),
    args = list(
      x = quote(double(1)),
      alpha = quote(double(1))
    ),
    outputType = quote(double())
  ),
  argTypes = c(x='double(1)', alpha='double(1)'),
  wrt = c('x', 'alpha'),
  inputs = list(record = list(x = c(log(.2/.5), log(.3/.5)), alpha = c(2, 4, 4)),
                test   = list(x = c(log(.4/.45), log(.15/.45)), alpha = c(3, 2, 5))))

dirch_test_out <- test_AD2(dirch_test_fixedlog) 

## Wishart
## Testing Wishart is tricky because both x and cholesky are
## matrices and also because not all elements are indendent
makeARcov <- function(n, rho, sigma) {
  s2 <- sigma^2
  ans <- matrix(nrow = n, ncol = n)
  for(i in 1:n) {
    for(j in 1:i) {
      ans[i,j] <- ans[j,i] <- s2 * rho^(abs(i-j))
    }
  }
  ans
}
set.seed(123)
cholRec <- chol(makeARcov(4, .6, 2))
cholTest <- chol(makeARcov(4, .55, 3))

wRec <- rwish_chol(1, cholRec, df = 7, scale_param = FALSE)
wTest <- rwish_chol(1, cholTest, df = 7, scale_param = FALSE)

cholRecTri <- cholRec[ upper.tri(cholRec, TRUE)]
cholTestTri <-cholTest[ upper.tri(cholTest, TRUE)]

wRecTri <- wRec[ upper.tri(wRec, TRUE) ]
wTestTri <- wRec[ upper.tri(wTest, TRUE) ]

wish_test_log <- make_AD_test2(
  op = list(
    name = "dwish_chol manual",
    opParam = list(name = "dwish_chol manual"),
    expr = quote({
      # populate 2D matrices from the vectors
      # created from upper triangular values.
      x2D <- nimMatrix(nrow = 4, ncol = 4, init=FALSE)
      chol2D <- nimMatrix(nrow = 4, ncol = 4, init = FALSE)
      i <- 1L
      j <- 1L
      indOrig <- 1L
      for(j in 1:4) {
        for(i in 1:j) {
          x2D[i,j] <- x2D[j,i] <- x[indOrig]
          chol2D[i,j] <- chol[indOrig]  # chol2D does not need lower triangular entries
          indOrig <- indOrig + 1
        }
      }
      out <- dwish_chol(x = x2D, cholesky=chol2D, df = df, log = log)
    }),
    args = list(
      x = quote(double(1)),
      chol = quote(double(1)),
      df = quote(double()),
      log = quote(double())
    ),
    outputType = quote(double())
  ),
  argTypes = c(x='double(1)', chol='double(1)', df = 'double()', log='double()'),
  wrt = c('x', 'chol'),
  inputs = list(record = list(x = wRecTri, chol = cholRecTri, df = 7, log = 0),
                test   = list(x = wTestTri, chol = cholTestTri, df = 8, log = 1)),
)

wish_test_out <- test_AD2(wish_test_log)

## Tests above each take 10-20 seconds.

#######################
## run all of the tests
#######################

## from AD_math_test_lists.R
## See notes there on Spring 2022 update to "version 2" with f(g(x)) and g(f(x)) test variantes as well
test_AD_batch(unaryOpTests2, testFun = test_AD2, knownFailures = AD_knownFailures, verbose = FALSE)
test_AD_batch(unaryOpTests2_inner, testFun = test_AD2, knownFailures = AD_knownFailures, verbose = FALSE)
test_AD_batch(unaryOpTests2_outer, testFun = test_AD2, knownFailures = AD_knownFailures, verbose = FALSE)

test_AD_batch(unaryReductionOpTests2, testFun = test_AD2, knownFailures = AD_knownFailures, verbose = FALSE)
test_AD_batch(unaryReductionOpTests2_inner, testFun = test_AD2, knownFailures = AD_knownFailures, verbose = FALSE)
test_AD_batch(unaryReductionOpTests2_outer, testFun = test_AD2, knownFailures = AD_knownFailures, verbose = FALSE)

test_AD_batch(binaryOpTests2, testFun = test_AD2, knownFailures = AD_knownFailures, verbose = FALSE)
test_AD_batch(binaryOpTests2_inner, testFun = test_AD2, knownFailures = AD_knownFailures, verbose = FALSE)
test_AD_batch(binaryOpTests2_outer, testFun = test_AD2, knownFailures = AD_knownFailures, verbose = FALSE)

test_AD_batch(powOpTests2, testFun = test_AD2, knownFailures = AD_knownFailures, verbose = FALSE)
test_AD_batch(powOpTests2_inner, testFun = test_AD2, knownFailures = AD_knownFailures, verbose = FALSE)
test_AD_batch(powOpTests2_outer, testFun = test_AD2, knownFailures = AD_knownFailures, verbose = FALSE)

test_AD_batch(pow_int_OpTests2, knownFailures = AD_knownFailures, verbose = FALSE)
test_AD_batch(pow_int_OpTests2_inner, testFun = test_AD2, knownFailures = AD_knownFailures, verbose = FALSE)
test_AD_batch(pow_int_OpTests2_outer, testFun = test_AD2, knownFailures = AD_knownFailures, verbose = FALSE)

# inprod
test_AD_batch(binaryReductionOpTests2, testFun = test_AD2, knownFailures = AD_knownFailures, verbose = FALSE)
ADtestEnv$RCrelTol[4] <- 1e-6
test_AD_batch(binaryReductionOpTests2_inner, testFun = test_AD2, knownFailures = AD_knownFailures, verbose = FALSE)
test_AD_batch(binaryReductionOpTests2_outer, testFun = test_AD2, knownFailures = AD_knownFailures, verbose = FALSE)
resetTols()

test_AD_batch(squareMatrixOpTests2, testFun = test_AD2, knownFailure = AD_knownFailures, verbose=FALSE)
test_AD_batch(squareMatrixOpTests[9], verbose = FALSE) ## trace has a knownFailures entry for compilation failure.  Do we even support it?

test_AD_batch(squareMatrixOpTests, knownFailures = AD_knownFailures, verbose = FALSE) ## trace has a knownFailures entry for compilation failure.  Do we even support it?
test_AD_batch(binaryMatrixOpTests, knownFailures = AD_knownFailures, verbose = FALSE)

# test_AD_batch does not work with these.  Check that out
# make lapply run line for distn_tests2
# test_AD_batch(distn_with_log_tests2, knownFailures = AD_knownFailures, verbose = TRUE)

nimbleOptions(enableDerivs = EDopt)
nimbleOptions(buildModelDerivs = BMDopt)
