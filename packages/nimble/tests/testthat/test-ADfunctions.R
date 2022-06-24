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
            print(cderivs)
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

#######################
## run all of the tests
#######################

## from AD_math_test_lists.R
## See notes there on Spring 2022 update to "version 2" with f(g(x)) and g(f(x)) test variantes as well
test_AD_batch(unaryOpTests2, testFun = test_AD2, knownFailures = AD_knownFailures)
ADtestEnv$RCrelTol <- c(1e-15, 1e-4, 1e-2) ## Loosen tols because there are more operations
test_AD_batch(unaryOpTests2_inner, testFun = test_AD2, knownFailures = AD_knownFailures)
resetTols()
ADtestEnv$RCrelTol <- c(1e-15, 1e-6, 1e-2)
test_AD_batch(unaryOpTests2_outer, testFun = test_AD2, knownFailures = AD_knownFailures)
resetTols()



test_AD_batch(unaryReductionOpTests2, testFun = test_AD2, knownFailures = AD_knownFailures)
ADtestEnv$RCrelTol <- c(1e-12, 1e-5, 1e-3)
test_AD_batch(unaryReductionOpTests2_inner, testFun = test_AD2, knownFailures = AD_knownFailures)
resetTols()
ADtestEnv$RCrelTol <- c(1e-12, 1e-5, 1e-3)
test_AD_batch(unaryReductionOpTests2_outer, testFun = test_AD2, knownFailures = AD_knownFailures)
resetTols()


test_AD_batch(binaryOpTests2, testFun = test_AD2, knownFailures = AD_knownFailures)
ADtestEnv$RCrelTol <- c(1e-15, 1e-4, 1e-2) ## Loosen tols because there are more operations
test_AD_batch(binaryOpTests2_inner, testFun = test_AD2, knownFailures = AD_knownFailures)
test_AD_batch(binaryOpTests2_outer, testFun = test_AD2, knownFailures = AD_knownFailures)

resetTols()

test_AD_batch(powOpTests2, testFun = test_AD2, knownFailures = AD_knownFailures)
## STOPPED HERE

test_AD_batch(pow_int_OpTests, knownFailures = AD_knownFailures)
test_AD_batch(binaryReductionOpTests, knownFailures = AD_knownFailures)
test_AD_batch(squareMatrixOpTests, knownFailures = AD_knownFailures) ## trace has a knownFailures entry for compilation failure.  Do we even support it?
test_AD_batch(binaryMatrixOpTests, knownFailures = AD_knownFailures)

## from AD_distribution_test_lists.R
test_AD_batch(distn_tests,  knownFailures = AD_knownFailures) ## dcat are knownFailures.  
debug(test_AD)

test_AD_batch(distn_tests[9])#,  knownFailures = AD_knownFailures) ## dcat are knownFailures.  
test_AD_batch(distn_tests[13:14],  knownFailures = AD_knownFailures) ## dmulti fails, perhaps in regular compilation too?
## test_AD_batch(distn_tests[15:148],  knownFailures = AD_knownFailures) ## dmulti fails, perhaps in regular compilation too?
## test_AD_batch(distn_tests[39:64]) ## ddexp was in knownFailures.  Now it's ok.
## test_AD_batch(distn_tests[65:148]) ## There is a problem in sqrtinvgamma with scale, vectorized
## test_AD_batch(distn_tests[85:148]) ## There is a problem in dt with vector x
## test_AD_batch(distn_tests[113:148])  ## There is a problem in dt_nonstandard with vector x
## test_AD_batch(distn_tests[129:148]) ## ddirch fails
## test_AD_batch(distn_tests[146]) 
## test_AD_batch(distn_tests[147]) 
## test_AD_batch(distn_tests[148]) ## dwish fails


test_AD_batch(distn_with_log_tests,  knownFailures = AD_knownFailures)

nimbleOptions(enableDerivs = EDopt)
nimbleOptions(buildModelDerivs = BMDopt)
