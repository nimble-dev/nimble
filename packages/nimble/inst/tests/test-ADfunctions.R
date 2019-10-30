source(system.file(file.path('tests', 'AD_test_utils.R'), package = 'nimble'))
source(system.file(file.path('tests', 'AD_math_test_lists.R'), package = 'nimble'))
source(system.file(file.path('tests', 'AD_distribution_test_lists.R'), package = 'nimble'))
source(system.file(file.path('tests', 'AD_knownFailures.R'), package = 'nimble'))
nimbleOptions(experimentalEnableDerivs = TRUE)
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
      ), enableDerivs = list('testMethod')
    )
    ADfunInst <- ADfun1()
    x <- matrix(c(2, -2))
    Rderivs <- ADfunInst$run(x)
    temporarilyAssignInGlobalEnv(ADfunInst)  
    cADfunInst <- compileNimble(ADfunInst)
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
                
              ), enableDerivs = list('testMethod', 'testMethod2')
            )
            ADfunInst <- ADfun2()
            x <- c(1, 1)
            Rderivs <- ADfunInst$run(x)
            temporarilyAssignInGlobalEnv(ADfunInst)  
            cADfunInst <- compileNimble(ADfunInst)
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
              ), enableDerivs = list('testMethod')
            )
            ADfunInst <- ADfun3()
            x <- c(1, 1)
            Rderivs <- ADfunInst$run(x)
            temporarilyAssignInGlobalEnv(ADfunInst)  
            cADfunInst <- compileNimble(ADfunInst)
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
                                     inits = list(x = rep(0, 2)),
                                     constants = list(zeroVec = rep(0, 2),
                                                      diagMat = diag(2)))
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
              ), enableDerivs = list('testMethod')
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



test_that('Derivatives of matrix exponentiation function correctly.',
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
              ), enableDerivs = list('testMethod')
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


## one for exp(mat)
## get order to owrk without c()!

#######################
## run all of the tests
#######################

## from AD_math_test_lists.R
test_AD_batch(unaryOpTests, knownFailures = AD_knownFailures)
test_AD_batch(unaryReductionOpTests, knownFailures = AD_knownFailures)
test_AD_batch(binaryOpTests, knownFailures = AD_knownFailures)
test_AD_batch(powOpTests, knownFailures = AD_knownFailures)
test_AD_batch(binaryReductionOpTests, knownFailures = AD_knownFailures)
test_AD_batch(squareMatrixOpTests, knownFailures = AD_knownFailures) ## These don't work and aren't in AD_knownFailures:chol, inverse (3 works, 4 fails), det (5 works, 6 fails), log det not sure (8 gets a NaN), trace (which is in knownFailures for compilation failure).
test_AD_batch(binaryMatrixOpTests, knownFailures = AD_knownFailures)

## from AD_distribution_test_lists.R
test_AD_batch(distn_tests,  knownFailures = AD_knownFailures) ## 56 is known failure that didn't fail. 91 failed, huge numbers.
nimbleOptions(pauseAfterWritingFiles = TRUE)
test_AD_batch(distn_with_log_tests[27], knownFailures = AD_knownFailures)
tempdir()
test_AD_batch(distn_tests[141])

test_AD_batch(distn_tests[141], knownFailures = AD_knownFailures)
test_AD_batch(distn_with_log_tests[141], knownFailures = AD_knownFailures)
debug(test_AD)
test_AD_batch(distn_with_log_tests[141])

idexp <- which(grepl("exp_R", names(distn_tests)))
test_AD_batch(distn_tests[idexp])
test_AD_batch(distn_with_log_tests[idexp])

idexp <- which(grepl("exp_nimble", names(distn_tests)))
test_AD_batch(distn_tests[idexp], knownFailures = AD_knownFailures)
test_AD_batch(distn_with_log_tests[idexp], knownFailures = AD_knownFailures)

idig <- which(grepl("invgamma", names(distn_tests)))
test_AD_batch(distn_tests[idig], knownFailures = AD_knownFailures)
test_AD_batch(distn_tests[idig[2]], knownFailures = AD_knownFailures)
test_AD_batch(distn_with_log_tests[idig], knownFailures = AD_knownFailures)

id <- which(grepl("sqrtinvgamma", names(distn_tests))) ## Something weird happening
test_AD_batch(distn_tests[id[1]])
test_AD_batch(distn_with_log_tests[id], knownFailures = AD_knownFailures)
test_AD_batch(distn_tests[id], knownFailures = AD_knownFailures)
test_AD_batch(distn_with_log_tests[id], knownFailures = AD_knownFailures)

id <- which(grepl("^gamma", names(distn_tests)))
test_AD_batch(distn_tests[id], knownFailures = AD_knownFailures)
test_AD_batch(distn_with_log_tests[id], knownFailures = AD_knownFailures)

id <- which(grepl("^t_R", names(distn_tests))) ## Review what happens here.
test_AD_batch(distn_tests[id])

id <- which(grepl("dlnorm", names(distn_tests))) ## Review what happens here.
test_AD_batch(distn_tests[id])
test_AD_batch(distn_tests[88])
test_AD_batch(distn_tests[89:92])
test_AD_batch(distn_with_log_tests[id])

id <- which(grepl("dlogis", names(distn_tests))) ## Review what happens here.
test_AD_batch(distn_tests[id])
test_AD_batch(distn_with_log_tests[id])

id <- which(grepl("dunif", names(distn_tests))) ## Review what happens here.
test_AD_batch(distn_tests[id])
undebug(test_AD)
test_AD_batch(distn_tests[136])
debug(test_AD)
test_AD_batch(distn_with_log_tests[id])
