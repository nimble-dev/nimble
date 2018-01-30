source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
nimbleOptions(experimentalEnableDerivs = TRUE)
nimbleOptions(showCompilerOutput = TRUE)
context("Testing of derivatives for distributions and dsl functions")

# untested:
# 'dbeta', 'dbin', 'dcat', 'dchisq', 'ddirch',
# 'dexp', 'dinvgamma', 'dlogis', 'dlnorm',
# 'dmulti', 'dmnorm', 'dmvt', 'dnegbin', 'dnorm', 'dpois', 
# 'dt', 'dunif', 'dweib', 'dwish')

# tested:
# 'dbeta' (although boundary at x=0 and x=1 not consistent between R and c++),
# 'dgamma',

distributionTests <- list()

distributionArgsList <- list()

distributionArgsList[['dbeta']] <- list(
  distnName = 'dbeta',
  args = list(x = quote(double(0)),
              shape1 = quote(double(0)),
              shape2 = quote(double(0))),
  argsValues = list(list(x = -1, shape1 = 1, shape2 = 1), ## x below support
                   # list(x = 0, shape1 = 1, shape2 = 1), ## x at boundary of support, currently not working
                    # list(x = 1, shape1 = 1, shape2 = 1), ## x at boundary of support, currently not working
                    list(x = 2, shape1 = 1,
                         shape2 = 1),
                    list(x = .3,  ## first param not valid
                         shape1 = -1,
                         shape2 = 1),
                    list(x = .3,  ## second param not valid
                         shape1 = 1,
                         shape2 = -1),
                    list(x = .9,  ## okay
                         shape1 = 12,
                         shape2 = .1),
                    list(x = .3,  ## okay
                         shape1 = 10,
                         shape2 = 1))
)
distributionArgsList[['dgamma']] <- list(
  distnName = 'dgamma',
  args = list(x = quote(double(0)),
              shape = quote(double(0)),
              rate = quote(double(0))),
  argsValues = list(
    list(x = -1, shape = 1, rate = 1),
    list(x = -1, shape = -1, rate = 1),
    list(x = -1, shape = 1, rate = -1),
    list(x = .1, shape = 1, rate = 1),
    list(x = 22.2, shape = 10, rate = 15.2))
)

runFun <- gen_runFunCore(makeADDistributionTestList(distributionArgsList[[1]]))
methodFun <- gen_runFunCore(makeADDistributionMethodTestList(distributionArgsList[[1]]))
thisNf <- nimbleFunction(setup = function(){},
                      run = runFun,
                      methods = list(
                        method1 = methodFun
                      ),
                      enableDerivs = list('method1'))
testADDistribution(thisNf, distributionArgsList[[1]]$argsValues, distributionArgsList[[1]]$distnName)

testFn <- function(x, shape1, shape2){
  - lgamma(shape1) - lgamma(shape2) ;
}
nimDerivs(testFn(.9, 12, .1))

Field "value":
  [1]   0.00000 -19.75502
Field "gradient":
  [,1]      [,2]     [,3]
[1,]    0  0.000000 0.000000
[2,]    0 -2.442662 7.981093