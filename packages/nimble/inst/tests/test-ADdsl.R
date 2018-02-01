source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
nimbleOptions(experimentalEnableDerivs = TRUE)
nimbleOptions(showCompilerOutput = TRUE)
context("Testing of derivatives for distributions and dsl functions")

# untested:
# 'dbin', 'dcat', 'dchisq', 'ddirch',
# 'dexp', 'dinvgamma', 'dlogis', 'dlnorm',
# 'dmulti', 'dmnorm', 'dmvt', 'dnegbin', 'dnorm', 'dpois', 
# 'dt', 'dunif', 'dweib', 'dwish')

# tested:
# 'dbeta' (although boundary at x=0 and x=1 not consistent between R and c++),
# 'dbinom' (works as long as wrt = prob and not x or size),
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
distributionArgsList[['dbinom']] <- list(
  distnName = 'dbinom',
  args = list(x = quote(double(0)),
              prob = quote(double(0)),
              size = quote(double(0))),
  argsValues = list(
    list(x = -1, prob = .1, size = 1),
    list(x = 1.5, prob = .1, size = 1),
    list(x = 1.5, prob = .1, size = 2),
    list(x = 2, prob = .1, size = 1),
    list(x = 2, prob = .1, size = 2.5),
    list(x = 2, prob = -1, size = 3),
    list(x = 2, prob = 10, size = 3),
    list(x = 2, prob = 0, size = 3),
    list(x = 2, prob = 1, size = 3),
    list(x = 2, prob = .1, size = 3),
    list(x = 25, prob = .5, size = 50)),
  WRT = c('prob')
)

runFun <- gen_runFunCore(makeADDistributionTestList(distributionArgsList[[3]]))
methodFun <- gen_runFunCore(makeADDistributionMethodTestList(distributionArgsList[[3]]))
thisNf <- nimbleFunction(setup = function(){},
                      run = runFun,
                      methods = list(
                        method1 = methodFun
                      ),
                      enableDerivs = list('method1'))
testADDistribution(thisNf, distributionArgsList[[3]]$argsValues, distributionArgsList[[3]]$distnName)

testFn <- function(x, shape1, shape2){
  - lgamma(shape1) - lgamma(shape2) ;
}
nimDerivs(dbinom(1,2,1),wrt = c('prob'))


testMod <- nimbleCode({
  a ~ dbin(.1, 5)
})

testMod1 <- nimbleModel(testMod)
ctestMod1 <- compileNimble(testMod1)

as.call(list(quote(dbinom), x = 1, size = 2, prob = .1))
test <- list(a = 1, b =2)
lapply(test, function(x){return(x[[1]])})
