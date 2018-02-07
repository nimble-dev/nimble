source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
nimbleOptions(experimentalEnableDerivs = TRUE)
nimbleOptions(showCompilerOutput = TRUE)
context("Testing of derivatives for distributions and dsl functions")

# untested:
# 'dcat', 'ddirch',
# 'dmulti', 'dmvt', 'dnegbin', 'dnorm', 'dpois', 
# 'dt', 'dunif', 'dweib', 'dwish')

# tested:
# 'dbeta' (although boundary at x=0 and x=1 not consistent between R and c++),
# 'dbinom' (works as long as wrt = prob and not x or size),
# 'dchisq' (inconsistencies at df = 0, incorrect R derivs at x = 0),
# 'dexp',
# 'dgamma',
# 'dinvgamma',
# 'dlogis',
# 'dlnorm',
# 'dmnorm',

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
  WRT = c('prob') ## only take derivs wrt 'prob' argument, otherwise warnings
                  ## from R, and inconsistencies between R and C++ will occur.
)
distributionArgsList[['dchisq']] <- list(
  distnName = 'dchisq',
  args = list(x = quote(double(0)),
              df = quote(double(0))),
  argsValues = list(
    list(x = -1, df = 1),
    # list(x = 0, df = 1),  #R derivs here are incorrect
    list(x = 1, df = -1),
    # list(x = 1, df = 0), #slight inconsistency between R and C++
    list(x = 1, df = 0.5),
    list(x = 1, df = 2),
    list(x = 12.5, df =22)
  )
)

distributionArgsList[['dexp']] <- list(
  distnName = 'dexp',
  args = list(x = quote(double(0)),
              rate = quote(double(0))),
  argsValues = list(
    list(x = -1, rate = 1),
    # list(x = -1, rate = 0), inconsistency between R and C++
    list(x = -1, rate = -1),
    # list(x = 0, rate = 1),  R derivs incorrect here
    # list(x = 0, rate = 0),  R derivs incorrect here
    # list(x = 0, rate = -1), R derivs incorrect here
    # list(x = 1, rate = 0), R derivs incorrect here
    list(x = 1, rate = 1),
    list(x = 12, rate = 13.2)
  )
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

distributionArgsList[['dinvgamma']] <- list(
  distnName = 'dinvgamma',
  args = list(x = quote(double(0)),
              shape = quote(double(0)),
              scale = quote(double(0))),
  argsValues = list(
    list(x = -1, shape = 1, scale = 1),
    list(x = -1, shape = -1, scale = 1),
    list(x = -1, shape = 1, scale = -1),
    list(x = .1, shape = 1, scale = 1),
    list(x = 22.2, shape = 10, scale = 15.2))
)

distributionArgsList[['dlogis']] <- list(
  distnName = 'dlogis',
  args = list(x = quote(double(0)),
              location = quote(double(0)),
              scale = quote(double(0))),
  argsValues = list(
    list(x = -1, location = 1, scale = 0),
    list(x = -1, location = -1, scale = -1),
    list(x = -1, location = 1, scale = 1),
    list(x = .1, location = 1, scale = 2),
    list(x = 22.2, location = 10, scale = 15.2))
)

distributionArgsList[['dlnorm']] <- list(
  distnName = 'dlnorm',
  args = list(x = quote(double(0)),
              meanlog = quote(double(0)),
              taulog = quote(double(0))),
  argsValues = list(
    # list(x = -1, meanlog = 1, taulog = 0),
    list(x = -1, meanlog = -1, taulog = -1),
    list(x = -1, meanlog = 1, taulog = 1),
    list(x = .1, meanlog = 1, taulog = 2),
    list(x = 22.2, meanlog = 10, taulog = 15.2))
)
runFun <- gen_runFunCore(makeADDistributionTestList(distributionArgsList[['dlnorm']]))
methodFun <- gen_runFunCore(makeADDistributionMethodTestList(distributionArgsList[['dlnorm']]))
thisNf <- nimbleFunction(setup = function(){},
                      run = runFun,
                      methods = list(
                        method1 = methodFun
                      ),
                      enableDerivs = list('method1'))
testADDistribution(thisNf, distributionArgsList[['dlnorm']]$argsValues, distributionArgsList[['dlnorm']]$distnName)

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
