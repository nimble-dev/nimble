source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
nimbleOptions(experimentalEnableDerivs = TRUE)
nimbleOptions(showCompilerOutput = TRUE)
context("Testing of derivatives for distributions and dsl functions")


'dbern', 'dbeta', 'dbin', 'dcat', 'dchisq', 'ddirch',
'dexp', 'dgamma', 'dinvgamma', 'dlogis', 'dlnorm',
'dmulti', 'dmnorm', 'dmvt', 'dnegbin', 'dnorm', 'dpois', 
'dt', 'dunif', 'dweib', 'dwish')


distributionTests <- list()

distributionArgsList <- list()
distributionArgsList[[1]] <- list(
  distnName = 'dbern',
  args = list(x = quote(double(0)),
              prob = quote(double(0)))
)
distributionArgsList[[2]] <- list(
  distnName = 'dbeta',
  args = list(x = quote(double(0)),
              shape1 = quote(double(0)),
              shape2 = quote(double(0)))
)



runFun <- gen_runFunCore(makeADDistributionTestList(distributionArgsList[[2]]))
methodFun <- gen_runFunCore(makeADDistributionMethodTestList(distributionArgsList[[2]]))
thisNf <- nimbleFunction(setup = function(){},
                      run = runFun,
                      methods = list(
                        method1 = methodFun
                      ),
                      enableDerivs = list('method1'))
RthisNf <- thisNf()
CthisNf <- compileNimble(RthisNf)

