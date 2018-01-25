source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
nimbleOptions(experimentalEnableDerivs = TRUE)
nimbleOptions(showCompilerOutput = TRUE)
context("Testing of derivatives for distributions and dsl functions")


# 'dbern', 'dbeta', 'dbin', 'dcat', 'dchisq', 'ddirch',
# 'dexp', 'dgamma', 'dinvgamma', 'dlogis', 'dlnorm',
# 'dmulti', 'dmnorm', 'dmvt', 'dnegbin', 'dnorm', 'dpois', 
# 'dt', 'dunif', 'dweib', 'dwish')


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
              shape2 = quote(double(0))),
  argsValues = list(x = .1,
                    shape1 = 1,
                    shape2 = 1)
)
distributionArgsList[[3]] <- list(
  distnName = 'dgamma',
  args = list(x = quote(double(0)),
              shape = quote(double(0)),
              rate = quote(double(0))),
  argsValues = list(x = .1,
                    shape = 1,
                    rate = 1)
)

runFun <- gen_runFunCore(makeADDistributionTestList(distributionArgsList[[3]]))
methodFun <- gen_runFunCore(makeADDistributionMethodTestList(distributionArgsList[[3]]))
thisNf <- nimbleFunction(setup = function(){},
                      run = runFun,
                      methods = list(
                        method1 = methodFun
                      ),
                      enableDerivs = list('method1'))
# testADDistribution(thisNf, distributionArgsList[[2]]$argsValues, distributionArgsList[[2]]$distnName)

RthisNf <- thisNf()
CthisNf <- compileNimble(RthisNf)
round(CthisNf$run(1, 1, 1)$hessian -RthisNf$run(1, 1, 1)$hessian, 2)


CthisNf$run(1, .2, 1)$hessian

[,1]       [,2]       [,3]
[1,]  0.3678794  0.3678794 -0.3678794
[2,]  0.1555337 -0.6051374  0.3678794
[3,] -0.3678794  0.3678794 -0.3678794

[,1]      [,2] [,3]
[1,]    0  1.000000   -1
[2,]    1 -1.644934    1
[3,]   -1  1.000000   -1



for(x in 1:50){
 if(abs(CthisNf$run(1, x, 1)$value - RthisNf$run(1, x, 1)$value) > .05) browser()
}

x <- .1
digamma(2) - digamma(1) + log(.1)

