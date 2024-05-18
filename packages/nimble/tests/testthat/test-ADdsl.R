source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
source(system.file(file.path('tests', 'testthat', 'AD_test_utils.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
CCopt <- nimbleOptions("useClearCompiledInADTesting")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)
nimbleOptions(useClearCompiledInADTesting = FALSE)

# tested:
# 'dbeta' (although boundary at x=0 and x=1 not consistent between R and c++),
# 'dbinom' (works as long as wrt = prob and not x or size),
# 'dchisq' (inconsistencies at df = 0, incorrect R derivs at x = 0),
# 'dcat' #not implemented
# 'dexp',
# 'dgamma',
# 'dinvgamma',
# 'dlogis',
# 'dlnorm',
# 'dmnorm_chol',
# 'dmulti',
# 'dnegbin'
# 'dnorm',
# 'dpois',
# 'dt',
# 'dt_nonstandard',
# 'dunif',
# 'dweib''
# 'ddirch',
# 'dmvt', 
#  'dwish'

distributionArgsList <- list()

distributionArgsList[['dbeta']] <- list(
  distnName = 'dbeta',
  args = list(x = quote(double(0)),
              shape1 = quote(double(0)),
              shape2 = quote(double(0))),
  argsValues = list(# list(x = -1, shape1 = 1, shape2 = 1), ## x below support
                   # list(x = 0, shape1 = 1, shape2 = 1), ## x at boundary of support, currently not working
                    # list(x = 1, shape1 = 1, shape2 = 1), ## x at boundary of support, currently not working
                    # list(x = 2, shape1 = 1,
                    #      shape2 = 1),
                    #list(x = .3,  ## first param not valid
                    #      shape1 = -1,
                    #      shape2 = 1),
                    # list(x = .3,  ## second param not valid
                    #      shape1 = 1,
                    #      shape2 = -1),
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
    # list(x = -1, prob = .1, size = 1),
    # list(x = 1.5, prob = .1, size = 1),
    # list(x = 1.5, prob = .1, size = 2),
    # list(x = 2, prob = .1, size = 1),
    # list(x = 2, prob = .1, size = 2.5),
    # list(x = 2, prob = -1, size = 3),
    # list(x = 2, prob = 10, size = 3),
    # list(x = 2, prob = 0, size = 3),
    # list(x = 2, prob = 1, size = 3),
    list(x = 2, prob = .1, size = 3),
    list(x = 25, prob = .5, size = 50)),
  WRT = c('prob') ## only take derivs wrt 'prob' argument, otherwise warnings
                  ## from R, and inconsistencies between R and C++ will occur.
)

# distributionArgsList[['dcat']] <- list(
#   distnName = 'dcat',
#   args = list(x = quote(double(0)),
#               prob = quote(double(1, 3))),
#   argsValues = list(
#     list(x = 3, prob = 1:3),
#     list(x = 2, prob = 1:3*2)
#   ),
#   WRT = c('prob')
# )

distributionArgsList[['dchisq']] <- list(
  distnName = 'dchisq',
  args = list(x = quote(double(0)),
              df = quote(double(0))),
  argsValues = list(
    # list(x = -1, df = 1),
    # list(x = 0, df = 1),  #R derivs here are incorrect
    # list(x = 1, df = -1),
    # list(x = 1, df = 0), #slight inconsistency between R and C++
    list(x = 1, df = 0.5),
    list(x = 1, df = 2),
    list(x = 12.5, df =22)
  )
)

distributionArgsList[['ddirch']] <- list(
  distnName = 'ddirch',
  args = list(x = quote(double(1, 3)),
              alpha = quote(double(1, 3))),
  argsValues = list(
    # list(x = c(0,0,0), alpha = c(0,0,0)),
    # list(x = c(0,-1,0), alpha = c(0,0,0)),
    # list(x = c(0,0,0), alpha = c(-1,0,0)),
    # list(x = c(.1,0,0), alpha = c(1,2,3)), # R and C Hessians don't match
    list(x = c(1/3,1/3,1/3), alpha = c(1,2,3))
  ),
  WRT = 'alpha'
)

distributionArgsList[['dexp']] <- list(
  distnName = 'dexp',
  args = list(x = quote(double(0)),
              rate = quote(double(0))),
  argsValues = list(
    # list(x = -1, rate = 1),
    # list(x = -1, rate = 0), inconsistency between R and C++
    # list(x = -1, rate = -1),
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
    # list(x = -1, shape = 1, rate = 1),
    # list(x = -1, shape = -1, rate = 1),
    # list(x = -1, shape = 1, rate = -1),
    list(x = .1, shape = 1, rate = 1),
    list(x = 22.2, shape = 10, rate = 15.2))
)

distributionArgsList[['dinvgamma']] <- list(
  distnName = 'dinvgamma',
  args = list(x = quote(double(0)),
              shape = quote(double(0)),
              scale = quote(double(0))),
  argsValues = list(
    # list(x = -1, shape = 1, scale = 1),
    # list(x = -1, shape = -1, scale = 1),
    # list(x = -1, shape = 1, scale = -1),
    list(x = .1, shape = 1, scale = 1),
    list(x = 22.2, shape = 10, scale = 15.2))
)

distributionArgsList[['dlogis']] <- list(
  distnName = 'dlogis',
  args = list(x = quote(double(0)),
              location = quote(double(0)),
              scale = quote(double(0))),
  argsValues = list(
    # list(x = -1, location = 1, scale = 0),
    # list(x = -1, location = -1, scale = -1),
    # list(x = -1, location = 1, scale = 1),
    list(x = .1, location = 1, scale = 2),
    list(x = 22.2, location = 10, scale = 15.2))
)

distributionArgsList[['dlnorm']] <- list(
  distnName = 'dlnorm',
  args = list(x = quote(double(0)),
              meanlog = quote(double(0)),
              sdlog = quote(double(0))),
  argsValues = list(
    # list(x = -1, meanlog = 1, taulog = 0),
    # list(x = -1, meanlog = -1, taulog = -1),
    # list(x = -1, meanlog = 1, taulog = 1),
    list(x = .1, meanlog = 1, sdlog = 2),
    list(x = 22.2, meanlog = 10, sdlog = 15.2))
)

distributionArgsList[['dmnorm_chol']] <- list(
  distnName = 'dmnorm_chol',
  args = list(x = quote(double(1, 2)),
              mean = quote(double(1, 2)),
              cholesky = quote(double(2, c(2, 2))),
              prec_param = quote(double())),
  argsValues = list(
    # list(x = numeric(2), mean = numeric(2), cholesky = -diag(2), prec_param = F),
    # list(x = numeric(2), mean = numeric(2), cholesky = -diag(2), prec_param = T),
    
    # list(x = numeric(2), mean = numeric(2), cholesky = matrix(c(0,0,0,0), nrow = 2)) R and C inconsistency
    list(x = numeric(2), mean = numeric(2), cholesky = diag(2), prec_param = F),
    list(x = numeric(2), mean = numeric(2), cholesky = diag(2), prec_param = T)
    
    # list(x = c(1.3, 4.1), mean = c(1,4), cholesky = chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2))) # inconsistency between R and C++ hessians, hypothesis is that R finite element diff. is giving incorrect results for cholesky param.
  ),
  WRT = c('x', 'mean', 'cholesky')
)

distributionArgsList[['dmulti']] <- list(
   distnName = 'dmulti',
   args = list(x = quote(double(1, 3)),
               size = quote(double(0)),
               prob = quote(double(1, 3))),
   argsValues = list(
#     list(x = numeric(3), size = 0, prob = numeric(3)),
#     list(x = numeric(3) - 1, size = -3, prob = numeric(3)),
#     list(x = numeric(3) - 1, size = -3, prob = numeric(3) + 1/3),
     list(x = 1:3, size = sum(1:3), prob = numeric(3) + 1/3)
   ),
   WRT = c('prob')
)
 

distributionArgsList[['dnbinom']] <- list(
  distnName = 'dnbinom',
  args = list(x = quote(double(0)),
              prob = quote(double(0)),
              size = quote(double(0))),
  argsValues = list(
    # list(x = -1, prob = .1, size = 1),
    # list(x = 0, prob = .1, size = 1),
    # list(x = 1.5, prob = .1, size = 1),
    # list(x = 1.5, prob = .1, size = 2),
    # list(x = 2, prob = .1, size = 1),
    # list(x = 2, prob = .1, size = 2.5),
    # list(x = 2, prob = -1, size = 3),
    # list(x = 2, prob = 10, size = 3),
    # list(x = 2, prob = 0, size = 3), #again issue with prob at boundaries
    # list(x = 2, prob = 1, size = 3),
    list(x = 2, prob = .1, size = 3),
    list(x = 25, prob = .5, size = 50)),
  WRT = c('prob') ## only take derivs wrt 'prob' argument, otherwise warnings
  ## from R, and inconsistencies between R and C++ will occur.
)
distributionArgsList[['dnorm']] <- list(
  distnName = 'dnorm',
  args = list(x = quote(double(0)),
              mean = quote(double(0)),
              sd = quote(double(0))),
  argsValues = list(
    # list(x = -1, mean = 1, tau = 0), # same issue with inconsistent R/C++ derivs of un-logged prob wrt tau
    # list(x = -1, mean = -1, tau = -1),
    list(x = -1, mean = 1, sd = 1),
    list(x = .1, mean = 1, sd = 2),
    list(x = 22.2, mean = 10, sd = 15.2))
)

distributionArgsList[['dpois']] <- list(
  distnName = 'dpois',
  args = list(x = quote(double(0)),
              lambda = quote(double(0))),
  argsValues = list(
    # list(x = -1, lambda = -1),
    # list(x = 0, lambda = -1),
    # list(x = 1, lambda = -1),
    # list(x = -1, lambda = 0),
    # list(x = 0, lambda = 0), #some R and C inconsistencies here
    # list(x = 1, lambda = 0),
    list(x = 1, lambda = 1),
    list(x = 14, lambda = 12)),
  WRT = c('lambda')
)

distributionArgsList[['dt']] <- list(
  distnName = 'dt',
  args = list(x = quote(double(0)),
              df = quote(double(0))),
  argsValues = list(
    # list(x = -1, df = -1),
    # list(x = -1, df = 0),
    list(x = .1, df = .5),
    list(x = 22.2, df = 10))
)

distributionArgsList[['dt_nonstandard']] <- list(
  distnName = 'dt_nonstandard',
  args = list(x = quote(double(0)),
              df = quote(double(0)),
              mu = quote(double(0)),
              sigma = quote(double(0))),
  argsValues = list(
    # list(x = -1, df = -1, mu = 0, sigma = 1),
    # list(x = -1, df = 0, mu = 0, sigma = 1),
    #list(x = -1, df = 0, mu = -1, sigma = 0), # R and C derivs don't align
    #list(x = -1, df = 1, mu = -1, sigma = 0), # R and C derivs don't align
    # list(x = -1, df = 2, mu = -2, sigma = 0), # R and C derivs don't align
    # list(x = -1, df = 2, mu = -2, sigma = -1),
    list(x = .1, df = .5, mu = 0, sigma = 1),
    list(x = 22.2, df = 10, mu = 0, sigma = 1))
)


distributionArgsList[['dunif']] <- list(
  distnName = 'dunif',
  args = list(x = quote(double(0)),
              min = quote(double(0)),
              max = quote(double(0))),
  argsValues = list(
    # list(x = -1, min = 0, max = 1),
    # list(x = -1, min = 0, max = 0),
    # list(x = 0, min = 0, max = 0),
    # list(x = 0, min = 0, max = -1),
    # list(x = 0, min = 0, max = 1), # R derivs problem at boundary
    list(x = .5, min = 0, max = 1))
)

distributionArgsList[['dweibull']] <- list(
  distnName = 'dweibull',
  args = list(x = quote(double(0)),
              shape = quote(double(0)),
              scale = quote(double(0))),
  argsValues = list(
    # list(x = -1, shape = 0, scale = 1),
    # list(x = -1, shape = 0, scale = 0),
    # list(x = 0, shape = 0, scale = 0), #R derivs at x = 0 incorrect
    # list(x = 0, shape = 1, scale = 1),
    # list(x = 0.5, shape = 0, scale = 0),
    list(x = 0.5, shape = 1, scale = 1),
    list(x = 2, shape = 3, scale = 2))
)



# distributionArgsList[['dmvt_chol']] <- list(
#   distnName = 'dmvt_chol',
#   args = list(x = quote(double(1, 2)),
#               mu = quote(double(1, 2)),
#               cholesky = quote(double(2, c(2, 2))),
#               df = quote(double()),
#               prec_param = quote(double())),
#   argsValues = list(
#     # list(x = numeric(2), mu = numeric(2), cholesky = -diag(2), df = 1, prec_param = 1),
#     # list(x = numeric(2), mean = numeric(2), cholesky = matrix(c(0,0,0,0), nrow = 2)) R and C inconsistency
#     # list(x = numeric(2), mu = numeric(2), cholesky = diag(2), df = 1, prec_param = 1),
#     list(x = c(1.3, 4.1), mu = c(1,4),
#          cholesky = chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2)),
#          df = 3, prec_param = 0),
#     list(x = c(12.1, 42.1), mu = c(10,40),
#          cholesky = chol(matrix(c(13.2, 2.14, 2.14, 2.73), nrow = 2)),
#          df = 3, prec_param = 0),
#     list(x = c(1.3, 4.1), mu = c(1,4),
#          cholesky = chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2)),
#          df = 3, prec_param = 1),
#     list(x = c(12.1, 42.1), mu = c(10,40),
#          cholesky = chol(matrix(c(13.2, 2.14, 2.14, 2.73), nrow = 2)),
#          df = 3, prec_param = 1)
#   )
# )

# distributionArgsList[['dwish_chol']] <- list(
#   distnName = 'dwish_chol',
#   args = list(x = quote(double(2, c(2,2))),
#               cholesky = quote(double(2, c(2, 2))),
#               df = quote(double()),
#               scale_param = quote(double())),
#   argsValues = list(
#     list(x = matrix(c(1,.8,.8,1), nrow = 2), 
#          cholesky = chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2)),
#          df = 3, scale_param = 1)
#   )
# )

lapply(distributionArgsList, function(x){
    test_that(paste0('AD for distribution ', x$distnName),
    {
        runFun <- gen_runFunCore(makeADDistributionTestList(x))
        methodFun <- gen_runFunCore(makeADDistributionMethodTestList(x))
        thisNf <- nimbleFunction(setup = function(){},
                                 run = runFun,
                                 methods = list(
                                     method1 = methodFun
                                 ),
                                 buildDerivs = list(method1 = list())) 
        test_ADDistribution(thisNf, x$argsValues,
                           x$distnName)
    })
})

nimbleOptions(enableDerivs = EDopt)
nimbleOptions(buildModelDerivs = BMDopt)
nimbleOptions(useClearCompiledInADTesting = CCopt)

do_test <- function(x){
    test_that(paste0('AD for distribution ', x$distnName),
    {
        runFun <- gen_runFunCore(makeADDistributionTestList(x))
        methodFun <- gen_runFunCore(makeADDistributionMethodTestList(x))
        thisNf <- nimbleFunction(setup = function(){},
                                 run = runFun,
                                 methods = list(
                                     method1 = methodFun
                                 ),
                                 buildDerivs = list(method1 = list()))
        test_ADDistribution(thisNf, x$argsValues,
                           x$distnName)
    })
}
