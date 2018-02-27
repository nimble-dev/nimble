source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
nimbleOptions(experimentalEnableDerivs = TRUE)
nimbleOptions(showCompilerOutput = TRUE)
context("Testing of derivatives for distributions and dsl functions")

  testNFgen <- nimbleFunction(setup = function(){
    cholMat = chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2))
  },
  run = function(x = double(2), df = double()){
    returnType(ADNimbleList())
    return(nimDerivs(ADmethod(x, cholMat, df), order = c(0,1,2)))
  },
  methods = list(
    ADmethod = function(x = double(2, c(2,2)), chol = double(2, c(2,2)), df = double()){
      returnType(double(0))
      return(dwish_chol(x, cholesky = chol, df = df, scale_param = TRUE))
    }
  ),
  enableDerivs = 'ADmethod')
  nimbleOptions(pauseAfterWritingFiles = TRUE)
  
  file.copy("C:\\Users\\iateb\\Documents\\GitHub\\nimble\\packages\\nimble\\inst\\include\\TMB\\distributions_R.hpp",
  "C:\\Users\\iateb\\Documents\\R\\win-library\\3.4\\nimble\\include\\TMB\\distributions_R.hpp",
  overwrite = TRUE)
  
  
  
  testNF <- testNFgen()
  ctestNF <- compileNimble(testNF)
nimbleOptions(pauseAfterWritingFiles = FALSE)

# testNF$ADmethod( matrix(c(1,.7,.7,1), nrow = 2),chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2)), 3)
ctestNF$ADmethod( matrix(c(1,.7,.7,1), nrow = 2),chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2)), 3)
ctestNF$run( matrix(c(1,.7,.7,1), nrow = 2), 3)


dwish_chol(matrix(c(1,.7,.7,1), nrow = 2),chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2)), 3)

# testNF$run( matrix(c(1,.7,.7,1), nrow = 2))

chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2))
chol(matrix(c(1,.7,.7,1), nrow = 2))
ADindependentVars[0] = 1;
ADindependentVars[1] = 0;
ADindependentVars[2] = 0;
ADindependentVars[3] = 1;
ADindependentVars[4] = 1;
ADindependentVars[5] = 0;
ADindependentVars[6] = 0;
ADindependentVars[7] = 1;
ADindependentVars[8] = 3;
ADindependentVars[9] = 1;

# untested:
# 'ddirch',
# 'dmvt', 
#  'dwish'

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
# 'dmulti', # not implemented
# 'dnegbin'
# 'dnorm',
# 'dpois',
# 'dt',
# 'dt_nonstandard',
# 'dunif',
# 'dweib'

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
    list(x = -1, df = 1),
    # list(x = 0, df = 1),  #R derivs here are incorrect
    list(x = 1, df = -1),
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
    list(x = c(0,0,0), alpha = c(0,0,0)),
    list(x = c(0,-1,0), alpha = c(0,0,0)),
    list(x = c(0,0,0), alpha = c(-1,0,0)),
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

distributionArgsList[['dmnorm_chol']] <- list(
  distnName = 'dmnorm_chol',
  args = list(x = quote(double(1, 2)),
              mean = quote(double(1, 2)),
              cholesky = quote(double(2, c(2, 2)))),
  argsValues = list(
    list(x = numeric(2), mean = numeric(2), cholesky = -diag(2)),
    # list(x = numeric(2), mean = numeric(2), cholesky = matrix(c(0,0,0,0), nrow = 2)) R and C inconsistency
    list(x = numeric(2), mean = numeric(2), cholesky = diag(2))
    # list(x = c(1.3, 4.1), mean = c(1,4), cholesky = chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2))) # inconsistency between R and C++ hessians, hypothesis is that R finite element diff. is giving incorrect results for cholesky param.
  )
)

#   distnName = 'dmulti',
#   args = list(x = quote(double(1, 3)),
#               size = quote(double(0)),
#               prob = quote(double(1, 3))),
#   argsValues = list(
#     list(x = numeric(3), size = 0, prob = numeric(3)),
#     list(x = numeric(3) - 1, size = -3, prob = numeric(3)),
#     list(x = numeric(3) - 1, size = -3, prob = numeric(3) + 1/3),
#     list(x = 1:3, size = sum(1:3), prob = numeric(3) + 1/3)
#   ),
#   WRT = c('prob')
# )

distributionArgsList[['dnbinom']] <- list(
  distnName = 'dnbinom',
  args = list(x = quote(double(0)),
              prob = quote(double(0)),
              size = quote(double(0))),
  argsValues = list(
    list(x = -1, prob = .1, size = 1),
    list(x = 0, prob = .1, size = 1),
    list(x = 1.5, prob = .1, size = 1),
    list(x = 1.5, prob = .1, size = 2),
    list(x = 2, prob = .1, size = 1),
    list(x = 2, prob = .1, size = 2.5),
    list(x = 2, prob = -1, size = 3),
    list(x = 2, prob = 10, size = 3),
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
              tau = quote(double(0))),
  argsValues = list(
    # list(x = -1, mean = 1, tau = 0), # same issue with inconsistent R/C++ derivs of un-logged prob wrt tau
    list(x = -1, mean = -1, tau = -1),
    list(x = -1, mean = 1, tau = 1),
    list(x = .1, mean = 1, tau = 2),
    list(x = 22.2, mean = 10, tau = 15.2))
)

distributionArgsList[['dpois']] <- list(
  distnName = 'dpois',
  args = list(x = quote(double(0)),
              lambda = quote(double(0))),
  argsValues = list(
    list(x = -1, lambda = -1),
    list(x = 0, lambda = -1),
    list(x = 1, lambda = -1),
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
    list(x = -1, df = -1),
    list(x = -1, df = 0),
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
    list(x = -1, df = -1, mu = 0, sigma = 1),
    list(x = -1, df = 0, mu = 0, sigma = 1),
    #list(x = -1, df = 0, mu = -1, sigma = 0), # R and C derivs don't align
    #list(x = -1, df = 1, mu = -1, sigma = 0), # R and C derivs don't align
    # list(x = -1, df = 2, mu = -2, sigma = 0), # R and C derivs don't align
    list(x = -1, df = 2, mu = -2, sigma = -1),
    list(x = .1, df = .5, mu = 0, sigma = 1),
    list(x = 22.2, df = 10, mu = 0, sigma = 1))
)


distributionArgsList[['dunif']] <- list(
  distnName = 'dunif',
  args = list(x = quote(double(0)),
              min = quote(double(0)),
              max = quote(double(0))),
  argsValues = list(
    list(x = -1, min = 0, max = 1),
    list(x = -1, min = 0, max = 0),
    list(x = 0, min = 0, max = 0),
    list(x = 0, min = 0, max = -1),
    # list(x = 0, min = 0, max = 1), # R derivs problem at boundary
    list(x = .5, min = 0, max = 1))
)

distributionArgsList[['dweibull']] <- list(
  distnName = 'dweibull',
  args = list(x = quote(double(0)),
              shape = quote(double(0)),
              scale = quote(double(0))),
  argsValues = list(
    list(x = -1, shape = 0, scale = 1),
    list(x = -1, shape = 0, scale = 0),
    # list(x = 0, shape = 0, scale = 0), #R derivs at x = 0 incorrect
    # list(x = 0, shape = 1, scale = 1),
    list(x = 0.5, shape = 0, scale = 0),
    list(x = 0.5, shape = 1, scale = 1),
    list(x = 2, shape = 3, scale = 2))
)


distributionArgsList[['dmvt_chol']] <- list(
  distnName = 'dmvt_chol',
  args = list(x = quote(double(1, 2)),
              mu = quote(double(1, 2)),
              cholesky = quote(double(2, c(2, 2))),
              df = quote(double()),
              prec_param = quote(double())),
  argsValues = list(
    # list(x = numeric(2), mu = numeric(2), cholesky = -diag(2), df = 1, prec_param = 1),
    # list(x = numeric(2), mean = numeric(2), cholesky = matrix(c(0,0,0,0), nrow = 2)) R and C inconsistency
    # list(x = numeric(2), mu = numeric(2), cholesky = diag(2), df = 1, prec_param = 1),
    list(x = c(1.3, 4.1), mu = c(1,4),
         cholesky = chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2)),
         df = 3, prec_param = 0),
    list(x = c(12.1, 42.1), mu = c(10,40),
         cholesky = chol(matrix(c(13.2, 2.14, 2.14, 2.73), nrow = 2)),
         df = 3, prec_param = 0),
    list(x = c(1.3, 4.1), mu = c(1,4),
         cholesky = chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2)),
         df = 3, prec_param = 1),
    list(x = c(12.1, 42.1), mu = c(10,40),
         cholesky = chol(matrix(c(13.2, 2.14, 2.14, 2.73), nrow = 2)),
         df = 3, prec_param = 1)
      )
)

distributionArgsList[['dmvt_chol']] <- list(
  distnName = 'dmvt_chol',
  args = list(x = quote(double(1, 2)),
              mu = quote(double(1, 2)),
              cholesky = quote(double(2, c(2, 2))),
              df = quote(double()),
              prec_param = quote(double())),
  argsValues = list(
    # list(x = numeric(2), mu = numeric(2), cholesky = -diag(2), df = 1, prec_param = 1),
    # list(x = numeric(2), mean = numeric(2), cholesky = matrix(c(0,0,0,0), nrow = 2)) R and C inconsistency
    # list(x = numeric(2), mu = numeric(2), cholesky = diag(2), df = 1, prec_param = 1),
    list(x = c(1.3, 4.1), mu = c(1,4),
         cholesky = chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2)),
         df = 3, prec_param = 0),
    list(x = c(12.1, 42.1), mu = c(10,40),
         cholesky = chol(matrix(c(13.2, 2.14, 2.14, 2.73), nrow = 2)),
         df = 3, prec_param = 0),
    list(x = c(1.3, 4.1), mu = c(1,4),
         cholesky = chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2)),
         df = 3, prec_param = 1),
    list(x = c(12.1, 42.1), mu = c(10,40),
         cholesky = chol(matrix(c(13.2, 2.14, 2.14, 2.73), nrow = 2)),
         df = 3, prec_param = 1)
  )
)

distributionArgsList[['dwish_chol']] <- list(
  distnName = 'dwish_chol',
  args = list(x = quote(double(2, c(2,2))),
              cholesky = quote(double(2, c(2, 2))),
              df = quote(double()),
              scale_param = quote(double())),
  argsValues = list(
    list(x = matrix(c(1,.7,.7,1), nrow = 2), 
         cholesky = chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2)),
         df = 3, scale_param = 1)
  )
)

nimbleOptions(pauseAfterWritingFiles = TRUE)

file.copy("C:\\Users\\iateb\\Documents\\GitHub\\nimble\\packages\\nimble\\inst\\include\\TMB\\distributions_R.hpp",
          "C:\\Users\\iateb\\Documents\\R\\win-library\\3.4\\nimble\\include\\TMB\\distributions_R.hpp",
          overwrite = TRUE)




  runFun <- gen_runFunCore(makeADDistributionTestList(distributionArgsList[['dwish_chol']]))
methodFun <- gen_runFunCore(makeADDistributionMethodTestList(distributionArgsList[['dwish_chol']]))
thisNf <- nimbleFunction(setup = function(){},
                      run = runFun,
                      methods = list(
                        method1 = methodFun
                      ),
                      enableDerivs = list('method1'))
testADDistribution(thisNf, distributionArgsList[['dwish_chol']]$argsValues,
                   distributionArgsList[['dwish_chol']]$distnName, debug = TRUE)

nimbleOptions(pauseAfterWritingFiles = FALSE)


ADindependentVars[0] = 1;
ADindependentVars[1] = 0;
ADindependentVars[2] = 0;
ADindependentVars[3] = 1;
ADindependentVars[4] = 1;
ADindependentVars[5] = 0;
ADindependentVars[6] = 0;
ADindependentVars[7] = 1;
ADindependentVars[8] = 3;
ADindependentVars[9] = 1;

testNFgen <- nimbleFunction(setup = function(){
                              cholMat = chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2))
                            },
                         run = function(x = double(2)){
                           nimDerivs(ADmethod(x, cholMat, 3))
                         },
                         methods = list(
                           ADmethod = function(x = double(2, c(2,2)), chol = double(2, c(2,2)), df = double()){
                             returnType(double(0))
                             return(dwish_chol(x, cholesky = chol, df = df, scale_param = TRUE))
                           }
                         ),
                         enableDerivs = 'ADmethod')


nimbleOptions(pauseAfterWritingFiles = TRUE)

testNF <- testNFgen()
ctestNF <- compileNimble(testNF)
nimbleOptions(pauseAfterWritingFiles = FALSE)

testNF$ADmethod( matrix(c(1,.7,.7,1), nrow = 2),chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2)), 3)
ctestNF$ADmethod( matrix(c(1,.7,.7,1), nrow = 2),chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2)), 3)

dwish_chol(x = matrix(c(1,.7,.7,1), nrow = 2), cholesky =  chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2)),
           df = 3, log = T, scale_param = T)

a <- rwish_chol( cholesky =  chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2)),
            df = 3)

x = matrix(c(12.1, 42.1), nrow = 2)
mu = matrix(c(10,40), nrow = 2)
cholesky = chol(matrix(c(13.2, 2.14, 2.14, 2.73), nrow = 2))
df = 3
n = 2
dmvt_chol(x, mu = mu, cholesky = cholesky, df = 3, prec_param = FALSE, TRUE)


# lapply(distributionArgsList, function(x){
#   runFun <- gen_runFunCore(makeADDistributionTestList(x))
#   methodFun <- gen_runFunCore(makeADDistributionMethodTestList(x))
#   thisNf <- nimbleFunction(setup = function(){},
#                            run = runFun,
#                            methods = list(
#                              method1 = methodFun
#                            ),
#                            enableDerivs = list('method1'))
#   testADDistribution(thisNf, x$argsValues,
#                      x$distnName)
#   
# })

x <- 0
shape <- .1
scale <- 1

(shape/scale)*(x/scale)^(shape-1)*exp(-(x/scale)^shape)
dweibull(0, .1, 1)

(gamma(k+r)/(factorial(k)*gamma(r)))*p^(r)*(1-p)^k
dnbinom(k,r,p)

testFn <- function(x, mu, chol){
  out <- dmvnc(x, mu, chol, TRUE)
  return(out)
}
nimDerivs(testFn(x = c(1.3, 4.1), mu = c(1,4), 
                 chol = chol(matrix(c(1.2, .14, .14, 2.7), nrow = 2))),wrt = c('chol'))


logres = 0
logres = logres + lgamma((df + n) / 2) - lgamma(df / (2)) - n * (log(sqrt(pi))) - n * log(df) / (2);
i = 1
logCholSum = 0
while(i <= n){
  logCholSum = logCholSum + log(cholesky[i, i]);
  i = i + 1
}
logres = logres -logCholSum;
eigenXcopy = x - mu;
eigenXcopy = solve(t(cholesky), eigenXcopy);
tmp = 0
for(i in 1:n)
  tmp = tmp + eigenXcopy[i] * eigenXcopy[i];

logres = logres -0.5 * (df + n) * log(1 + tmp / df);
logres




testMod1 <- nimbleModel(testMod)
ctestMod1 <- compileNimble(testMod1)

as.call(list(quote(dbinom), x = 1, size = 2, prob = .1))
test <- list(a = 1, b =2)
lapply(test, function(x){return(x[[1]])})
