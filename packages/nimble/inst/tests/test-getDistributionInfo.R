source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

## Testing functions that query distribution info based on distribution name or node/variable names

context('Testing distributions API')

# tests of distribution functions applied to distribution name

try(test_that("Test that requesting unknown distribution returns error",
              expect_error(getDistributionInfo('dfoobar'),
    info = "No error when 'dfoobar' provided to getDistributionInfo")))

try(test_that("Test that getDistributionInfo returns parameter dimension",
                 expect_identical(getDistributionInfo('dnorm')$types$tau$nDim, 0,
                 info = "getDistributionInfo doesn't correctly give dimension of tau in dnorm")))

try(test_that("Test that requesting unknown distribution returns error",
                               expect_error(isDiscrete('dfoobar'),
                                                   info = "No error when 'dfoobar' provided to isDiscrete")))

try(test_that("Test of isDiscrete",
                                  expect_identical(isDiscrete('dnorm'), FALSE,
                                                                       info = "isDiscrete says dnorm is discrete")))

try(test_that("Test of isDiscrete",
                                                   expect_identical(isDiscrete('dbin'), TRUE,
                                                                                                                                              info = "isDiscrete says dbin is not discrete")))

try(test_that("Test of isUserDefined",
                                                expect_identical(isUserDefined('dnorm'), FALSE,
                                               info = "isUserDefined says dnorm is user-defined")))

dmyexp <- nimbleFunction(
        run = function(x = double(0), rate = double(0, default = 1),
                               log = integer(0, default = 0)) {
                    returnType(double(0))
                            logProb <- log(rate) - x*rate
                            if(log) return(logProb)
                    else return(exp(logProb))
        })
assign('dmyexp', dmyexp, envir = .GlobalEnv)

rmyexp <- nimbleFunction(
        run = function(n = integer(0), rate = double(0, default = 1)) {
                    returnType(double(0))
                            if(n != 1) print("rmyexp only allows n = 1; using n = 1.")
                    dev <- runif(1)
                            return(-log(1-dev) / rate)
        })
assign('rmyexp', rmyexp, envir = .GlobalEnv)

try(test_that("Test of registerDistributions with character string",
              expect_output(registerDistributions('dmyexp'),
                            "Registering the following user-provided distributions: dmyexp.")))

try(test_that("Test of registerDistributions with character string for missing distribution",
              expect_error(registerDistributions('dfoo'))))

try(test_that("Test of isUserDefined",
                                                                 expect_identical(isUserDefined('dmyexp'), TRUE,
                                                                                                                                    info = "isUserDefined says dmyexp is not user-defined")))

try(test_that("Test of pqDefined",
                                 expect_identical(pqDefined('dmyexp'), FALSE,
                                 info = "pqDefined says pqDefined is TRUE for dmyexp")))

try(test_that("Test of pqDefined",
                                                  expect_identical(pqDefined('dnorm'), TRUE,
                                                       info = "isUserDefined says pqDefined is FALSE for dnorm")))


output <- 0; names(output) <- 'value'

try(test_that("Test of getDimension, value only",
              expect_identical(getDimension('dnorm'), output,
                               info = "incorrect value dimension for dnorm")))

output <- rep(0, 5); names(output) <- c('value', 'mean','sd', 'tau', 'var')
try(test_that("Test of getDimension, value and params",
              expect_equal(getDimension('dnorm', includeParams = TRUE), output,
                               info = "incorrect dimensions for dnorm")))

output <- c(1,0); names(output) <- c('prob', 'size')

try(test_that("Test of getDimension, params only",
              expect_equal(getDimension('dmulti', params = c('prob', 'size')), output,
                               info = "incorrect param dimensions for dmulti")))

try(test_that("Test of getDimension error checking",
              expect_error(getDimension('dmulti', params = c('foo')), 
                               info = "getDimension not detecting erroneous parameter names")))

    
try(test_that("Test of getParamNames",
              expect_identical(getParamNames('dnorm'), c('value','mean','sd','tau','var'),
                               info = "incorrect param names from getParamNames for dnorm")))

try(test_that("Test of getParamNames",
           expect_identical(getParamNames('dnorm', includeValue = FALSE), c('mean','sd','tau','var'),
                                       info = "incorrect param names without value from getParamNames for dnorm")))

# test of distribution functions applied to model variables/nodes

code <- nimbleCode({
    for(i in 1:10)
        y[i] ~ dnorm(mu, sd = sigma)
    mu ~ dnorm(0, 1)
    sigma ~ T(dgamma(a, 1), 0, 5)
    a ~ dgamma(1, 1)
    z <- mu + 3
    u ~ dmyexp(a)
    w <- u + 3
    x[1:2] ~ dmnorm(zeroes[1:2], prec[1:2, 1:2])
    yy ~ dpois(a)
    uu ~ dbinom(p, 1)
    p ~ dunif(0,1)
})

m <- nimbleModel(code, data = list(y = rnorm(10)),
                 constants = list(zeroes = rep(0,2), prec = diag(rep(1,2))), calculate=FALSE)

code2 <- nimbleCode({
    for(i in 1:10)
        y[i] ~ dnorm(mu, sd = sigma)
    mu ~ dnorm(0, 1)
    sigma ~ T(dgamma(a, 1), 0, 5)
    a ~ dgamma(1, 1)
    z <- mu + 3
    u ~ dmyexp(a)
    w <- u + 3
    x[1:2] ~ dmnorm(zeroes[1:2], prec[1:2, 1:2])
    yy ~ dpois(a)
})

m2 <- nimbleModel(code2, data = list(y = rnorm(10)),
                 constants = list(zeroes = rep(0,2), prec = diag(rep(1,2))), calculate=FALSE)

print('A gets here ok')

out <- c(rep(TRUE, 10), FALSE, TRUE)
names(out) <- m$expandNodeNames(c('y', 'mu', 'w'))

print('B.1 gets here ok')

try(test_that("Test of isEndNode",
              expect_identical(m$isEndNode(c('y', 'mu', 'w')), out,
                               info = "incorrect results from isEndNode")))

print('B.2 gets here ok')  ## this one ok


##try(
##    test_that("Test of isEndNode, unknown node",
##              expect_error(
##                  m$isEndNode(c('zzz', 'mu', 'w')),
##                  ##m$isEndNode('w')
##                  ##m$isEndNode('mu')
##                  ##m$isEndNode('zzz')
##                  out,
##                  info = "unknown node not detected in isEndNode"
##              )
##              )
##)
## 
## 
## 
## 
## 
##try(test_that("Test of isEndNode, unknown node",
##              expect_error(m$isEndNode(c('zzz', 'mu', 'w')), out,
##                               info = "unknown node not detected in isEndNode")))

print('B.3 gets here ok')  ## NEVER GETS HERE

vars <- c('y', 'mu', 'w', 'x')
out <- c(rep('dnorm', 11), NA, 'dmnorm')
names(out) <- m$expandNodeNames(vars)

print('B.4 gets here ok')

try(test_that("Test of getDistribution",
              expect_identical(m$getDistribution(vars), out,
                               info = "incorrect results from getDistribution")))

print('B.5 gets here ok')

vars <- c('y', 'yy', 'w', 'x')
out <- c(rep(FALSE, 10), TRUE, NA, FALSE)
names(out) <- m$expandNodeNames(vars)

print('C gets here ok')

try(test_that("Test of isDiscrete model method",
              expect_identical(m$isDiscrete(vars), out,
                               info = "incorrect results from isDiscrete")))

vars <- c('y', 'yy', 'w', 'x', 'uu')
out <- c(rep(FALSE, 10), FALSE, NA, FALSE, TRUE)
names(out) <- m$expandNodeNames(vars)

try(test_that("Test of isBinary model method",
              expect_identical(m$isBinary(vars), out,
                               info = "incorrect results from isBinary")))

vars <- c('y', 'yy', 'w', 'x')
out <- c(rep(FALSE, 10), FALSE, NA, FALSE)
names(out) <- m$expandNodeNames(vars)

try(test_that("Test of isBinary model method",
              expect_identical(m2$isBinary(vars), out,
                               info = "incorrect results from isBinary")))


print('D gets here ok')

vars <- c('y', 'yy', 'w', 'x', 'uu', 'z')
out <- !c(rep(FALSE, 10), FALSE, TRUE, FALSE, FALSE, TRUE)
names(out) <- m$expandNodeNames(vars)

try(test_that("Test of isStoch model method",
              expect_identical(m$isStoch(vars), out,
                               info = "incorrect results from isStoch")))

vars <- c('y', 'yy', 'w', 'x', 'uu', 'z')
out <- c(rep(FALSE, 10), FALSE, TRUE, FALSE, FALSE, TRUE)
names(out) <- m$expandNodeNames(vars)

try(test_that("Test of isDeterm model method",
              expect_identical(m$isDeterm(vars), out,
                               info = "incorrect results from isDeterm")))

vars <- c('y', 'sigma', 'w', 'x', 'uu', 'z')
out <- c(rep(FALSE, 10), TRUE, FALSE, FALSE, FALSE, FALSE)
names(out) <- m$expandNodeNames(vars)

try(test_that("Test of isTruncated model method",
              expect_identical(m$isTruncated(vars), out,
                               info = "incorrect results from isTruncated")))

vars <- c('y', 'sigma', 'w', 'x', 'uu', 'z')
out <- c(rep(FALSE, 10), FALSE, NA, TRUE, FALSE, NA)
names(out) <- m$expandNodeNames(vars)

print('E gets here ok')

try(test_that("Test of isUnivariate model method",
              expect_identical(m$isUnivariate(vars), out,
                               info = "incorrect results from isUnivariate")))
out <- 0
names(out) <- "value"
try(test_that("Test of getDimension, value only",
              expect_identical(m$getDimension('mu'), out,
                               info = "incorrect result from getDimension for value only, scalar")))
out <- 1
names(out) <- "value"
try(test_that("Test of getDimension, value only",
              expect_identical(m$getDimension('x'), out,
                               info = "incorrect result from getDimension for value only, vector")))

out <- c(1,1,2,0,2,2)
names(out) <- c("value","mean","cholesky","prec_param","prec","cov")
try(test_that("Test of getDimension, value and params",
              expect_identical(m$getDimension('x', includeParams = TRUE), out,
                               info = "incorrect result from getDimension for value and params, vector")))

print('F.1 gets here ok')

out <- c(2,1)
names(out) <- c('cov','mean')
try(test_that("Test of getDimension, specific params",
              expect_identical(m$getDimension('x', params = c('cov','mean')), out,
                               info = "incorrect result from getDimension for specific params, vector")))

print('F.2 gets here ok') ## reaches here

## this next test fails:


##m$getDimension('foo', includeParams = TRUE)
## 
##try(
## 
##    test_that("Test of getDimension, value only",
##              expect_error(m$getDimension('foo', includeParams = TRUE),
##                           info = "not detecting missing node in getDimension"))
## 
##)
## 
## 
##try(test_that("Test of getDimension, value only",
##              expect_error(m$getDimension('foo', includeParams = TRUE),
##                               info = "not detecting missing node in getDimension")))

print('F.3 gets here ok')  ## NOT HERE!!!

try(test_that("Test of getDimension, value only",
              expect_error(m$getDimension('x', params = 'foo'),
                               info = "not detecting missing param in getDimension")))

# test of user-defined distributions

print('G gets here ok')

ddirchmulti <- nimbleFunction(
            run = function(x = double(1), alpha = double(1), size = double(0),
                log = integer(0, default = 0)) {

                returnType(double(0))
                logProb <- lgamma(size) - sum(lgamma(x)) + 
                        lgamma(sum(alpha)) -
                        sum(lgamma(alpha)) + sum(lgamma(alpha + x)) - 
                        lgamma(sum(alpha) + size)
                if(log) return(logProb)
                        else return(exp(logProb))
})
assign('ddirchmulti', ddirchmulti, envir = .GlobalEnv)

print('H gets here ok')

rdirchmulti <- nimbleFunction(
            run = function(n = integer(0), alpha = double(1), 
                size = double(0)) {

                returnType(double(1))
                if(n != 1) print("rdirchmulti only allows n = 1; using n = 1.")
                p <- rdirch(1, alpha)
                return(rmulti(1, size = size, prob = p))
})
assign('rdirchmulti', rdirchmulti, envir = .GlobalEnv)

print('I gets here ok')

code <- nimbleCode({
    for(i in 1:n)
                                        # likelihood 
        y[i,1:K] ~ ddirchmulti(alpha[topic[i], 1:K], N[i])
                                        # priors for hyperparameters
    for(tp in 1:M)
        for(k in 1:K)
            alpha[tp, k] ~ dunif(0, 100)
})
const <- list(M = 2, K = 4, n = 5, N = rep(1000, 5),
      topic = c(1, 1, 1, 2, 2))
alphaInits <- rbind(c(10, 30, 100, 3), c(12, 15, 15, 8))

print('J gets here ok')

m <- nimbleModel(code, constants = const,
                  inits = list(alpha = alphaInits))

# won't be correct because haven't used registerDistributions such that it's noted as discrete
## try(test_that("Test of isDiscrete for user-defined ddirchmulti",
##                                                    expect_identical(isDiscrete('ddirchmulti'), TRUE,
##                                                                                                                                               info = "isDiscrete says ddirchmulti is not discrete")))


print('K gets here ok')



output <- c(1,1,0); names(output) <- c('value', 'alpha', 'size')

try(test_that("Test of getDimension, value and params, for user-defined ddirchmulti",
              expect_equal(getDimension('ddirchmulti', includeParams = TRUE), output,
                               info = "incorrect param dimensions for ddirchmulti")))


print('L gets here ok')

set.seed(0)
try(test_that("Test of use of ddirchmulti in R model",
              expect_silent(m$simulate('y'))))
try(test_that("Test of use of ddirchmulti in R model",
              expect_equal(m$y[5,4], 149, info = "unexpected result from use of ddirchmulti")))

print('M gets here ok')


set.seed(0)
cm <- compileNimble(m)

print('N gets here ok')


try(test_that("Test of use of ddirchmulti in R model",
              expect_silent(cm$simulate('y'))))
try(test_that("Test of use of user-defined ddirchmulti",
              expect_identical(m$y, cm$y, info = "R and compiled model values from ddirchmulti not the same")))

print('O gets here ok')


