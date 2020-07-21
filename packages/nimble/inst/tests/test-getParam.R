source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

context("Testing of getParam")


test_getParam(quote(dbern(prob = 0.2)))
test_getParam(quote(dbin(prob = 0.2, size = 3)))
test_getParam(quote(dbinom(prob = 0.2, size = 3)), dist = 'dbin')  
test_getParam(quote(dnegbin(prob = 0.2, size = 3)))
test_getParam(quote(dnbinom(prob = 0.2, size = 3)), dist = 'dnegbin')
test_getParam(quote(dpois(lambda = 2.5)))
test_getParam(quote(dbeta(shape1 = 1.5, shape2 = 2.5)))
test_getParam(quote(dbeta(mean = .6, sd = .05)))
test_getParam(quote(dchisq(df = 3)))
test_getParam(quote(dexp(rate = .3)))
test_getParam(quote(dexp(scale = 3)))
test_getParam(quote(dgamma(shape = 2, scale = 1.5)))
test_getParam(quote(dgamma(shape = 2, rate = 3.0)))
test_getParam(quote(dgamma(mean = 2.0, sd = 1.5)))
test_getParam(quote(dlnorm(meanlog = 2.0, taulog = 1.5)))
test_getParam(quote(dlnorm(meanlog = 2.0, sdlog = .8)))
test_getParam(quote(dlnorm(meanlog = 2.0, varlog = .6)))
test_getParam(quote(dlogis(location = 1.5, rate = .2)))
test_getParam(quote(dlogis(location = 1.5, scale = 5.0)))
test_getParam(quote(dnorm(mean = 10.5, sd = 1.5)))
test_getParam(quote(dnorm(mean = 10.5, var = 1.5)))
test_getParam(quote(dnorm(mean = 10.5, tau = 1.5)))
test_getParam(quote(dt(df = 3, mu = 1.5, tau = 0.9)))
test_getParam(quote(dt(df = 3, mu = 1.5, sigma2 = 1.1)))
test_getParam(quote(dt(df = 3, mu = 1.5, sigma = 1.2)))
test_getParam(quote(dunif(min = 1.2, max = 1.3)))
test_getParam(quote(dweib(shape = 1.2, scale = 1.3)))
test_getParam(quote(dweib(shape = 1.2, rate = 1.3)))
test_getParam(quote(dweib(shape = 1.2, lambda = 1.3)))

## We haven't written an extensive version of testing getParam for non-scalar parameters
## However the following covers testing that the size processing and eigenization steps work with getParam.

testCode <- nimbleCode({
    for(i in 1:3) x[i] ~ dnorm(0, 1)
    y[1:3] ~ dmnorm(x[1:3], mycov[1:3, 1:3])
})

y <- rnorm(3)
mycov <- diag(3)
testModel <- nimbleModel(testCode, data = list(y = y, mycov = diag(3)))
x <- rnorm(3)
testModel$x <- x

nf <- nimbleFunction(
    setup = function(model, mvNode){},
    run = function() {
        ans <- model$getParam(mvNode, 'mean')
        return(ans)
        returnType(double(1))
    },
    methods = list(
        test2 = function() {
            ans <- 1.1 + model$getParam(mvNode, 'mean')
            return(ans)
            returnType(double(1))
        },
        test3 = function(z = double(1)) {
            ans <- z + model$getParam(mvNode, 'mean')
            return(ans)
            returnType(double(1))
        })
)

nf1 <- nf(testModel, 'y[1:3]')
test_that('multivar 1', expect_equivalent(nf1$run(), testModel$x))
test_that('multivar 2', expect_equivalent(nf1$test2(), testModel$x + 1.1))
test_that('multivar 3', expect_equivalent(nf1$test3(11:13), testModel$x + 11:13))

Ctest <- compileNimble(testModel, nf1)
test_that('multivar 4', expect_equivalent(Ctest$nf1$run(), Ctest$testModel$x))
test_that('multivar 5', expect_equivalent(Ctest$nf1$test2(), Ctest$testModel$x + 1.1))
test_that('multivar 6', expect_equivalent(Ctest$nf1$test3(11:13), Ctest$testModel$x + 11:13))

# basic non-scalar test
code = nimbleCode({
    a[1:3] ~ dmnorm(mu[1:3],pr[1:3,1:3])
})
pr1 = diag(3)
pr1[1,2]=pr1[2,1]=.3
pr2 <- pr1
pr1[1,2]=pr1[2,1]=.5

m = nimbleModel(code, inits =list(mu=rep(1,3), pr = pr1))
cm = compileNimble(m)

cm$pr <- pr2
cm$calculate(cm$getDependencies('pr'))

test_that('non-scalar 1', expect_equal(pr1, m$getParam('a', 'prec')))
test_that('non-scalar 2', expect_equal(pr2, cm$getParam('a', 'prec')))
test_that('non-scalar 3', expect_equal(solve(pr1), m$getParam('a', 'cov')))
test_that('non-scalar 4', expect_equal(solve(pr2), cm$getParam('a', 'cov')))

test_that('getParam, three-dimensional', {
    dtest <- nimbleFunction(
        run = function(x = double(0), theta = double(3), log = integer(0, default = 0)) {
            returnType(double(0))
            return(0)
        })
    rtest <- nimbleFunction(
        run = function(n = integer(0), theta = double(3)) {
            returnType(double(0))
            return(0)
        })
    temporarilyAssignInGlobalEnv(dtest)
    temporarilyAssignInGlobalEnv(rtest)
    code <- nimbleCode({
        y ~ dtest(theta[1:2,1:3,1:4])
    })
    init <- array(as.numeric(1:24), c(2,3,4))
    m <- nimbleModel(code, inits = list(theta = init))
    cm <- compileNimble(m)
    expect_identical(m$getParam('y','theta'), NULL, 'getParam_3D in uncompiled model')
    expect_identical(m$getParam('y','theta'), NULL, 'getParam_3D in compiled model')
    
    mynf <- nimbleFunction(setup = function(model, node, param) {},
                           run = function() {
                               tmp <- model$getParam(node, param)
                           })
    rnf <- mynf(m, 'y', 'theta')
    expect_identical(rnf$run(), NULL, 'getParam_3D in uncompiled nf')
    expect_error(cnf <- compileNimble(rnf, project = m), 'Failed to create the shared library')
    deregisterDistributions('dtest')  ## so can use it again in next test
})

test_that("Testing invalid parameter name in getParam", {
    code <- nimbleCode({
        for(i in 1:3)
            y[i] ~ dnorm(0, 1)
    })
    m <- nimbleModel(code)
    expect_error(m$getParam('y[1]', 'mu'), "parameter 'mu' not found")
    mynf <- nimbleFunction(
        setup=function(model,nodes){},
        run=function(){
            a = 0
            for(i in seq_along(nodes))
                a <- a + model$getParam(nodes[i], 'mu')
            print(a)
    })
    rnf <- mynf(m, c('y[1]','y[2]','y[3]'))
    cm <- compileNimble(m)
    expect_error(cnf <- compileNimble(rnf, project = m), "parameter 'mu' not found")
})

test_that('getParam, user-defined integer-valued', {
    dtest <- nimbleFunction(
        run = function(x = integer(0), thetaInt = integer(0), thetaDbl = double(0), log = integer(0, default = 0)) {
            returnType(double(0))
            return(0)
        })
    rtest <- nimbleFunction(
        run = function(n = integer(0), thetaInt = integer(0), thetaDbl = double(0)) {
            returnType(double(0))
            return(0)
        })
    temporarilyAssignInGlobalEnv(dtest)
    temporarilyAssignInGlobalEnv(rtest)

    code <- nimbleCode({
        y ~ dtest(1, 1.3)
    })
    m <- nimbleModel(code, data = list(y = 0))
    cm <- compileNimble(m)
    expect_identical(cm$getParam('y','value'), m$getParam('y', 'value'), 'issue with getParam with value')
    expect_identical(cm$getParam('y','value'), 0, 'issue with getParam with value')
    expect_identical(cm$getParam('y','thetaInt'), m$getParam('y', 'thetaInt'), 'issue with getParam with integer parameter')
    expect_identical(cm$getParam('y','thetaInt'), 1, 'issue with getParam with integer parameter')
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
