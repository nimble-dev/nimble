source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

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
        run = function(x = double(0), thetaInt = integer(0), thetaDbl = double(0), log = integer(0, default = 0)) {
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

test_that("getParam works in various multi-instance combinations", {
  code <- nimbleCode({
    a ~ dnorm(0, 1)
    for(i in 1:3) x[i] ~ dnorm(a, 1)
    b ~ dbeta(a,3)
    for(i in 1:3) y[i] ~ dbeta(5, a)
  })
  
  m <- nimbleModel(code, inits = list(a = 1, x = 2:4, b = .5, y = c(0.6, 0.7, 0.8)))
  cm <- compileNimble(m)
  
  # two singletons
  getSomeParams1 <- nimbleFunction(
    setup = function(model, nodes, param = "value") {
      nodes <- model$expandNodeNames(nodes)
      numNodes <- length(nodes)    
    },
    run = function() {
      v <- numeric(numNodes)
      for(i in 1:numNodes) {
        v[i] <- model$getParam(nodes[i], param)
      }
      return(v);
      returnType(double(1))
    }
  )
  gsp1_1 <- getSomeParams1(m, "b", "mean")
  gsp1_2 <- getSomeParams1(m, "x[1]", "mean")
  cgsp1 <- compileNimble(gsp1_1, gsp1_2, project = m)
  expect_identical(gsp1_1$run(), cgsp1$gsp1_1$run())
  expect_identical(gsp1_2$run(), cgsp1$gsp1_2$run())

  # one singleton and one doublet, each from a single declaration
  getSomeParams2 <- nimbleFunction(
    setup = function(model, nodes, param = "value") {
      nodes <- model$expandNodeNames(nodes)
      numNodes <- length(nodes)
    },
    run = function() {
      v <- numeric(numNodes)    
      for(i in 1:numNodes) {
        v[i] <- model$getParam(nodes[i], param)
      }
      return(v);
      returnType(double(1))
    }
  )
  gsp2_1 <- getSomeParams2(m, "b", "mean")
  gsp2_2 <- getSomeParams2(m, "x[2:3]", "mean")
  cgsp2 <- compileNimble(gsp2_1, gsp2_2, project = m)
  expect_identical(gsp2_1$run(), cgsp2$gsp2_1$run())
  expect_identical(gsp2_2$run(), cgsp2$gsp2_2$run())

  # two doublets, each from a single declaration
  getSomeParams3 <- nimbleFunction(
    setup = function(model, nodes, param = "value") {
      nodes <- model$expandNodeNames(nodes)
      numNodes <- length(nodes)
    },
    run = function() {
      v <- numeric(numNodes)    
      for(i in 1:numNodes) {
        v[i] <- model$getParam(nodes[i], param)
      }
      return(v);
      returnType(double(1))
    }
  )
  gsp3_1 <- getSomeParams3(m, "x[1:3]", "mean")
  gsp3_2 <- getSomeParams3(m, "y[1:3]", "mean")
  cgsp3 <- compileNimble(gsp3_1, gsp3_2, project = m)
  expect_identical(gsp3_1$run(), cgsp3$gsp3_1$run())
  expect_identical(gsp3_2$run(), cgsp3$gsp3_2$run())

  # two doublets, with both from a multiple declarations
  getSomeParams4 <- nimbleFunction(
    setup = function(model, nodes, param = "value") {
      nodes <- model$expandNodeNames(nodes)
      numNodes <- length(nodes)
    },
    run = function() {
      v <- numeric(numNodes)    
      for(i in 1:numNodes) {
        v[i] <- model$getParam(nodes[i], param)
      }
      return(v);
      returnType(double(1))
    }
  )
  gsp4_1 <- getSomeParams4(m, c("b", "x[1:3]"), "mean")
  gsp4_2 <- getSomeParams4(m, c("b", "y[1:3]"), "mean")
  cgsp4 <- compileNimble(gsp4_1, gsp4_2, project = m)
  expect_identical(gsp4_1$run(), cgsp4$gsp4_1$run())
  expect_identical(gsp4_2$run(), cgsp4$gsp4_2$run())

  # two doublets, with one from a multiple declarations and one from a single declaration
  getSomeParams5 <- nimbleFunction(
    setup = function(model, nodes, param = "value") {
      nodes <- model$expandNodeNames(nodes)
      numNodes <- length(nodes)
    },
    run = function() {
      v <- numeric(numNodes)    
      for(i in 1:numNodes) {
        v[i] <- model$getParam(nodes[i], param)
      }
      return(v);
      returnType(double(1))
    }
  )
  gsp5_1 <- getSomeParams5(m, c("x[1:3]"), "mean")
  gsp5_2 <- getSomeParams5(m, c("b", "y[1:3]"), "mean")
  cgsp5 <- compileNimble(gsp5_1, gsp5_2, project = m)
  expect_identical(gsp5_1$run(), cgsp5$gsp5_1$run()) ## Failure in 0.10-0
  expect_identical(gsp5_2$run(), cgsp5$gsp5_2$run())
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
