# This file contains tests of various modes of creating and invoking AD.

source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)

RCrelTol0 <- 1e-12
RCrelTol1 <- 1e-6
RCrelTol2 <- 1e-4

test_that('AD works from character buildDerivs with setup = TRUE', {
  adfun <- nimbleFunction(
    setup = TRUE,
    run = function(x = double(1)) {
      return(sum(exp(x)))
      returnType(double())
    },
    methods = list(
      derivsRun = function(x = double(1)) {
        ans <- nimDerivs(run(x), wrt = 1:length(x), order = 0:2)
        return(ans)
        returnType(ADNimbleList())
      }),
    buildDerivs = 'run')
  temporarilyAssignInGlobalEnv(adfun)

  Radfun1 <- adfun()
  Cadfun1 <- compileNimble(Radfun1)
  xRec <- c(.3, -1.2, 2.1)
  xTest <- c(-.4, 1.9, -0.8)
  Rans <- Radfun1$derivsRun(xTest)
  Cadfun1$derivsRun(xRec)
  Cans <- Cadfun1$derivsRun(xTest)
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)
})

test_that('AD works from empty list buildDerivs with setup = TRUE', {
  adfun <- nimbleFunction(
    setup = TRUE,
    run = function(x = double(1)) {
      return(sum(exp(x)))
      returnType(double())
    },
    methods = list(
      derivsRun = function(x = double(1)) {
        ans <- nimDerivs(run(x), wrt = 1:length(x), order = 0:2)
        return(ans)
        returnType(ADNimbleList())
      }),
    buildDerivs = list(run = list())
  )
  temporarilyAssignInGlobalEnv(adfun)
  Radfun1 <- adfun()
  Cadfun1 <- compileNimble(Radfun1)
  xRec <- c(.3, -1.2, 2.1)
  xTest <- c(-.4, 1.9, -0.8)
  Rans <- Radfun1$derivsRun(xTest)
  Cadfun1$derivsRun(xRec)
  Cans <- Cadfun1$derivsRun(xTest)
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)
})

test_that('AD works from list buildDerivs with empty ignore entry with setup = TRUE', {
  adfun <- nimbleFunction(
    setup = TRUE,
    run = function(x = double(1)) {
      return(sum(exp(x)))
      returnType(double())
    },
    methods = list(
      derivsRun = function(x = double(1)) {
        ans <- nimDerivs(run(x), wrt = 1:length(x), order = 0:2)
        return(ans)
        returnType(ADNimbleList())
      }),
    buildDerivs = list(run = list(ignore = character()))
  )
  temporarilyAssignInGlobalEnv(adfun)
  Radfun1 <- adfun()
  Cadfun1 <- compileNimble(Radfun1)
  xRec <- c(.3, -1.2, 2.1)
  xTest <- c(-.4, 1.9, -0.8)
  Rans <- Radfun1$derivsRun(xTest)
  Cadfun1$derivsRun(xRec)
  Cans <- Cadfun1$derivsRun(xTest)
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)
})

test_that('AD works from list with an ignore entry buildDerivs with setup = TRUE', {
  adfun <- nimbleFunction(
    setup = TRUE,
    run = function(x = double(1)) {
      ans <- 0
      for(i in 1:length(x))
        ans <- ans + exp(x[i])
      return(ans)
      returnType(double())
    },
    methods = list(
      derivsRun = function(x = double(1)) {
        ans <- nimDerivs(run(x), wrt = 1:length(x), order = 0:2)
        return(ans)
        returnType(ADNimbleList())
      }),
    buildDerivs = list(run = list(ignore = 'i'))
  )
  temporarilyAssignInGlobalEnv(adfun)
  Radfun1 <- adfun()
  Cadfun1 <- compileNimble(Radfun1)
  xRec <- c(.3, -1.2, 2.1)
  xTest <- c(-.4, 1.9, -0.8)
  Rans <- Radfun1$derivsRun(xTest)
  Cadfun1$derivsRun(xRec)
  Cans <- Cadfun1$derivsRun(xTest)
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)
})

test_that('AD works with calls to other methods and RCfunctions', {
  adhelper <- nimbleFunction(
    run = function(x = double(1)) {
      ans <- sum(exp(x))
      return(ans)
      returnType(double(0))
    },
    buildDerivs = TRUE
  )
  temporarilyAssignInGlobalEnv(adhelper)
  adfun <- nimbleFunction(
    setup = TRUE,
    run = function(x = double(1)) {
      ans <- 2*mysqrt(adhelper(x))
      return(ans)
      returnType(double())
    },
    methods = list(
      mysqrt = function(x = double()) {
        return(sqrt(x))
        returnType(double())
      },
      derivsRun = function(x = double(1)) {
        ans <- nimDerivs(run(x), wrt = 1:length(x), order = 0:2)
        return(ans)
        returnType(ADNimbleList())
      }),
    buildDerivs = c('run', 'mysqrt')
  )
  temporarilyAssignInGlobalEnv(adfun)
  Radfun1 <- adfun()
  Cadfun1 <- compileNimble(Radfun1) # This fails if adhelper is included.
  xRec <- c(.3, -1.2, 2.1)
  xTest <- c(-.4, 1.9, -0.8)
  Rans <- Radfun1$derivsRun(xTest)
  Cadfun1$derivsRun(xRec)
  Cans <- Cadfun1$derivsRun(xTest)
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)
})

test_that('AD works with calls to other methods and RCfunctions specified in a list', {
  adhelper <- nimbleFunction(
    run = function(x = double(1)) {
      ans <- 0
      for(i in 1:length(x))
        ans <- ans + exp(x[i])
      return(ans)
      returnType(double(0))
    },
    buildDerivs = list(run = list(ignore = 'i'))
  )
  temporarilyAssignInGlobalEnv(adhelper)
  adfun <- nimbleFunction(
    setup = TRUE,
    run = function(x = double(1)) {
      ans <- 2*mysqrt(adhelper(x))
      return(ans)
      returnType(double())
    },
    methods = list(
      mysqrt = function(x = double()) {
        return(sqrt(x))
        returnType(double())
      },
      derivsRun = function(x = double(1)) {
        ans <- nimDerivs(run(x), wrt = 1:length(x), order = 0:2)
        return(ans)
        returnType(ADNimbleList())
      }),
    buildDerivs = list(run = list(),
                       mysqrt = list())
  )
  temporarilyAssignInGlobalEnv(adfun)
  Radfun1 <- adfun()
  Cadfun1 <- compileNimble(Radfun1) # This fails if adhelper is included.
  xRec <- c(.3, -1.2, 2.1)
  xTest <- c(-.4, 1.9, -0.8)
  Rans <- Radfun1$derivsRun(xTest)
  Cadfun1$derivsRun(xRec)
  Cans <- Cadfun1$derivsRun(xTest)
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)
})

test_that('AD works with calls to other methods and RCfunctions specified in a list', {
  adhelper <- nimbleFunction(
    run = function(x = double(1)) {
      ans <- exp(x)
      return(ans)
      returnType(double(1))
    },
    buildDerivs = 'run'
  )
  temporarilyAssignInGlobalEnv(adhelper)
  adfun <- nimbleFunction(
    setup = TRUE,
    run = function(x = double(1)) {
      expx <- adhelper(x)
      ans <- 0
      for(i in 1:length(expx))
        ans <- ans + expx[i]
      ans <- 2*mysqrt(ans)
      return(ans)
      returnType(double())
    },
    methods = list(
      mysqrt = function(x = double()) {
        return(sqrt(x))
        returnType(double())
      },
      derivsRun = function(x = double(1)) {
        ans <- nimDerivs(run(x), wrt = 1:length(x), order = 0:2)
        return(ans)
        returnType(ADNimbleList())
      }),
    buildDerivs = list(run = list(ignore = 'i'),
                       mysqrt = list())
  )
  temporarilyAssignInGlobalEnv(adfun)
  Radfun1 <- adfun()
  Cadfun1 <- compileNimble(Radfun1) # This fails if adhelper is included.
  xRec <- c(.3, -1.2, 2.1)
  xTest <- c(-.4, 1.9, -0.8)
  Rans <- Radfun1$derivsRun(xTest)
  Cadfun1$derivsRun(xRec)
  Cans <- Cadfun1$derivsRun(xTest)
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)
})

test_that('AD works with a model', {
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:3)
        x[i] ~ dnorm(0, 1)
    }),
    buildDerivs = TRUE
  )
  temporarilyAssignInGlobalEnv(m)
  
  adfun <- nimbleFunction(
    setup = function(model, nodes) {
      infoNodes <- makeDerivsInfo(model, wrtNodes = nodes, calcNodes = nodes)
      updateNodes <- infoNodes$updateNodes
      # ignore constant nodes as both are empty anyway
    },
    run = function(x = double(1)) {
      values(model, nodes) <<- x
      ans <- model$calculate(nodes)
      return(ans)
      returnType(double())
    },
    methods = list(
      derivsRun = function(x = double(1)) {
        ans <- nimDerivs(run(x), wrt = 1:length(x), order = 0:2,
                         model = m, updateNodes = updateNodes)
        return(ans)
        returnType(ADNimbleList())
      }),
    buildDerivs = 'run'
  )
  temporarilyAssignInGlobalEnv(adfun)
  Radfun1 <- adfun(m, 'x[1:3]')
  Cm <- compileNimble(m)
  Cadfun1 <- compileNimble(Radfun1, project = m)
  xRec <- c(.3, -1.2, 2.1)
  xTest <- c(-.4, 1.9, -0.8)
  Rans <- Radfun1$derivsRun(xTest)
  Cadfun1$derivsRun(xRec)
  Cans <- Cadfun1$derivsRun(xTest)
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)
})

test_that('AD works with a model compiled in one call', {
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:3)
        x[i] ~ dnorm(0, 1)
    }),
    buildDerivs = TRUE
  )
  temporarilyAssignInGlobalEnv(m)
  
  adfun <- nimbleFunction(
    setup = function(model, nodes) {
      infoNodes <- makeDerivsInfo(model, wrtNodes = nodes, calcNodes = nodes)
      updateNodes <- infoNodes$updateNodes
      # ignore constant nodes as both are empty anyway
    },
    run = function(x = double(1)) {
      values(model, nodes) <<- x
      ans <- model$calculate(nodes)
      return(ans)
      returnType(double())
    },
    methods = list(
      derivsRun = function(x = double(1)) {
        ans <- nimDerivs(run(x), wrt = 1:length(x), order = 0:2,
                         model = m, updateNodes = updateNodes)
        return(ans)
        returnType(ADNimbleList())
      }),
    buildDerivs = 'run'
  )
  temporarilyAssignInGlobalEnv(adfun)
  Radfun1 <- adfun(m, 'x[1:3]')
  comp <- compileNimble(Radfun1, m)
  Cadfun1 <- comp[[1]]
  xRec <- c(.3, -1.2, 2.1)
  xTest <- c(-.4, 1.9, -0.8)
  Rans <- Radfun1$derivsRun(xTest)
  Cadfun1$derivsRun(xRec)
  Cans <- Cadfun1$derivsRun(xTest)
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)
})

test_that('AD works with a model calling a custom dist and custom fun', {
  foo <- nimbleFunction(
    run = function(x = double()) {
      return(exp(x))
      returnType(double())
    },
    buildDerivs = 'run')
  temporarilyAssignInGlobalEnv(foo)

  dmynorm <- nimbleFunction(
    run = function(x = double(), mu = double(), log = double(0, default = 0)) {
      return(dnorm(x, mu, 1, log = log))
      returnType(double())
    },
    buildDerivs = 'run')
  temporarilyAssignInGlobalEnv(dmynorm)

  m <- nimbleModel(
    nimbleCode({
      for(i in 1:3) {
        a[i] ~ dbeta(2, 2)
        x[i] ~ dmynorm(mu = foo(a[i]))
      }
    }),
    inits = list(a = c(.1, .2, .3), x = c(.2, .3, .4)),
    buildDerivs = TRUE
  )
  temporarilyAssignInGlobalEnv(m)
  
  adfun <- nimbleFunction(
    setup = function(model, wrt, calc) {
      infoNodes <- makeDerivsInfo(model, wrtNodes = wrt, calcNodes = calc)
      updateNodes <- infoNodes$updateNodes
      constantNodes <- infoNodes$constantNodes
    },
    run = function(x = double(1)) {
      values(model, wrt) <<- x
      ans <- model$calculate(calc)
      return(ans)
      returnType(double())
    },
    methods = list(
      derivsRun = function(x = double(1)) {
        ans <- nimDerivs(run(x), wrt = 1:length(x), order = 0:2,
                         model = m, updateNodes = updateNodes, constantNodes = constantNodes)
        return(ans)
        returnType(ADNimbleList())
      }),
    buildDerivs = 'run'
  )
  temporarilyAssignInGlobalEnv(adfun)

  wrt <- 'a[1:3]'
  Radfun1 <- adfun(m, wrt, m$getDependencies(wrt))
  Cm <- compileNimble(m)
  Cadfun1 <- compileNimble(Radfun1, project = m)
  xRec <- c(.3, .6, .9)
  xTest <- c(.8, .1, .4)
  Rans <- Radfun1$derivsRun(xTest)
  Cadfun1$derivsRun(xRec)
  Cans <- Cadfun1$derivsRun(xTest)
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)
})
