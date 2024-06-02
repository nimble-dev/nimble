# This file contains tests of various modes of creating and invoking AD.

source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)

RCrelTol0 <- 1e-12
RCrelTol1 <- 1e-6
RCrelTol2 <- 1e-4

test_that("AD error-trapping of unsupported dists and fxns in models works", {
  ## No easy way to check that there is output but that it doesn't include a given message.
  ## test_that disables `capture.output`. Can run these manually without use of test_that.
    
  ## origOption <- nimbleOptions('doADerrorTraps')

  ## nimbleOptions(buildModelDerivs = FALSE)
  ## expect_message(
  ##   m <- nimbleModel(
  ##   quote({
  ##     a ~ dcat(p[1:2])
  ##   })
  ##   ), "Defining modelBuilding model")

  ## expect_silent(
  ##   m <- nimbleModel(
  ##   quote({
  ##     a <- b %% c
  ##   })
  ## ))
  ## nimbleOptions(buildModelDerivs = TRUE)
  
  expect_message(
    m <- nimbleModel(
    quote({
      a <- b %% c
    }),
    buildDerivs = TRUE
  ), "is not supported for derivatives")

  expect_message(
    m <- nimbleModel(
    quote({
      a <- exp(b %% c)
    }),
    buildDerivs = TRUE
  ), "is not supported for derivatives")

  ## nimbleOptions(doADerrorTraps = FALSE)
  
  ## expect_silent(
  ##   m <- nimbleModel(
  ##   quote({
  ##     a ~ dcat(p[1:2])
  ##   }),
  ##   buildDerivs = TRUE
  ## ))

  ## expect_silent(
  ##   m <- nimbleModel(
  ##   quote({
  ##     a <- b %% c
  ##   }),
  ##   buildDerivs = TRUE
  ## ))

  ## nimbleOptions(doADerrorTraps = origOption)

})

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

test_that('AD works with ADbreak with setup = TRUE', {
  adfun <- nimbleFunction(
    setup = TRUE,
    run = function(x = double(1)) {
      ans <- 0
      for(i in 1:length(x))
        ans <- ans + exp(x[i])
      ans1 <- ADbreak(ans) # ans1 is CppAD::AD but assigned only value
      ans2 <- ADbreak(ans) # ans2 is double
      myint <- 0L
      myint <- ADbreak(ans) # assignment to integer also requires ADbreak
      ans <- ans2          # can assign from double to CppAD::AD
      return(ans)
      returnType(double())
    },
    methods = list(
      derivsRun = function(x = double(1)) {
        ans <- nimDerivs(run(x), wrt = 1:length(x), order = 0:2)
        return(ans)
        returnType(ADNimbleList())
      }),
    buildDerivs = list(run = list(ignore = c('i', 'ans2')))
  )
  temporarilyAssignInGlobalEnv(adfun)
  Radfun1 <- adfun()
  Cadfun1 <- compileNimble(Radfun1)
  xRec <- c(.3, -1.2, 2.1)
  xTest <- c(-.4, 1.9, -0.8)
  Rans <- Radfun1$derivsRun(xTest) # does not respect ADbreak
  CansRec <- Cadfun1$derivsRun(xRec) # xRec gets baked in due to ADbreak
  Cans <- Cadfun1$derivsRun(xTest)
  nim_expect_equal(CansRec$value, Cans$value, RCrelTol0)
  nim_expect_equal(Cans$jacobian, matrix(rep(0, 3), ncol = 3), RCrelTol1)
  nim_expect_equal(Cans$hessian, array(rep(0, 9), dim=c(3,3,1)), RCrelTol2)
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

test_that('AD works with model$calulate in a method', {
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
      infoNodes <- makeModelDerivsInfo(model, wrtNodes = nodes, calcNodes = nodes)
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

test_that('AD works with nimDerivs(model$calulate(...), ...)', {
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:3)
        x[i] ~ dnorm(0, 1)
    }),
    buildDerivs = TRUE
  )
  temporarilyAssignInGlobalEnv(m)

  adfun <- nimbleFunction(
    setup = function(model, wrt, calcNodes) {},
    methods = list(
      derivsRun = function() {
        ans <- nimDerivs(model$calculate(calcNodes), wrt = wrt, order = 0:2)
        return(ans)
        returnType(ADNimbleList())
      })
  )
  temporarilyAssignInGlobalEnv(adfun)

  Radfun1 <- adfun(m,
                   wrt = 'x[1:3]',
                   calcNodes = 'x[1:2]') # Deliberately using only x[1:2] as a tweak
  Cm <- compileNimble(m)
  Cadfun1 <- compileNimble(Radfun1, project = m)
  xRec <- c(.3, -1.2, 2.1)
  xTest <- c(-.4, 1.9, -0.8)
  m$x <- xTest
  Rans <- Radfun1$derivsRun()
  Cm$x <- xRec
  Cadfun1$derivsRun()
  Cm$x <- xTest
  Cans <- Cadfun1$derivsRun()
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
      infoNodes <- makeModelDerivsInfo(model, wrtNodes = nodes, calcNodes = nodes)
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

  ## Something goes wrong finding the automatic version of this in test_that evaluation,
  ## so we provide it explicitly
  rmynorm <- nimbleFunction(
    run = function(n = integer(), mu = double()) {
      return(rnorm(1, mu))
      returnType(double())
    })
  temporarilyAssignInGlobalEnv(rmynorm)

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
      infoNodes <- makeModelDerivsInfo(model, wrtNodes = wrt, calcNodes = calc)
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
  deregisterDistributions("dmynorm")
})

test_that('AD works with reset', {
  adfun <- nimbleFunction(
    setup = TRUE,
    run = function(x = double(1)) {
      return(sum(exp(x)))
      returnType(double())
    },
    methods = list(
      derivsRun = function(x = double(1),
                           reset = logical(0, default = FALSE)) {
        ans <- nimDerivs(run(x), wrt = 1:length(x), order = 0:2, reset = reset)
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

  xRec2 <- c(.1, .2, .3, .4)
  xTest2 <- xRec2 + .5
  Rans <- Radfun1$derivsRun(xTest2)
  Cadfun1$derivsRun(xRec2, reset = TRUE)
  Cans <- Cadfun1$derivsRun(xTest2)
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)
})

test_that('AD works with reset with empty model args', {
  adfun <- nimbleFunction(
    setup = TRUE,
    run = function(x = double(1)) {
      return(sum(exp(x)))
      returnType(double())
    },
    methods = list(
      derivsRun = function(x = double(1),
                           reset = logical(0, default = FALSE)) {
        ans <- nimDerivs(run(x), wrt = 1:length(x), order = 0:2, reset = reset,
                         model = NA, updateNodes = NA, constantNodes = NA)
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

  xRec2 <- c(.1, .2, .3, .4)
  xTest2 <- xRec2 + .5
  Rans <- Radfun1$derivsRun(xTest2)
  Cadfun1$derivsRun(xRec2, reset = TRUE)
  Cans <- Cadfun1$derivsRun(xTest2)
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)
})

test_that('AD works with a model with do_update', {
  m <- nimbleModel(
    nimbleCode({
      mu ~ dnorm(0, 1)
      for(i in 1:3)
        x[i] ~ dnorm(mu, 1)
    }),
    buildDerivs = TRUE
  )
  temporarilyAssignInGlobalEnv(m)

  adfun <- nimbleFunction(
    setup = function(model, nodes) {
      infoNodes <- makeModelDerivsInfo(model, wrtNodes = nodes, calcNodes = nodes)
      updateNodes <- infoNodes$updateNodes
      # constantNodes should be empty
    },
    run = function(x = double(1)) {
      values(model, nodes) <<- x
      ans <- model$calculate(nodes)
      return(ans)
      returnType(double())
    },
    methods = list(
      derivsRun = function(x = double(1),
                           do_update = logical(0, default = TRUE),
                           reset = logical(0, default = FALSE)) {
        ans <- nimDerivs(run(x), wrt = 1:length(x), order = 0:2,
                         do_update = do_update, reset = reset,
                         model = m, updateNodes = updateNodes)
        return(ans)
        returnType(ADNimbleList())
      }),
    buildDerivs = 'run'
  )

  temporarilyAssignInGlobalEnv(adfun)
  m$mu <- 0.85
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

  expect_equal(Radfun1$updateNodes, "mu")
  Cm$mu <- 1.15
  Cans <- Cadfun1$derivsRun(xTest, do_update = FALSE) # should match above with mu = 0.85
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)

  m$mu <- 1.15
  Rans <- Radfun1$derivsRun(xTest)
  Cans <- Cadfun1$derivsRun(xTest, do_update = TRUE) # should now use mu = 1.15
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)
})

test_that('AD works with a model with do_update', {
  m <- nimbleModel(
    nimbleCode({
      mu ~ dnorm(0, 1)
      for(i in 1:3)
        x[i] ~ dnorm(mu, 1)
    }),
    buildDerivs = TRUE
  )
  temporarilyAssignInGlobalEnv(m)

  adfun <- nimbleFunction(
    setup = function(model, nodes) {
      infoNodes <- makeModelDerivsInfo(model, wrtNodes = nodes, calcNodes = nodes)
      updateNodes <- infoNodes$updateNodes
      # constantNodes should be empty
    },
    run = function(x = double(1)) {
      values(model, nodes) <<- x
      ans <- model$calculate(nodes)
      return(ans)
      returnType(double())
    },
    methods = list(
      derivsRun = function(x = double(1),
                           do_update = logical(0, default = TRUE),
                           reset = logical(0, default = FALSE)) {
        ans <- nimDerivs(run(x), wrt = 1:length(x), order = 0:2,
                         do_update = do_update, reset = reset,
                         model = m, updateNodes = updateNodes)
        return(ans)
        returnType(ADNimbleList())
      }),
    buildDerivs = 'run'
  )

  temporarilyAssignInGlobalEnv(adfun)
  m$mu <- 0.85
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

  ## Changing mu in the model but then not updating should not change derivs.
  expect_equal(Radfun1$updateNodes, "mu")
  Cm$mu <- 1.15
  Cans <- Cadfun1$derivsRun(xTest, do_update = FALSE) # should match above with mu = 0.85
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)

  ## Changing my in the model and then updating should change derivs.
  m$mu <- 1.15
  Rans <- Radfun1$derivsRun(xTest)
  Cans <- Cadfun1$derivsRun(xTest, do_update = TRUE) # should now use mu = 1.15
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)
})

test_that('AD works with a model with do_update in double-taping', {
  m <- nimbleModel(
    nimbleCode({
      mu ~ dnorm(0, 1)
      for(i in 1:3)
        x[i] ~ dnorm(mu, 1)
    }),
    buildDerivs = TRUE
  )
  temporarilyAssignInGlobalEnv(m)

  adfun <- nimbleFunction(
    setup = function(model, nodes) {
      infoNodes <- makeModelDerivsInfo(model, wrtNodes = nodes, calcNodes = nodes)
      updateNodes <- infoNodes$updateNodes
      # constantNodes should be empty
    },
    run = function(x = double(1)) {
      values(model, nodes) <<- x
      ans <- model$calculate(nodes)
      return(ans)
      returnType(double())
    },
    methods = list(
      valueRun = function(x = double(1),
                           do_update = logical(0, default = TRUE),
                           reset = logical(0, default = FALSE)) {
        ans <- nimDerivs(run(x), wrt = 1:length(x), order = 0,
                         do_update = do_update, reset = reset,
                         model = m, updateNodes = updateNodes)
        return(ans$value[1])
        returnType(double())
      },
      derivsValueRun = function(x = double(1),
                                do_update = logical(0, default = TRUE),
                                inner_do_update = logical(0, default = TRUE),
                                reset = logical(0, default = FALSE),
                                inner_reset = logical(0, default = FALSE)) {
        ans <- nimDerivs(valueRun(x, inner_do_update, inner_reset),
                         wrt = 1:length(x),
                         do_update = do_update, reset = reset,
                         model = m, updateNodes = updateNodes)
        return(ans)
        returnType(ADNimbleList())
      }
      ),
    buildDerivs = c('run', 'valueRun')
  )
  temporarilyAssignInGlobalEnv(adfun)

  m$mu <- 0.85
  Radfun1 <- adfun(m, 'x[1:3]')
  Cm <- compileNimble(m)
  Cadfun1 <- compileNimble(Radfun1, project = m)
  xRec <- c(.3, -1.2, 2.1)
  xTest <- c(-.4, 1.9, -0.8)

  # Basic check
  Rans <- Radfun1$derivsValueRun(xTest)
  Cadfun1$derivsValueRun(xRec)
  Cans <- Cadfun1$derivsValueRun(xTest)
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)

  # Changing mu in the model and turning off inner update should still do outer update and change derivs.
  expect_equal(Radfun1$updateNodes, "mu")
  Cm$mu <- 1.15
  Cans <- Cadfun1$derivsValueRun(xTest, inner_do_update = FALSE, do_update = TRUE)
  m$mu <- 1.15
  Rans <- Radfun1$derivsValueRun(xTest)
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)

  # Changing mu in the model and turning off outer update should not change derivs.
  Cm$mu <- 1.05
  Cans <- Cadfun1$derivsValueRun(xTest, inner_do_update = FALSE, do_update = FALSE)
  m$mu <- 1.15
  Rans <- Radfun1$derivsValueRun(xTest)
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)

  # Changing mu in the model and turning off outer update should not change derivs even if inner update is on.
  Cm$mu <- 1.05
  Cans <- Cadfun1$derivsValueRun(xTest, inner_do_update = TRUE, do_update = FALSE)
  m$mu <- 1.15
  Rans <- Radfun1$derivsValueRun(xTest)
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)

  # Changing mu and turning on outer update with inner update on should also change derivs
  Cm$mu <- 1.05
  Cans <- Cadfun1$derivsValueRun(xTest, inner_do_update = TRUE, do_update = TRUE)
  m$mu <- 1.05
  Rans <- Radfun1$derivsValueRun(xTest)
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)

  # Resetting inner tape even without do_update should get new update values
  Cm$mu <- 0.95
  Cans <- Cadfun1$derivsValueRun(xTest, inner_do_update = FALSE, inner_reset = TRUE, do_update = TRUE)
  m$mu <- 0.95
  Rans <- Radfun1$derivsValueRun(xTest)
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)

  # Resetting inner tape even without do_update should get previous update values from outer tape if outer do_update if FALSE
  Cm$mu <- 0.65
  Cans <- Cadfun1$derivsValueRun(xTest, inner_do_update = FALSE, inner_reset = TRUE, do_update = FALSE)
  m$mu <- 0.95
  Rans <- Radfun1$derivsValueRun(xTest)
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)

  ## Resetting outer tape but not inner tape is really weird.
  ## The outer tape is cleared, and the inner tape used to re-record with the values it
  ## had last time it was recorded, which was 1.15 above.
  ## This is hard to follow but could be almost undefined behavior.
  ## This test is commented out because it is basically weird and unsupported
  ##   on the uncompiled side.
  ## Cm$mu <- 0.8
  ## Cans <- Cadfun1$derivsValueRun(xTest, inner_do_update = FALSE, inner_reset = FALSE,
  ##                                do_update = TRUE, reset = TRUE)
  ## m$mu <- 1.15
  ## Rans <- Radfun1$derivsValueRun(xTest)
  ## nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  ## nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  ## nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)

  ## This case is also weird:
  ## Cm$mu <- 0.9
  ## Cans <- Cadfun1$derivsValueRun(xTest, inner_do_update = FALSE, inner_reset = TRUE,
  ##                               do_update = FALSE, reset = TRUE)

  ## Resetting both tapes resets all updates
  Cm$mu <- 0.9
  Cans <- Cadfun1$derivsValueRun(xTest, inner_do_update = TRUE, inner_reset = TRUE,
                                 do_update = TRUE, reset = TRUE)
  m$mu <- 0.9
  Rans <- Radfun1$derivsValueRun(xTest)
  nim_expect_equal(Rans$value, Cans$value, RCrelTol0)
  nim_expect_equal(Rans$jacobian, Cans$jacobian, RCrelTol1)
  nim_expect_equal(Rans$hessian, Cans$hessian, RCrelTol2)
})


## Following are some draft tests relating to how we handle
## non-AD functions within an AD-enabled model (and functions)

test_that('non-AD function in AD-enabled model', {
    dfoo <- nimbleFunction(
        run = function(x = double(), k = double(), log = integer(0, default = 0)) {
            return(-k * x)
            returnType(double())
        },
        buildDerivs= FALSE
    )
    temporarilyAssignInGlobalEnv(dfoo)
    mc <- nimbleCode({
        a ~ dnorm(0, 1)
        z ~ dfoo(a)
    })

    expect_message(m <- nimbleModel(mc, data = list(z = 2.3), buildDerivs = TRUE),
                   "Distribution dfoo does not appear to support derivatives")
})


test_that("derivatives work if a non-supported user-defined call is not included", {
  dfoo2 <- nimbleFunction(
    run = function(x = double(), k = double(), log = integer(0, default = 0)) {
      return(-k * x)
      returnType(double())
    },
    buildDerivs= FALSE
  )
  temporarilyAssignInGlobalEnv(dfoo2)
  mc <- nimbleCode({
    a ~ dnorm(0, 1)
    z ~ dfoo2(a)
  })

  m <- nimbleModel(mc, data = list(z = 2.3), inits = list(a = 1.5),
                   buildDerivs = TRUE)
  temporarilyAssignInGlobalEnv(rfoo2)
  cm <- compileNimble(m)

  test <- nimbleFunction(
    setup = function(m, nodes) {},
    run = function() {
      ans <- nimDerivs(m$calculate(nodes), wrt = "a", order = 0:1 )
      return(ans)
      returnType(ADNimbleList())
    }
  )

  test1 <- test(m, "a")
  ctest1 <- compileNimble(test1, project = m)
  ans <- ctest1$run()
  expect_equal( round(ans$value, digits = 2), -2.04 )
  expect_equal( round(ans$jacobian, digits = 2), matrix(-1.5) )
})

test_that("C++ run-time error if non-supported user-defined call is included in model$calculate", {
  dfoo3 <- nimbleFunction(
    run = function(x = double(), k = double(), log = integer(0, default = 0)) {
      return(-k * x)
      returnType(double())
    },
    buildDerivs= FALSE
  )
  temporarilyAssignInGlobalEnv(dfoo3)

  mc <- nimbleCode({
    a ~ dnorm(0, 1)
    z ~ dfoo3(a)
  })

  m <- nimbleModel(mc, data = list(z = 2.3), inits = list(a = 1.5),
                   buildDerivs = TRUE)
  temporarilyAssignInGlobalEnv(rfoo3)
  cm <- compileNimble(m)

  test <- nimbleFunction(
    setup = function(m, nodes) {},
    run = function() {
      ans <- nimDerivs(m$calculate(nodes), wrt = "a", order = 0:1 )
      return(ans)
      returnType(ADNimbleList())
    }
  )

  test1 <- test(m, c("a", "z"))
  ctest1 <- compileNimble(test1, project = m)
  expect_error(ans <- ctest1$run(), "a model node is being included in derivatives that does not support them")
  deregisterDistributions("dfoo3")
})

test_that("C++ run-time error if non-supported user-defined call is included in model$calculate", {
  mc <- nimbleCode({
    a ~ dnorm(0, 1)
    z ~ dconstraint(a > 1)
  })

  m <- nimbleModel(mc, data = list(z = 2.3), inits = list(a = 1.5),
                   buildDerivs = TRUE)
  cm <- compileNimble(m)

  test <- nimbleFunction(
    setup = function(m, nodes) {},
    run = function() {
      ans <- nimDerivs(m$calculate(nodes), wrt = "a", order = 0:1 )
      return(ans)
      returnType(ADNimbleList())
    }
  )

  test1 <- test(m, c("a", "z"))
  ctest1 <- compileNimble(test1, project = m)
  expect_error(ctest1$run(), "a model node is being included in derivatives that does not support them")
})


test_that("warn if non-supported user-defined call is included in model$calculate", {
    nimbleOptions(unsupportedDerivativeHandling='warn')
    mc <- nimbleCode({
        a ~ dnorm(0, 1)
        z ~ dconstraint(a > 1)
    })

    m <- nimbleModel(mc, data = list(z = 2.3), inits = list(a = 1.5),
                     buildDerivs = TRUE)
    cm <- compileNimble(m)
    nimbleOptions(unsupportedDerivativeHandling='error')

    test <- nimbleFunction(
        setup = function(m, nodes) {},
        run = function() {
            ans <- nimDerivs(m$calculate(nodes), wrt = "a", order = 0:1 )
            return(ans)
            returnType(ADNimbleList())
        }
    )

    test1 <- test(m, c("a", "z"))
    ctest1 <- compileNimble(test1, project = m)
    expect_output(ctest1$run(),
                   "a model node is being included in derivatives that does not support them")
    nimbleOptions(unsupportedDerivativeHandling='error')
})
