source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

context("Testing of expandNodeNames")

test_that("expandNodeNames works for various cases, including going beyond extent of variable", {
    
   code <- nimbleCode({
    for(i in 1:4)
        mu[i] ~ dnorm(0,1)
    for(i in 1:3)
        for(j in 1:3)
            theta[i,j] ~ dnorm(0,1)
    p[1:4] ~ ddirch(alpha[1:4])
   })
   m <- nimbleModel(code, inits = list(alpha = rep(1, 4)))

   ## vector variable
   expect_equal(m$expandNodeNames("mu"), c("mu[1]","mu[2]","mu[3]","mu[4]"))
   expect_equal(m$expandNodeNames("mu[3:5]"), c("mu[3]","mu[4]"))
   expect_equal(m$expandNodeNames("mu[5:7]"), character(0))

   ## matrix variable
   expect_equal(m$expandNodeNames("theta"), c("theta[1, 1]","theta[2, 1]","theta[3, 1]","theta[1, 2]","theta[2, 2]","theta[3, 2]","theta[1, 3]","theta[2, 3]","theta[3, 3]"))
   expect_equal(m$expandNodeNames("theta[3:5,1:2]"), c("theta[3, 1]","theta[3, 2]"))
   expect_equal(m$expandNodeNames("theta[1:2,3:5]"), c("theta[1, 3]","theta[2, 3]"))
   expect_equal(m$expandNodeNames("theta[4:6,5]"), character(0))

   ## multiple inputs, mixed
   expect_equal(m$expandNodeNames(c("theta[1, 7]", "mu")), c("mu[1]","mu[2]","mu[3]","mu[4]"))
   expect_equal(m$expandNodeNames(c("theta[1, 3:5]", "mu[3:5]")), c("theta[1, 3]", "mu[3]", "mu[4]"))
   expect_equal(m$expandNodeNames(c("theta[1, 7]", "mu[5]")), character(0))
   expect_equal(m$expandNodeNames(c("mu[1, 1]", "theta[1, 1]")), "theta[1, 1]")
                
   ## multiple inputs, mixed, not unique
   expect_equal(m$expandNodeNames(c("mu[3:5]", "mu[3:9]"), unique = FALSE),
                c("mu[3]","mu[4]","mu[3]","mu[4]"))
   expect_equal(m$expandNodeNames(c("theta[1:3, 3:5]", "theta[3:5, 1:3]"), unique = FALSE),
                c("theta[1, 3]","theta[2, 3]","theta[3, 3]","theta[3, 1]","theta[3, 2]","theta[3, 3]"))
                                  
   ## indexing with a variable
   expect_error(m$expandNodeNames("theta[[a]]"), "variable was found in the indexing")

   ## multivariate node
   expect_equal(m$expandNodeNames("p[1:4]"), "p[1:4]")
   expect_equal(m$expandNodeNames("p[1:2]"), "p[1:4]")
   expect_equal(m$expandNodeNames("p[1:5]"), "p[1:4]")
   expect_equal(m$expandNodeNames("p[5:7]"), character(0))
   expect_equal(m$expandNodeNames(c("p[1:5]", "mu[3:5]")), c("p[1:4]", "mu[3]", "mu[4]"))
   expect_equal(m$expandNodeNames(c("p[1:5]", "p"), unique = FALSE), c("p[1:4]", "p[1:4]"))              
})

test_that("Node index 100000 is not handled as 'x[1e+5]'", {
  set.seed(1)
  mc <- nimbleCode({
    for(i in 99998:100002) {
      x[i] ~ dnorm(0,1)
    }
  })
  m <- nimbleModel(mc, inits = list(x = 1:100002), calculate = FALSE)
  expect_identical(m$getNodeNames()[3], "x[100000]")
  expect_identical(m$expandNodeNames("x[1e+5]"), "x[100000]")
  Cm <- compileNimble(m)

  nf <- nimbleFunction(
    setup = function(model) {
      node <- 'x[100000]'
      nodeSci <- 'x[1e+5]'
      mv <- modelValues(model)
    },
    run = function() {},
    methods = list(
      getVal1 = function() {return(model[[node]]); returnType(double())},
      getVal2 = function() {return(model[[nodeSci]]); returnType(double())},
      getVal3 = function() {return(model[['x[1e+5]']]); returnType(double())},
      getVal4 = function() {v <- values(model, node); return(v[1]); returnType(double())},
      #    getVal5 = function() {v <- values(model, nodeSci); return(v[1]); returnType(double())}, # not supported in compilation
      #    getVal6 = function() {v <- values(model, 'x[1e+5]'); return(v[1]); returnType(double())}, # not supported in compilation
      calc1 = function() {return(model$calculate(node)); returnType(double())},
      calc2 = function() {return(model$calculate(nodeSci)); returnType(double())},
      calc3 = function() {return(model$calculate('x[1e+5]')); returnType(double())},
      calc4 = function() {return(model$calculate(node)); returnType(double())},
      copy1 = function() {nimCopy(from = model, to = mv, rowTo = 1, nodes = node, logProb = TRUE)}
      #    copy2 = function() {nimCopy(from = model, to = mv, rowTo = 1, nodes = nodeSci, logProb = TRUE)}, # not supported in compilation
      #    copy3 = function() {nimCopy(from = model, to = mv, rowTo = 1, nodes = 'x[1e+5]', logProb = TRUE)} # not supported in compilation
    )
  )

  nf1 <- nf(m)
  Cnf1 <- compileNimble(nf1, project = m)

  expect_identical(m$getDependencies("x[1e+5]"), 'x[100000]')
  expect_identical(m$expandNodeNames("x[1e+5]"), 'x[100000]')

  m$x[100000] <- 1
  Cm$x[100000] <- 2
  for(ver in 1:2) {
    if(ver == 1) {
      mod <- m
      fun <- nf1
    } else {
      mod <- Cm
      fun = Cnf1
    }
    i <- 0
    i <- i+1; mod$x[100000] <- i; expect_identical(fun$getVal1(), i)
    i <- i+1; mod$x[100000] <- i; expect_identical(fun$getVal2(), i)
    i <- i+1; mod$x[100000] <- i; expect_identical(fun$getVal3(), i)
    i <- i+1; mod$x[100000] <- i; expect_identical(fun$getVal4(), i)
    i <- 0
    i <- i+1; mod$x[100000] <- i; expect_identical(fun$calc1(), dnorm(i, 0, 1, log = TRUE))
    i <- i+1; mod$x[100000] <- i; expect_identical(fun$calc2(), dnorm(i, 0, 1, log = TRUE))
    i <- i+1; mod$x[100000] <- i; expect_identical(fun$calc3(), dnorm(i, 0, 1, log = TRUE))
    i <- i+1; mod$x[100000] <- i; expect_identical(fun$calc4(), dnorm(i, 0, 1, log = TRUE))
    i <- 0
    i <- i+1; mod$x[100000] <- i; mod$calculate(); fun$copy1(); expect_identical(i, fun$mv['x', 1][100000])
  }
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)

   

