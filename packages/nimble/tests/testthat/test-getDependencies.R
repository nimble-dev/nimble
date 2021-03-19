source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

context("Testing of getDependencies")

## note this testing is not intended to blur into general model processing testing.
## It assumes the model processing is ok and really tests traversal of the graph to get dependencies

## also self=FALSE is not a deep testing need because this is processed in R after the deeper processing (graph traversal) in C++

## first model, with no criss-crossing dependencies.
test_that("getDependencies in first model with no criss-crossing dependencies", {
    c1 <- nimbleCode({
        a1 ~ dnorm(0,1)
        a2 ~ dnorm(0,1)
        b1 ~ dnorm(a1, 1) ## only a1 
        b2 ~ dnorm(a1, sd = a2) ## a1 and a2
        c1 <- a1 + 1 
        d1 ~ dnorm(c1, 1) ## only a1 via c1
        c2 <- a1 + 2
        c3 <- c2^2
        e1 ~ dnorm(c3, 1) ## only a1 via c2 and d1
        c4 <- a2 + 3
        f1 ~ dnorm(c3, sd = c4) ## a1 via c1 and d1; a2 via c3
        g1 ~ dnorm(f1, 1)
    })
    
    m1 <- nimbleModel(c1)
    ans1 <- m1$topologicallySortNodes(c('a1','b1','b2','c1','d1','c2','c3','e1','f1'))
    
    ## basic cases from a1
    expect_identical(m1$getDependencies('a1'), ans1)
    expect_identical(m1$getDependencies(c('a1','a1')), ans1)
    expect_identical(m1$getDependencies(c('a1','c1')), ans1)
    expect_identical(m1$getDependencies(c('c1','a1')), ans1)
    expect_identical(m1$getDependencies(c('a1', 'f1')), c(ans1, 'g1'))
    expect_identical(m1$getDependencies(c('f1', 'a1')), c(ans1, 'g1'))
    expect_identical(m1$getDependencies(c('a1'), downstream=TRUE), c(ans1, 'g1'))
    
    ## basic cases from a1 and a2
    ans2 <- m1$topologicallySortNodes(c('a1','a2','b1','b2','c4','c1','d1','c2','c3','e1','f1'))
    expect_identical(m1$getDependencies(c('a1','a2')), ans2)
    expect_identical(m1$getDependencies(c('a2','a1')), ans2)
    
    ## some omit cases
    ans3 <- m1$topologicallySortNodes(c('a1','b1','b2','c1','d1','c2'))
    expect_identical(m1$getDependencies('a1', omit = 'c3'), ans3)
    
    ans4 <- m1$topologicallySortNodes(c('a1','a2','b1','b2','c1','d1','c2','c3','e1','f1'))
    expect_identical(m1$getDependencies(c('a2','a1'), omit = 'c4', downstream = TRUE), c(ans4, 'g1'))
    
    ## some cases starting from deterministic
    expect_identical(m1$getDependencies('c1'), c('c1','d1'))
    expect_identical(m1$getDependencies('c2'), c('c2','c3','e1','f1'))
    expect_identical(m1$getDependencies(c('c1','c2')), c('c1','c2','d1','c3','e1','f1'))
})
    
################
    ## second model, with some criss-crossing dependencies
test_that("getDependencies in second model with some criss-crossing dependencies", {
    c2 <- nimbleCode({
        a1 <- 1
        b1 ~ dnorm(a1, 1)
        c1 <- b1 + 1
        c2 <- b1 + 2
        d1 ~ dnorm(c1, sd = c2)
        e1 ~ dnorm(d1, 1)
        f1 ~ dnorm(d1, sd = c1)
        g1 ~ dnorm(e1, sd = c1)
        c3 <- c1 + c2
        h1 ~ dnorm(c3, 1)
        h2 ~ dnorm(c3, sd = c1)
    })

    m2 <- nimbleModel(c2)

    ans1 <- m2$topologicallySortNodes(c('a1','b1'))
    expect_identical(m2$getDependencies('a1'), ans1)

    ans2 <- m2$topologicallySortNodes(c('b1','c1','c2','d1','f1','g1','c3','h1','h2'))
    expect_identical(m2$getDependencies('b1'), ans2)
    expect_identical(m2$getDependencies(c('b1','c1')), ans2)
    expect_identical(m2$getDependencies(c('c1','b1')), ans2)
    expect_identical(m2$getDependencies(c('b1','c1','f1')), ans2)
    expect_identical(m2$getDependencies(c('f1','c1','b1')), ans2)
    expect_identical(m2$getDependencies(c('f1','b1','c1')), ans2)

    ans3 <- m2$topologicallySortNodes(c('b1','c1','c2','d1','e1','f1','g1','c3','h1','h2'))
    expect_identical(m2$getDependencies(c('b1','c1','d1')), ans3)
    expect_identical(m2$getDependencies(rev(ans2)), ans3)
    expect_identical(m2$getDependencies(c('b1','c1','e1')), ans3)
    expect_identical(m2$getDependencies(c('e1','c1','b1')), ans3)

    ans4 <- m2$topologicallySortNodes(c('c2','d1','c3','h1','h2'))
    expect_identical(m2$getDependencies(c('c2')), ans4)
    expect_identical(m2$getDependencies(c('c3','c2')), ans4)
    expect_identical(m2$getDependencies(c('c2','h1','c3')), ans4)
})

test_that("old getParentNodes works", {
    code=nimbleCode({
        for(j in 1:J) {
            for(i in 1:I)
                y[j,i] ~ dnorm(theta[j], sigma)
            theta[j] ~ dnorm(mu, sd = tau)
        }
        mu ~ dnorm(mu0,1)
        mu0 <- mu00
        sigma ~ dunif(sigma0,1)
        sigma1 ~ dunif(0,1)
        tau ~ dunif(0,1)
    })
    I <- 4
    J <- 3
    constants <- list(I = I, J = J)
    m <- nimbleModel(code,
                     data = list(y = matrix(rnorm(I*J), J, I)),
                     constants = constants)
    expect_identical(nimble:::getParentNodes('mu', m, stochOnly =  TRUE), character(0))
    expect_identical(nimble:::getParentNodes('mu', m), c('mu0', 'mu00'))
    expect_identical(nimble:::getParentNodes('y', m, stochOnly = TRUE),
                     c('theta[1]','theta[2]','theta[3]', 'sigma'))
    expect_identical(nimble:::getParentNodes('y', m),
                     c('lifted_d1_over_sqrt_oPsigma_cP', 'theta[1]','theta[2]','theta[3]', 'sigma'))
    expect_identical(nimble:::getParentNodes('y[2, 1:3]', m, stochOnly = TRUE),
                     c('theta[2]', 'sigma'))
    expect_identical(nimble:::getParentNodes('theta', m),
                     c('tau', 'mu'))
    
})


## First model from test-getDependencies
test_that("getParents works in model with no criss-crossing dependencies", {
  c1 <- nimbleCode({
    a1 ~ dnorm(0,1)
    a2 ~ dnorm(0,1)
    b1 ~ dnorm(a1, 1) ## only a1 
    b2 ~ dnorm(a1, sd = a2) ## a1 and a2
    c1 <- a1 + 1 
    d1 ~ dnorm(c1, 1) ## only a1 via c1
    c2 <- a1 + 2
    c3 <- c2^2
    e1 ~ dnorm(c3, 1) ## only a1 via c2 and d1
    c4 <- a2 + 3
    f1 ~ dnorm(c3, sd = c4) ## a1 via c1 and d1; a2 via c3
    g1 ~ dnorm(f1, 1)
    h1 ~ dnorm(lho, 1) # left-hand-side only
  })
  m1 <- nimbleModel(c1)

  expect_identical(m1$getParents("f1"), c("a1", "a2", "c2", "c4", "c3"))
  expect_identical(m1$getParents("f1", immediateOnly = TRUE), c("c4", "c3"))
  expect_identical(m1$getParents(c("f1", "g1"), stochOnly = TRUE), c("a1", "a2", "f1"))
  expect_identical(m1$getParents(c("f1", "g1"), immediateOnly = TRUE), c("c4", "c3", "f1"))
  expect_identical(m1$getParents(c("g1", "f1"), stochOnly = TRUE), c("a1", "a2", "f1"))
  expect_identical(m1$getParents(c("f1", "g1"), stochOnly = TRUE, self = TRUE), c("a1", "a2", "f1", "g1"))
  expect_identical(m1$getParents(c("c3", "c2", "e1"), stochOnly = FALSE), c("a1", "c2", "c3"))
  expect_identical(m1$getParents(c("c3", "c2", "e1"), immediateOnly = TRUE, stochOnly = FALSE), c("a1", "c2", "c3"))
  expect_identical(m1$getParents("h1", includeRHSonly = TRUE, stochOnly = FALSE), c("lho"))
})

test_that("getParents works in model with some criss-crossing dependencies", {
  ## Second model from test-getDependencies
  c2 <- nimbleCode({
    a1 <- 1
    b1 ~ dnorm(a1, 1)
    c1 <- b1 + 1
    c2 <- b1 + 2
    d1 ~ dnorm(c1, sd = c2)
    e1 ~ dnorm(d1, 1)
    f1 ~ dnorm(d1, sd = c1)
    g1 ~ dnorm(e1, sd = c1)
    c3 <- c1 + c2
    h1 ~ dnorm(c3, 1)
    h2 ~ dnorm(c3, sd = c1)
  })
  m2 <- nimbleModel(c2)

  expect_identical(m2$getParents("h2", stochOnly = TRUE), "b1")
  expect_identical(m2$getParents("h2", stochOnly = FALSE), c("b1", "c1", "c2", "c3"))
  expect_identical(m2$getParents("h2", stochOnly = FALSE, determOnly = TRUE), c("c1", "c2", "c3"))
  expect_identical(m2$getParents("h2", stochOnly = TRUE, self = TRUE), c("b1", "h2"))
  expect_identical(m2$getParents("g1", stochOnly = TRUE), c("b1", "e1"))
  expect_identical(m2$getParents("g1", stochOnly = FALSE), c("b1", "c1", "e1"))
  expect_identical(m2$getParents("b1", stochOnly = TRUE), character())
  expect_identical(m2$getParents("b1", stochOnly = FALSE), c("a1"))
  expect_identical(m2$getParents("f1", stochOnly = TRUE), c("b1", "d1"))
  expect_identical(m2$getParents("f1", stochOnly = FALSE), c("b1", "c1", "d1"))
  expect_identical(m2$getParents(c("e1", "h1"), stochOnly = TRUE), c("b1", "d1"))
  expect_identical(m2$getParents("d1", stochOnly = TRUE), c("b1"))
  expect_identical(m2$getParents("d1", stochOnly = TRUE, upstream = TRUE ), c("b1"))
  expect_identical(m2$getParents("d1", upstream = TRUE, stochOnly = FALSE ),
                   c("a1", "b1","c1","c2"))
  expect_identical(m2$getParents("e1", upstream = TRUE, stochOnly = FALSE ),
                   c("a1", "b1", "c1", "c2", "d1"))
  expect_identical(m2$getParents("e1", omit = "d1", upstream = TRUE, stochOnly = FALSE ),
                   character()) 
  expect_identical(m2$getParents("e1", omit = "d1", stochOnly = FALSE ),
                   character()) 
  expect_identical(m2$getParents("e1", omit = c("c1", "c2"), upstream = TRUE, stochOnly = FALSE ),
                   "d1") 
  m2$setData(list(d1 = 5))
  expect_identical(m2$getParents("g1", upstream = TRUE, includeData = FALSE, stochOnly = TRUE ),
                   c("b1", "e1"))
})

test_that("getParents works in model with LHSinferred (aka split) nodes", {
  c3 <- nimbleCode({
    a[1:3] ~ dmnorm( mu[1:3], cov = cov[1:3, 1:3])
    b[1:3] <- a[1] + c(1, 2, 3)
    sig <- sqrt(var)
    var ~ dunif(0, 1)
    c[1] ~ dnorm(b[1], sd = sig)
    a2 ~ dnorm(0,1)
    sig2 ~ dunif(0,1)
    c[2] ~ dnorm(a[2], sd = sig2)
  })
  m3 <- nimbleModel(c3, inits = list(mu = 1:3, cov = diag(3)))

  expect_identical(m3$getParents("c[1]", stochOnly = TRUE), c("var", "a[1:3]"))
  expect_identical(m3$getParents("c[2]", stochOnly = TRUE), c("sig2", "a[1:3]"))
  expect_identical(m3$getParents("c[1:2]", stochOnly = TRUE), c("var", "sig2", "a[1:3]"))
  expect_identical(m3$getParents("c[1:2]", immediateOnly = TRUE, stochOnly = TRUE), c("sig2", "a[1:3]"))
  expect_identical(m3$getParents("sig", stochOnly = TRUE), c("var"))
  expect_identical(m3$getDependencies("a[1]"), c("a[1:3]", "b[1:3]", "c[1]"))
  expect_identical(m3$getDependencies("a[2]"), c("a[1:3]", "c[2]"))
  expect_identical(m3$getDependencies("b[2]"), c("b[1:3]"))
  expect_identical(m3$getDependencies("b[1]"), c("b[1:3]", "c[1]"))
  expect_identical(m3$getDependencies("b[1]", returnScalarComponents = TRUE), c("b[1]", "b[2]", "b[3]", "c[1]"))
})

## Test getConditionallyIndependentSets
test_that("getConditionallyIndependentSets works in model with a couple of sets", {
  mc <- nimbleCode({
    mu ~ dnorm(0,1)
    for(i in 1:2) {
      x[i] ~ dnorm(mu, 1)
      y[i] ~ dnorm(x[i], 1)
      z[i] ~ dnorm(y[i], 1)
    }
  })
  m <- nimbleModel(mc, data = list(z = 1:2))

  expect_identical(m$getConditionallyIndependentSets(), list(c('x[1]', 'y[1]'), c('x[2]', 'y[2]')))
  expect_identical(m$getConditionallyIndependentSets('y[2]', inputType = "param"), list('y[2]'))
  expect_identical(m$getConditionallyIndependentSets('x[1:2]', inputType = "data"), list(c('x[1]'), c('x[2]')))
  expect_true(nimble:::testConditionallyIndependentSets(m, m$getConditionallyIndependentSets()))
  expect_identical(m$getConditionallyIndependentSets(omit = 'y[2]'), list(c('x[1]', 'y[1]'), c('x[2]')))
  expect_identical(m$getConditionallyIndependentSets(omit = 5), list(c('x[1]', 'y[1]'), c('x[2]')))
})

test_that("getConditionallyIndependentSets works in model with a couple of sets", {
  # Two state-space models, one of which has data at the end.
  mc <- nimbleCode({
    x[1] ~ dnorm(0, 1)
    w[1] ~ dnorm(0, 1)
    for(i in 1:3) {
      x[i+1] ~ dnorm(x[i], 1)
      y[i+1] ~ dnorm(x[i+1], 1)
      w[i+1] ~ dnorm(w[i], 1)
      z[i+1] ~ dnorm(w[i+1], 1)    
    }
  })
  m <- nimbleModel(mc, data = list(y = 1:4))
  expect_identical(m$getConditionallyIndependentSets(), list(c("x[2]", "x[3]", "x[4]"),
                                                             c("w[2]", "w[3]", "z[2]", "w[4]", "z[3]", "z[4]")))
  expect_identical(m$getConditionallyIndependentSets("x[1:2]", givenNodes = c("w[1]", "y")), list(c("x[1]", "x[2]", "x[3]", "x[4]")))
  expect_identical(m$getConditionallyIndependentSets("x[2]", givenNodes = c("w[1]", "y")), list(c("x[1]", "x[2]", "x[3]", "x[4]")))
  expect_identical(m$getConditionallyIndependentSets("x[1]", givenNodes = c("w[1]", "y")), list(c("x[1]", "x[2]", "x[3]", "x[4]")))
  expect_identical(m$getConditionallyIndependentSets(givenNodes = c("y")),
                   list(c("x[1]", "x[2]", "x[3]", "x[4]"),
                        c("w[1]", "w[2]", "w[3]", "z[2]", "w[4]", "z[3]", "z[4]")))
  expect_identical(m$getConditionallyIndependentSets('z[1]'), list())
  expect_identical(m$getConditionallyIndependentSets('y[2]', inputType = "param"), list())
  expect_true(nimble:::testConditionallyIndependentSets(m, m$getConditionallyIndependentSets()))
})

test_that("getConditionallyIndependentSets works in model with diamond shape", {
  # diamond graph shape
  mc <- nimbleCode({
    mu ~ dnorm(0,1)
    for(i in 1:2) {
      x[i] ~ dnorm(mu, 1)
    }
    y ~ dnorm(x[1] + x[2], 1)
  })
  m <- nimbleModel(mc, data = list(y = 1))

  expect_identical(m$getConditionallyIndependentSets(), list(c("x[1]","x[2]")))
  nimble:::testConditionallyIndependentSets(m, m$getConditionallyIndependentSets())
})

test_that("getConditionallyIndependentSets works in double-state state-space model", {
  # two stae-space chains of latent states with one data set that depends on both
  mc <- nimbleCode({
    x[1] ~ dnorm(0, 1)
    w[1] ~ dnorm(0, 1)
    for(i in 1:3) {
      x[i+1] ~ dnorm(x[i], 1)
      w[i+1] ~ dnorm(w[i], 1)
      y[i+1] ~ dnorm(x[i+1] + w[i+1], 1)
    }
  })
  m <- nimbleModel(mc, data = list(y = 1:4))

  expect_identical(m$getConditionallyIndependentSets(),
                   list(c("x[2]", "w[2]", "x[3]", "w[3]", "x[4]", "w[4]")))
  expect_identical(m$getConditionallyIndependentSets(omit = "w[2]"),
                   list(c("x[2]", "x[3]", "w[3]", "x[4]", "w[4]")))
  expect_identical(m$getConditionallyIndependentSets(givenNodes = c("y", "w[3]")),
                   list(c("x[1]", "w[1]", "x[2]", "w[2]", "x[3]", "x[4]", "w[4]")))
  expect_identical(m$getConditionallyIndependentSets(givenNodes = c("y", "x[3]", "w[3]")),
                   list(c("x[1]", "w[1]", "x[2]", "w[2]"), c("x[4]", "w[4]")))
  expect_true(nimble:::testConditionallyIndependentSets(m, m$getConditionallyIndependentSets()))
})


test_that("getConditionallyIndependentSets works in model with LHSinferred (aka split) nodes", {
  c3 <- nimbleCode({
    a[1:3] ~ dmnorm( mu[1:3], cov = cov[1:3, 1:3])
    b[1:3] <- a[1] + c(1, 2, 3)
    sig <- sqrt(var)
    var ~ dunif(0, 1)
    c[1] ~ dnorm(b[1], sd = sig)
    a2 ~ dnorm(0,1) ## isolated node
    sig2 ~ dunif(0,1)
    c[2] ~ dnorm(a[2], sd = sig2)
  })
  m3 <- nimbleModel(c3, inits = list(mu = 1:3, cov = diag(3)))

  expect_identical(m3$getConditionallyIndependentSets(), list()) # Note there are no latent nodes
  expect_identical(m3$getConditionallyIndependentSets("a[1:3]", givenNodes = "sig2"),
                   list(c("var", "a[1:3]", "c[2]", "c[1]")))
  expect_identical(m3$getConditionallyIndependentSets("a[1]", givenNodes = "sig2"),
                   list(c("var", "a[1:3]", "c[2]", "c[1]")))
  expect_identical(m3$getConditionallyIndependentSets("a[2]", givenNodes = "sig2"),
                   list(c("var", "a[1:3]", "c[2]", "c[1]")))
  expect_identical(m3$getConditionallyIndependentSets("b[1]", givenNodes = "sig2"),
                   list(c("var", "a[1:3]", "c[2]", "c[1]")))
  expect_identical(m3$getConditionallyIndependentSets(m3$getNodeNames(stochOnly = TRUE), givenNodes = "sig2"),
                   list(c("var", "a[1:3]" , "c[2]", "c[1]"), c("a2")))
})

test_that("getConditionallyIndependentSets works for tweaked pump model", {
  pumpCode <- nimbleCode({ 
    # This pump code is tweaked so that theta[1] and theta[2] are
    # in a conditionally independent set together.
    # Other theta[i]s are in their own set.
    for (i in 3:N){
      theta[i] ~ dgamma(alpha, beta) 
      lambda[i] <- theta[i] * t[i]
      x[i] ~ dpois(lambda[i]) 
    }
    for (i in 1:2){
      theta[i] ~ dgamma(alpha, beta) 
      lambda[i] <- 0.5 * (theta[1] * t[1] + theta[2] * t[2])
      x[i] ~ dpois(lambda[i]) 
    }    
    alpha ~ dexp(1.0) 
    beta ~ dgamma(0.1, 1.0) 
  })

  pumpConsts <- list(N = 10,
                     t = c(94.3, 15.7, 62.9, 126, 5.24,
                           31.4, 1.05, 1.05, 2.1, 10.5))
  
  pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))
  
  pumpInits <- list(alpha = 0.1, beta = 0.1,
                    theta = rep(0.1, pumpConsts$N))
  
  
  ## Create the model
  pump <- nimbleModel(code = pumpCode, name = "pump", constants = pumpConsts,
                      data = pumpData, inits = pumpInits)
  
  expect_identical(pump$getConditionallyIndependentSets(),
                   list('theta[3]', 'theta[4]', 'theta[5]', 'theta[6]',
                        'theta[7]', 'theta[8]', 'theta[9]', 'theta[10]',
                        c('theta[1]', 'theta[2]')))  
})
    
options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
