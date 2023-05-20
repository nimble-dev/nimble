# This file tests both setupMargNodes and model$getConditionallyIndependentSets
# There is some overlap with test-ADlaplace, which relies on these features.

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

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
  expect_identical(m$getConditionallyIndependentSets('y[2]', explore = "down"), list('y[2]'))
  expect_identical(m$getConditionallyIndependentSets('x[1:2]', explore = "up"), list(c('x[1]'), c('x[2]')))
  expect_true(nimble:::testConditionallyIndependentSets(m, m$getConditionallyIndependentSets()))
  expect_identical(m$getConditionallyIndependentSets(omit = 'y[2]'), list(c('x[1]', 'y[1]'), c('x[2]')))
  expect_identical(m$getConditionallyIndependentSets(omit = 5), list(c('x[1]', 'y[1]'), c('x[2]')))
  expect_identical(m$getConditionallyIndependentSets('x[1]'), list(c('x[1]')))
  expect_identical(m$getConditionallyIndependentSets('x[1]', unknownAsGiven=FALSE), list(c('x[1]', 'y[1]')))

  SMN <- setupMargNodes(m)
  expect_identical(SMN$paramNodes, "mu")
  expect_identical(SMN$randomEffectsNodes, c("x[1]", "x[2]", "y[1]", "y[2]"))

  SMN <- setupMargNodes(m, randomEffectsNodes = c("y[1]"))
  expect_identical(SMN$paramNodes, "x[1]")
  expect_identical(SMN$calcNodes, c("y[1]", "z[1]"))

  SMN <- setupMargNodes(m, paramNodes = character())
  expect_identical(SMN$randomEffectsNodes, character())
  expect_identical(SMN$calcNodes, character())

  SMN <- setupMargNodes(m, calcNodes = c("y[1]", "z[1]"))
  expect_identical(SMN$paramNodes, "x[1]")
  expect_identical(SMN$randomEffectsNodes, "y[1]")

  SMN <- setupMargNodes(m, calcNodes = c("x[1]", "y[1]", "z[1]"))
  expect_identical(SMN$paramNodes, "mu")
  expect_identical(SMN$randomEffectsNodes, c("x[1]","y[1]"))

  SMN <- setupMargNodes(m, randomEffectsNodes = "mu")
  expect_identical(SMN$paramNodes, character())
  expect_identical(SMN$calcNodes, c("mu","x[1]","x[2]"))

  expect_warning(SMN <- setupMargNodes(m, paramNodes = character(),
                                       randomEffectsNodes = "x[1]", calcNodes = c("x[1]", "y[1]")))
  expect_warning(SMN <- setupMargNodes(m, paramNodes = c("mu"),
                                       randomEffectsNodes = "y[1]", calcNodes = c("y[1]", "z[1]")))

})

test_that("setupMargNodes/GCIS works in model with an extra edge among random effects", {
  mc <- nimbleCode({
    mu ~ dnorm(0,1)
    for(i in 1:2) {
      x[i] ~ dnorm(mu, 1)
      y[i] ~ dnorm(x[i], 1)
      z[i] ~ dnorm(y[i] + x[i], 1)
    }
  })
  m <- nimbleModel(mc, data = list(z = 1:2))
  SMN <- setupMargNodes(m)
  expect_identical(SMN$paramNodes, "mu")
  expect_identical(SMN$randomEffectsNodes, c("x[1]", "x[2]", "y[1]", "y[2]"))

  expect_warning(SMN <- setupMargNodes(m, randomEffectsNodes = "z[1]"))
  expect_warning(SMN <- setupMargNodes(m, randomEffectsNodes = "x[1]", calcNodes = c("y[1]", "z[1]")))
  expect_warning(SMN <- setupMargNodes(m, randomEffectsNodes = "x[1]", calcNodes = c("x[1]","y[1]", "z[1]")))
  expect_warning(SMN <- setupMargNodes(m, paramNodes = "mu", randomEffectsNodes = "y[1]"))
  SMN <- setupMargNodes(m, paramNodes = "x[1]", randomEffectsNodes = "y[1]")
  expect_identical(SMN$calcNodes,c('y[1]','lifted_y_oBi_cB_plus_x_oBi_cB_L5[1]','z[1]'))
})

test_that("setupMargNodes/GCIS works in model with an extra edge from param to random effects", {
  mc <- nimbleCode({
    mu ~ dnorm(0,1)
    for(i in 1:2) {
      x[i] ~ dnorm(mu, 1)
      y[i] ~ dnorm(x[i], 1)
      z[i] ~ dnorm(y[i] + mu, 1)
    }
  })
  m <- nimbleModel(mc, data = list(z = 1:2))
  SMN <- setupMargNodes(m)
  expect_identical(SMN$paramNodes, "mu")
  expect_identical(SMN$randomEffectsNodes, c("x[1]", "x[2]", "y[1]", "y[2]"))

  SMN <- setupMargNodes(m, randomEffectsNodes = c("y[1]", "x[2]"))
  expect_identical(SMN$paramNodes, c("mu", "x[1]"))

  # This gives a spurious warning that we live with
  SMN <- setupMargNodes(m, randomEffectsNodes = c("mu", "x[1]"))
  expect_identical(SMN$paramNodes, character())
  expect_identical(SMN$givenNodes, c('x[2]','y[1]','z[1]','z[2]'))
})

test_that("setupMargNodes/GCIS catches discrete randomEffectsNode", {
  code <- nimbleCode({
    p ~ dnorm(0,1)
    re1 ~ dnorm(p, 1)
    re2 ~ dbern(re1)    ## discrete node here
    re3 ~ dnorm(re2, 1)
    y ~ dnorm(re3, 1)
  })
  Rmodel <- nimbleModel(code, data = list(y = 1))

  expect_warning(SMN <- setupMargNodes(Rmodel, 'p'))
  expect_identical(SMN$randomEffectsNodes, c('re1','re3'))
})

test_that("getConditionallyIndependentSets works in model with one set and deterministic intermediates", {
  m <- nimbleModel(
    nimbleCode({
      P1 ~ dnorm(0,1)
      D1 <- P1 + 1

      REA1 ~ dnorm(D1, 1)
      D2 <- REA1 + 1
      REA2 ~ dnorm(D2, 1)
      D3 <- REA2 + 1

      REB1 ~ dnorm(D1, 1)
      C2 <- REB1 + 1
      REB2 ~ dnorm(C2, 1)
      C3 <- REB2 + 1

      D3C3 <- D3 + C3
      Y1 ~ dnorm(D3C3, 1)
    }),
    data = list(Y1 = 1)
  )
  # All sets
  expect_identical(m$getConditionallyIndependentSets(), list(c("REA1", "REB1", "REA2", "REB2")))
  # first-stage latents stay separated
  expect_identical(m$getConditionallyIndependentSets(c("REA1", "REB1")), list("REA1", "REB1"))
  # first-stage latents are combined if unknownAsGiven is FALSE
  expect_identical(m$getConditionallyIndependentSets(c("REA1", "REB1"), unknownAsGiven=FALSE), list(c("REA1", "REB1", "REA2", "REB2")))
  # second-stage latents are connected
  expect_identical(m$getConditionallyIndependentSets(c("REA2", "REB2")), list(c("REA2", "REB2")))
  # Using determinstics as given works
  expect_identical(m$getConditionallyIndependentSets(givenNodes = c("D1", "D3C3")), list(c("REA1", "REB1", "REA2", "REB2")))
  # Using determinstics as given works
  expect_identical(m$getConditionallyIndependentSets(givenNodes = c("D1", "D3", "C3")), list(c("REA1", "REA2"), c("REB1", "REB2")))
  expect_identical(m$getConditionallyIndependentSets(givenNodes = c("D1", "D3")),
                   list(c("REA1", "REA2")))
  expect_identical(m$getConditionallyIndependentSets(givenNodes = c("D3")),
                   list(c("REA1", "REA2")))
  expect_identical(m$getConditionallyIndependentSets(givenNodes = c("D3"), unknownAsGiven=FALSE),
                   list(c("P1", "REA1", "REB1", "REA2", "REB2", "Y1")))

  SMN <- setupMargNodes(m)
  expect_identical(SMN$paramNodes, "P1")
  expect_identical(SMN$randomEffectsNodes, c('REA1','REB1','REA2','REB2'))
  expect_identical(SMN$randomEffectsSets, list(c('REA1','REB1','REA2','REB2')))

  SMN <- setupMargNodes(m, paramNodes = "REA1")
  expect_identical(SMN$randomEffectsNodes, c('REA2'))
  expect_identical(SMN$calcNodes, c('REA2','D3','D3C3','Y1'))

  SMN <- setupMargNodes(m, calcNodes = m$getDependencies("REB2"))
  expect_identical(SMN$paramNodes, "REB1")

  SMN <- setupMargNodes(m, randomEffectsNodes = c("P1", "REA1","REA2", "REB1","REB2"))
  expect_identical(SMN$paramNodes, character())
  expect_identical(SMN$calcNodes, c('P1','D1','REA1','REB1','D2','C2','REA2','REB2','D3','C3','D3C3','Y1'))
}

test_that("getConditionallyIndependentSets works in state-space model with a couple of sets", {
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
  expect_identical(getConditionallyIndependentSets(m, endAsGiven=TRUE), list(c("x[2]", "x[3]", "x[4]"),
                                                                             c("w[2]", "w[3]", "w[4]")))
  expect_identical(m$getConditionallyIndependentSets(), list(c("x[2]", "x[3]", "x[4]")))
  expect_identical(m$getConditionallyIndependentSets("x[1:2]", givenNodes = c("w[1]", "y"), unknownAsGiven=FALSE),
                   list(c("x[1]", "x[2]", "x[3]", "x[4]")))
  expect_identical(m$getConditionallyIndependentSets("x[2]", givenNodes = c("w[1]", "y"), unknownAsGiven=FALSE),
                   list(c("x[1]", "x[2]", "x[3]", "x[4]")))
  expect_identical(m$getConditionallyIndependentSets("x[1]", givenNodes = c("w[1]", "y"), unknownAsGiven=FALSE),
                   list(c("x[1]", "x[2]", "x[3]", "x[4]")))
  expect_identical(m$getConditionallyIndependentSets(givenNodes = c("y"), unknownAsGiven=FALSE),
                   list(c("x[1]", "x[2]", "x[3]", "x[4]")))
  expect_identical(m$getConditionallyIndependentSets('z[1]'), list())
  expect_true(nimble:::testConditionallyIndependentSets(m, m$getConditionallyIndependentSets()))

  SMN <- setupMargNodes(m)
  expect_identical(SMN$paramNodes, "x[1]")
  expect_identical(SMN$randomEffectsNodes, c("x[2]","x[3]","x[4]"))

  SMN <- setupMargNodes(m, randomEffectsNodes = 'x[1:4]')
  expect_identical(SMN$paramNodes, character())
  expect_identical(SMN$randomEffectsNodes, c("x[1]","x[2]","x[3]","x[4]"))
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
  expect_true(nimble:::testConditionallyIndependentSets(m, m$getConditionallyIndependentSets()))

  SMN <- setupMargNodes(m)
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
  expect_identical(m$getConditionallyIndependentSets(givenNodes = c("y", "w[3]"), unknownAsGiven=FALSE),
                   list(c("x[1]", "w[1]", "x[2]", "w[2]", "x[3]", "x[4]", "w[4]")))
  expect_identical(m$getConditionallyIndependentSets(givenNodes = c("y", "x[3]", "w[3]"), unknownAsGiven=FALSE),
                   list(c("x[1]", "w[1]", "x[2]", "w[2]"), c("x[4]", "w[4]")))
  expect_true(nimble:::testConditionallyIndependentSets(m, m$getConditionallyIndependentSets()))

  SMN <- setupMargNodes(m)
  expect_identical(SMN$randomEffectsSets,
                   list(c('x[2]','w[2]','x[3]','w[3]','x[4]','w[4]')))
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
  expect_identical(m3$getConditionallyIndependentSets("a[1:3]", givenNodes = "sig2", unknownAsGiven=FALSE),
                   list(c("var", "a[1:3]", "c[2]", "c[1]")))
  expect_identical(m3$getConditionallyIndependentSets("a[1]", givenNodes = "sig2", unknownAsGiven=FALSE),
                   list(c("var", "a[1:3]", "c[2]", "c[1]")))
  expect_identical(m3$getConditionallyIndependentSets("a[2]", givenNodes = "sig2", unknownAsGiven=FALSE),
                   list(c("var", "a[1:3]", "c[2]", "c[1]")))
# Revisit
  #  expect_identical(m3$getConditionallyIndependentSets("b[1]", givenNodes = "sig2", unknownAsGiven=FALSE),
#                   list(c("var", "a[1:3]", "c[2]", "c[1]")))
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

test_that("getConditionallyIndependentSets works with unknownAsGiven=TRUE or FALSE", {
  # This test is from NCT #405
  m <- nimbleModel(
  nimbleCode({
    mu ~ dnorm(0,1)
    for(i in 1:4) a[i] ~ dnorm(mu, 1)
    y[1] ~ dnorm(a[1]+a[2], 1)
    y[2] ~ dnorm(a[1]-a[2], 1)
    y[3] ~ dnorm(a[3]+a[4], 1)
    y[4] ~ dnorm(a[3]-a[4], 1)
  }),
  data = list(y = c(1, 1, 1, 1)),
  inits = list(mu = 1))

  expect_identical(
    m$getConditionallyIndependentSets("a", givenNodes = c("y"), unknownAsGiven=FALSE),
    list(m$expandNodeNames(c("mu","a"))))

  expect_identical(
    m$getConditionallyIndependentSets("a", givenNodes = c("y"), unknownAsGiven=TRUE),
    list(c('a[1]','a[2]'), c('a[3]','a[4]')))
})
