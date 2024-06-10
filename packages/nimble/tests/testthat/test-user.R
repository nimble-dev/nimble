source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = TRUE)

goldFileName <- 'userTestLog_Correct.Rout'
tempFileName <- 'userTestLog.Rout'
generatingGoldFile <- !is.null(nimbleOptions('generateGoldFileForUserTesting'))
outputFile <- if(generatingGoldFile) file.path(nimbleOptions('generateGoldFileForUserTesting'), goldFileName) else tempFileName

sink(outputFile)

nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)

## Note on scoping and where objects will be found in test_that evaluation.
## registerDistributions, deregisterDistributions, and nimbleModel
## all now grab the calling environment and use it sensibly.
## Therefore there is no reason to use temporarilyAssignInGlobalEnv
##   if the test is limited to uncompiled and does not use NF objects.
## compileNimble(model) will need to find user-defined functions in
##   .GlobalEnv, so the temporarilyAssignInGlobalEnv trick will be needed.
## Also it appears that if a nimbleFunction object (with `$`, as of 1.2.0)
##   is used to provide a user-defined distribution or function, then
## nimbleModel will need the object to be in .GlobalEnv.
test_that("User-supplied functions", {
    dbl <- nimbleFunction(
        run = function(x = double(0)) {
            returnType(double(0))
            return(2*x)
        }
    )
    temporarilyAssignInGlobalEnv(dbl)
    ## not working at moment as enclosing env of fun is getting messed up


    ## two arguments to the nimbleFunction
    dblSum <- nimbleFunction(
        run = function(x = double(0), y = double(0)) {
            returnType(double(0))
            return(2*(x+y))
        }
    )
    temporarilyAssignInGlobalEnv(dblSum)

    ## vector input and output
    vecdbl <- nimbleFunction(
        run = function(x = double(1)) {
            returnType(double(1))
            return(2*x)
        }
    )
    temporarilyAssignInGlobalEnv(vecdbl)

                                        # function that will allow testing of arg matching by name
    mypow <- nimbleFunction(
        run = function(x = double(0), y = double(0)) {
            returnType(double(0))
            return(pow(x, y))
        }
    )
    temporarilyAssignInGlobalEnv(mypow)

    code <- nimbleCode({
        x ~ dnorm(0, 1)
        dx ~ dnorm(dbl(x), sd = .01)
        liftedmu[1:K] <- vecdbl(mu[1:K])
        y[1:K] ~ dmnorm(liftedmu[1:K], cov = eps[1:K,1:K])
        mu[1:K] ~ dmnorm(zeros[1:K], cov = I[1:K, 1:K])
        z ~ dnorm(0, 1)
        dz ~ dnorm(dblSum(x, z), sd = .01)

                                        # use of args out of order
        out <- mypow(y = yin, x = xin)
        
                                        # vectorized fun applied to scalar nodes-based variable
        for(i in 1:K) {
            theta[i] ~ dnorm(0, 1)
        }
        liftedmu2[1:K] <- vecdbl(theta[1:K])
        w[1:K] ~ dmnorm(liftedmu2[1:K], cov = eps[1:K, 1:K])
    })

    K <- 3
    m <- nimbleModel(code, inits = list(yin = 3, xin = 2, x = 0.25, y = 1:K, mu = 1:K,
                                        z = 0.5, theta = rep(.5, K), w = rep(1, K)),
                     constants = list(K = K, zeros = rep(0, K), I = diag(K),
                                      eps = diag(rep(.0001, K))))

    cm <- compileNimble(m)

    set.seed(0)
    simulate(m)
    set.seed(0)
    simulate(cm)

                                        # need c() in here because R nodes are arrays
    for(var in c('dx', 'y', 'dz', 'w')) {
        expect_equal(c(get(var, m)), (get(var, cm)),
                     info = paste0("Test that R and C models agree with user-supplied functions:", var, " values differ"))
    }
    expect_lt(abs(2*cm$x - cm$dx), (.03)) # "Test that values based on user-supplied functions are correct (x and dx consistent)"
    expect_lt(max(abs(2*cm$mu - cm$y)), (.03)) # "Test that values based on user-supplied functions are correct (mu and y consistent)"
    expect_lt(abs(2*(cm$x + cm$z) - cm$dz), (.03)) # "Test that values based on user-supplied functions are correct (x plus z and dz consistent)"
    expect_lt(max(abs(2*cm$theta - cm$w)), (.03)) # "Test that values based on user-supplied functions are correct (theta and w consistent)"
    expect_identical(c(m$out), (8),
                     info = "incorrect arg matching by name in R model")
    expect_identical(cm$out, (8),
                     info = "incorrect arg matching by name in C model")
})

## User-supplied distributions

test_that("Use of user-supplied distributions and registerDistributions", {

    dmyexp <- nimbleFunction(
        run = function(x = double(0), rate = double(0, default = 1), log = integer(0, default = 0)) {
            returnType(double(0))
            logProb <- log(rate) - x*rate
            if(log) {
                return(logProb)
            } else {
                return(exp(logProb))
            }
        })
    temporarilyAssignInGlobalEnv(dmyexp)


    rmyexp <- nimbleFunction(
        run = function(n = integer(0), rate = double(0, default = 1)) {
            returnType(double(0))
            if(n != 1) nimPrint("rmyexp only allows n = 1; using n = 1.")
            dev <- runif(1, 0, 1)
            return(-log(1-dev) / rate)
        }
    )
    temporarilyAssignInGlobalEnv(rmyexp)

    pmyexp <- nimbleFunction(
        run = function(q = double(0), rate = double(0, default = 1), lower.tail = integer(0, default = 1), log.p = integer(0, default = 0)) {
            returnType(double(0))
            if(!lower.tail) {
                logp = -rate * q
                if(log.p) {
                    return(logp)
                } else {
                    return(exp(logp))
                }
            } else {
                p = 1 - exp(-rate * q)
                if(!log.p) {
                    return(p)
                } else {
                    return(log(p))
                }
            }
        }
    )
    temporarilyAssignInGlobalEnv(pmyexp)

    qmyexp <- nimbleFunction(
        run = function(p = double(0), rate = double(0, default = 1), lower.tail = integer(0, default = 1), log.p = integer(0, default = 0)) {
            returnType(double(0))
            if(log.p) {
                p = exp(p)
            }
            if(!lower.tail) {
                p = 1 - p
            }
            return(-log(1 - p) / rate)
        }
    )
    temporarilyAssignInGlobalEnv(qmyexp)


    ddirchmulti <- nimbleFunction(
        run = function(x = double(1), alpha = double(1), size = double(0), log = integer(0)) {
            returnType(double(0))
            logProb <- lgamma(size) - sum(lgamma(x)) + lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum(lgamma(alpha + x)) - lgamma(sum(alpha) + size)
            if(log) {
                return(logProb)
            } else {
                return(exp(logProb))
            }
        })
    temporarilyAssignInGlobalEnv(ddirchmulti)

    rdirchmulti <- nimbleFunction(
        run = function(n = integer(0), alpha = double(1), size = double(0)) {
            returnType(double(1))
            if(n != 1) nimPrint("rdirchmulti only allows n = 1; using n = 1.")
            p <- rdirch(1, alpha)
            return(rmulti(1, size = size, prob = p))
        })
    temporarilyAssignInGlobalEnv(rdirchmulti)

    expect_warning(
        registerDistributions(list(
            dmyexp = list(
                BUGSdist = "dmyexp(rate, scale)",
                Rdist = "dmyexp(rate = 1/scale)",
                altParams = "scale = 1/rate",
                pqAvail = TRUE),
            ddirchmulti = list(
                BUGSdist = "ddirchmulti(alpha, size)",
                types = c('x = double(1)', 'alpha = double(1)', 'size = double(0)'))
        )), "Found 'x' in 'types'")

    code1 <- nimbleCode({
        for(i in 1:n1) {
            y1[i] ~ dmyexp(rate = r1)
            y2[i] ~ dmyexp(scale = s2)
        }
        r1 <- 1 / s1
        s1 ~ dunif(0, 100)
        s2 ~ dunif(0, 100)
    })

    code2 <- nimbleCode({
        for(i in 1:n2) {
            y3[i] ~ dpois(lambda)
        }
        lambda ~ T(dmyexp(scale = 5), 0, upper)
    })

    code3 <- nimbleCode({
        for(i in 1:m) {
            y[i, 1:P] ~ ddirchmulti(alpha[1:P], sz)
        }
        for(i in 1:P) {
            alpha[i] ~ dunif(0, 1000) # dgamma(.001, .001);
        }
    })

    set.seed(0)
    mn = 3
    n1 <- 1000
    y1 <- rexp(n1, rate = 1/mn)
    y2 <- rexp(n1, rate = 1/mn)
    
    lambda = 2.5
    n2 <- 50
    y3 <- rpois(n2, lambda)
    
    upper <-  3
    
    sz <- 100
    alpha <- 10*c(1,2,3)
    m <- 40
    P <- length(alpha)
    y <- p <- matrix(0, nrow = m, ncol = P)
    for( i in 1:m ) {
        p[i, ] <- rdirch(1, alpha)
        y[i, ] <- rmultinom(1, size = sz, prob  = p[i,])
    }
    
    data1 <- list(y1 = y1, y2 = y2, n1 = n1)
    data2 <- list(y3 = y3, n2 = n2, upper = upper)
    data3 <- list(y = y, m = m, P = P, sz = sz)
    
    inits1 <- list(s1 = 1, s2 = 1)
    inits2 <- list(lambda = 1)
    inits3 <- list(alpha = rep(30, P))
    
    m2 <- nimbleModel(code2, data = data2['y3'], constants = data2[c('upper', 'n2')],
                      inits = inits2)
    cm2 <- compileNimble(m2)
    simulate(m2, 'lambda')
    simulate(cm2, 'lambda')
    expect_lt(max(c(m2$lambda, cm2$lambda)), (upper)) # "parameter not exceed upper bound"
    m2$lambda <- cm2$lambda <- upper + 1
    expect_identical(max(c(calculate(m2, 'lambda'), calculate(cm2, 'lambda'))), (-Inf), info = "calculation on out-of-bounds not return -Inf")

    testBUGSmodel(code1, example = 'user1', dir = "", data = data1, inits = inits1, useInits = TRUE)
    testBUGSmodel(code2, example = 'user2', dir = "", data = data2, inits = inits2, useInits = TRUE)
    testBUGSmodel(code3, example = 'user3', dir = "", data = data3, inits = inits3, useInits = TRUE)

    test_mcmc(model = code1, data = data1, inits = inits1,
              results = list(mean = list(s1 = mn, s2 = mn)),
              resultsTolerance = list(mean = list(s1 = .2, s2 = .2)))
    
    test_mcmc(model = code3, data = data3, inits = inits3,
              results = list(mean = list(alpha = alpha)),
              resultsTolerance = list(mean = list(alpha = c(4, 6, 8))),
              numItsC_results = 50000)

    m <- nimbleModel(code2, constants = data2, inits = inits2)
    cm <- compileNimble(m)
    
    spec <- configureMCMC(m)
    Rmcmc <- buildMCMC(spec)
    Cmcmc <- compileNimble(Rmcmc, project = m)
    
    Cmcmc$run(5000)
    smp <- as.matrix(Cmcmc$mvSamples)
    
    expect_lt(max(smp[ , 'lambda']), (upper))
})


test_that("Test that user-defined distributions with 'lower' and 'upper' params don't cause problems", {
    dmyunif <- nimbleFunction(
    run = function(x = double(0), lower = double(0), upper = double(0), log = integer(0, default = 0)) {
        if(x > lower & x < upper)
            return(1) else return(0)
        returnType(double())
    })
    rmyunif <- nimbleFunction(
        run = function(n = double(0), lower = double(0), upper = double(0)) {
            return(0)
            returnType(double(0))
        })
    temporarilyAssignInGlobalEnv(dmyunif)
    temporarilyAssignInGlobalEnv(rmyunif)

    code = nimbleCode({
        y ~ dmyunif(lower = 3, upper = 7)
    })
    m <- nimbleModel(code, inits = list(y=4))
    cm <- compileNimble(m)
    expect_equal(m$getParam('y','lower'), 3)
    expect_equal(cm$getParam('y','lower'), 3)
    expect_equal(m$getBound('y','lower'), -Inf)
    expect_equal(cm$getBound('y','lower'), -Inf)
})



test_that("Test that deregistration of user-supplied distributions works", {
    deregisterDistributions('ddirchmulti')
    expect_true(is.null(nimble:::nimbleUserNamespace$distributions[['ddirchmulti']]),
                info = "ddirchmulti has not been deregistered")
})


# Test user-defined without r function: note registerDistributions()
# puts r function in calling frame.
test_that("User-defined dist without 'r' works when registered by nimbleModel", {
  dfooABC = nimbleFunction(
    run = function(x = double(1), log = integer(0, default = 0)) {
      returnType(double(0))
      return(0)
    })
  expect_true(exists('dfooABC', inherits=FALSE))

  code <- nimbleCode({
      v[1:3] ~ dfooABC()
  })

  m <- nimbleModel(code)
  ## temporarilyAssignInGlobalEnv(dfooABC)
  ## temporarilyAssignInGlobalEnv(rfooABC)
  ## cm <- compileNimble(m)
  expect_true(exists('rfooABC', inherits=FALSE))
  deregisterDistributions('dfooABC')
  expect_false(exists('rfooABC'))
  expect_false(exists('rfooABC'))
})

## Test user-defined without r function, with use of registerDistributions()
## Here we explicitly put it in .GlobalEnv
test_that("Test that user-defined distributions without 'r' function doesn't cause problems", {
    ## scalar density
    dfooDEF = nimbleFunction(
        run = function(x = double(0), log = integer(0, default = 0)) {
            returnType(double(0))
            return(0)
        })
    
    temporarilyAssignInGlobalEnv(dfooDEF)

    registerDistributions(list(
        dfooDEF = list(
            BUGSdist = "dfooDEF()")), userEnv = .GlobalEnv)
    
    code <- nimbleCode({
        v ~ dfooDEF()
    })
    ## Can't use expect_message as message doesn't appear when run Travis for some reason.
    m <- nimbleModel(code)
    cm <- compileNimble(m)
    expect_output(print(m), "Rmodel object")

    ## non-scalar density
    dfooDEF2 = nimbleFunction(
        run = function(x = double(1), log = integer(0, default = 0)) {
            returnType(double(0))
            return(0)
        })
    
    temporarilyAssignInGlobalEnv(dfooDEF2)

    registerDistributions(list(
        dfooDEF2 = list(
            BUGSdist = "dfooDEF2()",
            types = c('value = double(1)'))), userEnv = .GlobalEnv)
    
    code <- nimbleCode({
        v[1:3] ~ dfooDEF2()
    })
    m <- nimbleModel(code)
    cm <- compileNimble(m)
    expect_output(print(m), "Rmodel object")
    expect_true(exists("dfooDEF", envir=.GlobalEnv))
    expect_true(exists("rfooDEF", envir=.GlobalEnv))
    expect_true(exists("dfooDEF2", envir=.GlobalEnv))
    expect_true(exists("rfooDEF2", envir=.GlobalEnv))
    expect_silent(fooinfo <- getDistributionInfo("dfooDEF"))
    expect_silent(fooinfo <- getDistributionInfo("dfooDEF2"))
    deregisterDistributions("dfooDEF", .GlobalEnv)
    deregisterDistributions("dfooDEF2", .GlobalEnv)
    expect_error(fooinfo <- getDistributionInfo("dfooDEF"))
    expect_error(fooinfo <- getDistributionInfo("dfooDEF2"))
    expect_true(exists("dfooDEF", envir=.GlobalEnv))
    expect_false(exists("rfooDEF", envir=.GlobalEnv))
    expect_true(exists("dfooDEF2", envir=.GlobalEnv))
    expect_false(exists("rfooDEF2", envir=.GlobalEnv))
})

test_that("Test that non-scalar integer in user-defined distributions is trapped", {
    dfoo3 = nimbleFunction(run = function(x = integer(1), log=logical(0)) {
        returnType(double(0))
        return(0)
    })
    rfoo3 = nimbleFunction(run = function(n = integer()) {
        returnType(integer(1))
        return(rep(0,2))
    })
#    temporarilyAssignInGlobalEnv(dfoo3)
#    temporarilyAssignInGlobalEnv(rfoo3)
    
    code <- nimbleCode({
        y[1:3] ~ dfoo3()
    })
    expect_error(m <- nimbleModel(code), "Non-scalar integer or logical found")
    
    dfoo4 = nimbleFunction(run = function(x = double(1), theta = integer(1), log=logical(0)) {
        tmp <- sum(theta)
        returnType(double(0))
        return(0)
    })
    rfoo4 = nimbleFunction(run = function(n = integer(), theta = integer(1)) {
        returnType(double(1))
        return(rep(0,2))
    })
#    temporarilyAssignInGlobalEnv(dfoo4)
#    temporarilyAssignInGlobalEnv(rfoo4)

    code <- nimbleCode({
        y[1:3] ~ dfoo4(theta[1:3])
    })
    expect_error(m <- nimbleModel(code), "Non-scalar integer or logical found")

    dfoo5 = nimbleFunction(run = function(x = integer(0), theta = integer(0), log=logical(0)) {
        returnType(double(0))
        return(0)
    })
    rfoo5 = nimbleFunction(run = function(n = integer(), theta = integer(0)) {
        returnType(integer())
        return(0)
    })
#    temporarilyAssignInGlobalEnv(dfoo5)
#    temporarilyAssignInGlobalEnv(rfoo5)
    
    code <- nimbleCode({
        y  ~ dfoo5(theta)
    })
    m <- nimbleModel(code)
    expect_silent(dfoo5info <- getDistributionInfo("dfoo5"))
    deregisterDistributions("dfoo5")
    expect_error(dfoo5info <- getDistributionInfo("dfoo5"))
})

test_that("Test that missing/mismatched returnType in 'r' function is trapped", {
  dDist <- nimbleFunction(
        run = function(x = double(1), log = integer(0, default = 0)) {
            returnType(double())
            return(0)
        }
    )
    
    rDist <- nimbleFunction(
        run = function(n = integer()) {
            ##returnType(double(1))   ## FAILED TO PROVIDE
            return(rep(1,2))
        })
    
    code <- nimbleCode({ x[1,1:2] ~ dDist()})
    
    expect_error(Rmodel <- nimbleModel(code), "missing `returnType`")
    
    dDist <- nimbleFunction(
        run = function(x = double(1), log = integer(0, default = 0)) {
            returnType(double())
            return(0)
        }
    )
    
    rDist <- nimbleFunction(
        run = function(n = integer()) {
            returnType(double(1))
            return(rep(1,2))
        })
    
    code <- nimbleCode({ x[1,1:2] ~ dDist()})
    
    Rmodel <- nimbleModel(code)
    expect_true(exists('Rmodel', inherits=FALSE))
  expect_silent(stuff <- getDistributionInfo("dDist"))
  deregisterDistributions('dDist')
  expect_error(stuff <- getDistributionInfo("dDist"))
    rm(Rmodel)
    
    dDist <- nimbleFunction(
        run = function(x = double(1), log = integer(0, default = 0)) {
            returnType(double())
            return(0)
        }
    )
    
    rDist <- nimbleFunction(
        run = function(n = integer()) {
            returnType(double(0))
            return(rep(1,2))
        })
    
    code <- nimbleCode({ x[1,1:2] ~ dDist()})
    
    expect_error(Rmodel <- nimbleModel(code), "missing `returnType`")

    ## double() vs. double(0) should be fine
    dDist <- nimbleFunction(
        run = function(x = double(0), log = integer(0, default = 0)) {
            returnType(double())
            return(0)
        }
    )
    
    rDist <- nimbleFunction(
        run = function(n = integer()) {
            returnType(double())
            return(1)
        })
    
    code <- nimbleCode({ x ~ dDist()})
    
    Rmodel <- nimbleModel(code)
    expect_true(exists('Rmodel', inherits=FALSE))
    expect_silent(stuff <- getDistributionInfo("dDist"))
    deregisterDistributions('dDist', .GlobalEnv)
    expect_error(stuff <- getDistributionInfo("dDist"))
})

test_that("default argument for `x` in density is trapped", {
  ddist <- nimbleFunction(
    run = function(x = double(0, default = 0), log = integer(default = 1)) {
        returnType(double())
        return(0)
    }
  )
  expect_error(registerDistributions('ddist'), "not allowed to have a default")
})

test_that("Trap case where simulate function is removed", {
    dfooGHI <- nimbleFunction(
        run = function(x=double(), log=integer(0, default=0)) {
            return(0)
            returnType(double())
        }
    )
    temporarilyAssignInGlobalEnv(dfooGHI, replace=TRUE)

    m <- nimbleModel(
        nimbleCode({
            x ~ dfooGHI()
        })
    )
    temporarilyAssignInGlobalEnv(rfooGHI, replace=TRUE)

    cm <- compileNimble(m)

    rm('rfooGHI', pos = .GlobalEnv)
    dfooGHI <- nimbleFunction(
        run = function(x=double(), log=integer(0, default=0)) {
            return(0)
            returnType(double())
        }
    )
    temporarilyAssignInGlobalEnv(dfooGHI, replace=TRUE)

    m <- nimbleModel(
        nimbleCode({
            x ~ dfooGHI()
        })
    )

    expect_error(cm <- compileNimble(m), "is not available")
    expect_silent(stuff <- getDistributionInfo("dfooGHI"))
    deregisterDistributions('dfooGHI')
    expect_error(stuff <- getDistributionInfo("dfooGHI"))
})

test_that("Trap case where simulate function should be removed", {
    dfooJKL <- nimbleFunction(
        run = function(x = double(), p = double(), log = integer(default = 0)) {
            ans <- x * log(p)
            return(ans)
            returnType(double())
        }
    )
    m <- nimbleModel(quote({
        param <- 1.2
        x ~ dfooJKL(p = param)
    }))

    expect_true(exists('rfooJKL', inherits=FALSE))
    expect_silent(stuff <- getDistributionInfo("dfooJKL"))
    deregisterDistributions('dfooJKL')
    expect_error(stuff <- getDistributionInfo("dfooJKL"))
    expect_false(exists('rfooJKL', inherits=FALSE))
})

sink(NULL)

if(!generatingGoldFile) {
    test_that("Log file matches gold file", {
        trialResults <- readLines(tempFileName)
        correctResults <- readLines(system.file(file.path('tests', 'testthat', goldFileName), package = 'nimble'))
        compareFilesByLine(trialResults, correctResults)
    })
}

###############################################
## Tests for nimbleFunctions with setup code ##
###############################################

currentOption <- getNimbleOption("allowNFobjInModel")
nimbleOptions(allowNFobjInModel = TRUE)

test_that("user-defined dists with setup code with d only",{
  myDist <- nimbleFunction(
    setup = function() {mean <- 1},
    methods = list(
      dmynorm = function(x = double(0), sd = double(), log = integer(0, default = 0)) {
        return(dnorm(x, mean = mean, sd = sd, log = log))
        returnType(double())
      })
  )

  myDist1 <- myDist()
  temporarilyAssignInGlobalEnv(myDist1)

  mc <- nimbleCode({
    x ~ myDist1$dmynorm(2)
  })

  m <- nimbleModel(mc, data = list(y = 3), inits = list(x = .5))
  temporarilyAssignInGlobalEnv(rmynorm_myDist1)
  cm <- compileNimble(m)
  expect_equal(m$calculate(), cm$calculate())
  expect_silent(stuff <- getDistributionInfo("myDist1$dmynorm"))
  expect_true(exists("rmynorm_myDist1", inherits=FALSE))
  deregisterDistributions('myDist1$dmynorm')
  expect_false(exists("rmynorm_myDist1", inherits=FALSE))
  expect_error(stuff <- getDistributionInfo("myDist1$dmynorm"))
})

test_that("user-defined dists with setup code errors out for bad name",{
  myDist <- nimbleFunction(
    setup = function() {mean <- 1},
    methods = list(
      # missing "d" prefix
      mynorm = function(x = double(0), sd = double(), log = integer(0, default = 0)) {
        return(dnorm(x, mean = mean, sd = sd, log = log))
        returnType(double())
      })
  )

  myDist1 <- myDist()
  temporarilyAssignInGlobalEnv(myDist1)

  mc <- nimbleCode({
    x ~ myDist1$mynorm(2)
  })

  expect_error(
    m <- nimbleModel(mc, data = list(y = 3), inits = list(x = .5)),
    "must begin with 'd'")
  expect_warning(deregisterDistributions('myDist1$dmynorm'))
})

test_that("user-defined dists with setup code: trap non-existing object",{
  mc <- nimbleCode({
    x ~ myGhost1$dmynorm(2)
  })
  expect_error(
    m <- nimbleModel(mc, debug = FALSE, data = list(y = 3), inits = list(x = .5)),
    "could not find"
  )
})

test_that("user-defined dists with setup code: trap non-existing method",{
  myGhost <- nimbleFunction(
    setup = function() {mean <- 1},
    methods = list(
      dmynorm = function(x = double(0), sd = double(), log = integer(0, default = 0)) {
        return(dnorm(x, mean = mean, sd = sd, log = log))
        returnType(double())
      })
  )

  myGhost1 <- myGhost()
  temporarilyAssignInGlobalEnv(myGhost1)

  mc <- nimbleCode({
    x ~ myGhost1$dmyexp(2) # wrong name, needs to be trapped
  })

  expect_error(
    m <- nimbleModel(mc, debug = FALSE, data = list(y = 3), inits = list(x = .5)),
    "is not a method"
  )
  expect_error(stuff <- getDistributionInfo("myGhost1$dmyexp"))
})

test_that("user-defined dists with setup code with d and r",{
  myDistDEF <- nimbleFunction(
    setup = function() {mean <- 1},
    methods = list(
      dmynorm = function(x = double(0), sd = double(), log = integer(0, default = 0)) {
        return(dnorm(x, mean = mean, sd = sd, log = log))
        returnType(double())
      },
      rmynorm = function(n = integer(), sd = double()) {
        return(rnorm(1, mean = mean, sd = sd))
        returnType(double())
      })
  )

  myDistDEF1 <- myDistDEF()
  temporarilyAssignInGlobalEnv(myDistDEF1)

  mc <- nimbleCode({
    x ~ myDistDEF1$dmynorm(2)
  })

  m <- nimbleModel(mc, debug = FALSE, data = list(y = 3), inits = list(x = .5))
  set.seed(1)
  m$simulate()
  expect_false(exists("rmynorm_myDistDEF1", inherits=FALSE))

  cm <- compileNimble(m)
  set.seed(1)
  cm$simulate()
  expect_equal(m$x, cm$x)
  expect_silent(stuff <- getDistributionInfo("myDistDEF1$dmynorm"))
  deregisterDistributions('myDistDEF1$dmynorm')
  expect_error(stuff <- getDistributionInfo("myDistDEF1$dmynorm"))
})

test_that("user-defined dists with setup code: trap error of simulating with dummy r",{
  myDistGHI <- nimbleFunction(
    setup = function() {mean <- 1},
    methods = list(
      dmynorm = function(x = double(0), sd = double(), log = integer(0, default = 0)) {
        return(dnorm(x, mean = mean, sd = sd, log = log))
        returnType(double())
      })
  )

  myDistGHI1 <- myDistGHI()
  temporarilyAssignInGlobalEnv(myDistGHI1)

  mc <- nimbleCode({
    x ~ myDistGHI1$dmynorm(2)
  })

  m <- nimbleModel(mc, debug = FALSE, data = list(y = 3), inits = list(x = .5))
  temporarilyAssignInGlobalEnv(rmynorm_myDistGHI1)
  expect_error(m$simulate(), "provided without random generation")

  cm <- compileNimble(m)
  expect_error(cm$simulate(), "provided without random generation")
  expect_equal(m$calculate(), cm$calculate())
  expect_silent(stuff <- getDistributionInfo("myDistGHI1$dmynorm"))
  deregisterDistributions('myDistGHI1$dmynorm')
  expect_error(stuff <- getDistributionInfo("myDistGHI1$dmynorm"))
})

test_that("deregister user-defined dists with setup code",{
  myDistJKL <- nimbleFunction(
    setup = function() {mean <- 1},
    methods = list(
      dmynorm = function(x = double(0), sd = double(), log = integer(0, default = 0)) {
        return(dnorm(x, mean = mean, sd = sd, log = log))
        returnType(double())
      })
  )

  myDistJKL1 <- myDistJKL()
  temporarilyAssignInGlobalEnv(myDistJKL1)

  mc <- nimbleCode({
    x ~ myDistJKL1$dmynorm(2)
  })

  m <- nimbleModel(mc, debug = FALSE, data = list(y = 3), inits = list(x = .5))

  expect_true(exists("rmynorm_myDistJKL1", inherits = FALSE))
  expect_silent(stuff <- getDistributionInfo("myDistJKL1$dmynorm"))
  deregisterDistributions('myDistJKL1$dmynorm')
  expect_error(stuff <- getDistributionInfo("myDistJKL1$dmynorm"))
  expect_false(exists("rmynorm_myDistJKL1", inherits = FALSE))
})

test_that("user-defined dists with setup code with reduced name",{
  myDistMNO <- nimbleFunction(
    setup = function() {mean <- 1},
    methods = list(
      d = function(x = double(0), sd = double(), log = integer(0, default = 0)) {
        return(dnorm(x, mean = mean, sd = sd, log = log))
        returnType(double())
      })
  )

  myDistMNO1 <- myDistMNO()
  temporarilyAssignInGlobalEnv(myDistMNO1)

  mc <- nimbleCode({
    x ~ myDistMNO1$d(2)
  })

  m <- nimbleModel(mc, debug = FALSE, data = list(y = 3), inits = list(x = .5))
  temporarilyAssignInGlobalEnv(r_myDistMNO1)

  cm <- compileNimble(m)
  expect_equal(m$calculate(), cm$calculate())
  expect_silent(stuff <- getDistributionInfo("myDistMNO1$d"))
  deregisterDistributions('myDistMNO1$d')
  expect_error(stuff <- getDistributionInfo("myDistMNO1$d"))
})

test_that("user-defined dists with setup code with r and multiple dists",{
  myDistPQR <- nimbleFunction(
    setup = function() {mean <- 1; lambda <- 2},
    methods = list(
      d1 = function(x = double(0), sd = double(), log = integer(0, default = 0)) {
        return(dnorm(x, mean = mean, sd = sd, log = log))
        returnType(double())
      },
      d2 = function(x = double(0), log = integer(0, default = 0)) {
        return(dpois(x, lambda, log = log))
        returnType(double())
      },
      r1 = function(n = integer(), sd = double()) {
        return(rnorm(1, mean = mean, sd = sd))
        returnType(double())
      })
  )

  myDistPQR1 <- myDistPQR()
  temporarilyAssignInGlobalEnv(myDistPQR1)

  mc <- nimbleCode({
    x ~ myDistPQR1$d1(2)
    y ~ myDistPQR1$d2()
  })

  m <- nimbleModel(mc, debug = FALSE, data = list(y = 3), inits = list(x = .5, y = 2))
  temporarilyAssignInGlobalEnv(r2_myDistPQR1)

  cm <- compileNimble(m)
  expect_equal(m$calculate(), cm$calculate())
  expect_silent(stuff <- getDistributionInfo("myDistPQR1$d1"))
  expect_silent(stuff <- getDistributionInfo("myDistPQR1$d2"))
  deregisterDistributions('myDistPQR1$d1')
  deregisterDistributions('myDistPQR1$d2')
  expect_error(stuff <- getDistributionInfo("myDistPQR1$d1"))
  expect_error(stuff <- getDistributionInfo("myDistPQR1$d2"))
})

test_that("user-defined functions with setup code and dist in same NF",{
  fooXYZ <- nimbleFunction(
    setup = function() {f <- 2},
    run = function(z = double(1)) {
      return(f*z)
      returnType(double(1))
    },
    methods = list(
      go = function() {
        return(f)
        returnType(double())
      },
      dbar = function(x = double(), mean = double(), log = integer(0, default=0)) {
        return(dnorm(x, mean = mean, sd = f, log = log))
        returnType(double())
      }
    )
  )

  fooXYZ1 <- fooXYZ()
  temporarilyAssignInGlobalEnv(fooXYZ1)

  mc <- nimbleCode({
    for(i in 1:2) z[i] ~ dnorm(0, sd = 1)
    y[1:2] <- fooXYZ1$run(z[1:2])
    x[1:2] <- exp(fooXYZ1$run(y[1:2])) # check that it works nexted
    w ~ dnorm(sum(y[1:2]), sd = 1)
    q <- fooXYZ1$go()
    r ~ fooXYZ1$dbar(q)
  })

  m <- nimbleModel(mc, debug = FALSE, data = list(y = c(1.5, 2.5)), inits = list(z = c(.5, .7)))
  temporarilyAssignInGlobalEnv(rbar_fooXYZ1)

  cm <- compileNimble(m)
  expect_equal(m$calculate(), cm$calculate())
  expect_silent(stuff <- getDistributionInfo("fooXYZ1$dbar"))
  deregisterDistributions('fooXYZ1$dbar')
  expect_error(stuff <- getDistributionInfo("fooXYZ1$dbar"))
})

test_that("user-defined functions: trap non-existent object",{
  mc <- nimbleCode({
    q <- Ghost2$go()
  })

  expect_error(
    m <- nimbleModel(mc, debug = FALSE, data = list(y = c(1.5, 2.5)), inits = list(z = c(.5, .7))),
    "not found"
  )
})

test_that("user-defined functions: trap (message only) non-existent method",{
  fooUVW <- nimbleFunction(
    setup = function() {f <- 2},
    methods = list(
      go = function() {
        return(f)
        returnType(double())
      }
    )
  )

  fooUVW1 <- fooUVW()
  temporarilyAssignInGlobalEnv(fooUVW1) # not sure why needed here

  mc <- nimbleCode({
    q <- fooUVW1$ghost()
  })

  expect_message(
    m <- nimbleModel(mc, debug = FALSE, data = list(y = c(1.5, 2.5)), inits = list(z = c(.5, .7))),
    "not a valid field or method"
  )
})

test_that("user-defined distributions: p and q found and used", {
  myDistUVW <- nimbleFunction(
    setup = function() {mean <- 1},
    methods = list(
      dmynorm = function(x = double(0), sd = double(), log = integer(0, default = 0)) {
        return(dnorm(x, mean = mean, sd = sd, log = log))
        returnType(double())
      },
      qmynorm = function(p = double(), sd = double(), lower.tail = integer(0, default=1), log.p = integer(0, default = 0)) {
        return(qnorm(p, mean = mean, sd = sd, lower.tail = lower.tail, log.p = log.p))
        returnType(double())
      },
      pmynorm = function(q = double(), sd = double(), lower.tail = integer(0, default=1), log.p = integer(0, default = 0)) {
        return(pnorm(q, mean = mean, sd = sd, lower.tail = lower.tail, log.p = log.p))
        returnType(double())
      }
    )
  )

  myDistUVW1 <- myDistUVW()
  temporarilyAssignInGlobalEnv(myDistUVW1)

  mc <- nimbleCode({
    x ~ T(myDistUVW1$dmynorm(2), 0, Inf)
#    x ~ T(dnorm(2, sd = 2), 0, Inf)
  })

  m <- nimbleModel(mc, debug = FALSE, data = list(y = 3), inits = list(x = .5))
  cm <- compileNimble(m)
  expect_true(exists('cm',inherits=FALSE))
  expect_silent(stuff <- getDistributionInfo("myDistUVW1$dmynorm"))
  deregisterDistributions("myDistUVW1$dmynorm")
  expect_error(stuff <- getDistributionInfo("myDistUVW1$dmynorm"))
})

test_that("user-defined in NF object example from User Manual works", {
  code <- nimbleCode({
    beta0 ~ dnorm(0, sd = 100)
    beta1 ~ dnorm(0, sd = 100)
    sigma ~ dhalfflat()
    for(i in 1:n) {
      y[i] ~ dnorm(beta0 + beta1*x[i], sd = sigma)
    }
  })
  # Simulate some data:
  set.seed(1)
  x <- runif(10)
  y <- rnorm(10, 3 + 0.5 * x, sd = 0.2)
  model <- nimbleModel(code, constants = list(n = 10, x = x),
                     data = list(y = y),
                     inits = list(beta0=3.1,beta1=0.4,sigma=0.3))
  llh1 <- model$calculate()

  linear_reg <- nimbleFunction(
    setup = function(predictor) {},
    methods = list(
      dcalc = function(x = double(1),
                       beta0=double(), beta1=double(), sigma=double(),
                       log=integer(0, default=0)) {
        llh <- sum(dnorm(x, beta0+beta1*predictor, sd=sigma, log=TRUE))
        if(log) return(llh)
        return(exp(llh))
        returnType(double())
      }))
  my_linear_reg <- linear_reg(x)
  code2 <- nimbleCode({
    beta0 ~ dnorm(0, sd = 100)
    beta1 ~ dnorm(0, sd = 100)
    sigma ~ dhalfflat()
    y[1:n] ~ my_linear_reg$dcalc(beta0, beta1, sigma)
  })
  temporarilyAssignInGlobalEnv(my_linear_reg)
  model <- nimbleModel(code2, constants = list(n = 10),
                       data = list(y = y),
                       inits = list(beta0=3.1,beta1=0.4,sigma=0.3))
  temporarilyAssignInGlobalEnv(rcalc_my_linear_reg)
  llh2 <- model$calculate()
  expect_equal(llh1, llh2)
  deregisterDistributions("my_linear_reg$dcalc")
})

nimbleOptions(allowNFobjInModel = currentOption)

nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)
