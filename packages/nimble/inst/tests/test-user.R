source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of user-supplied distributions and functions in BUGS code")

nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

goldFileName <- 'userTestLog_Correct.Rout'
tempFileName <- 'userTestLog.Rout'
generatingGoldFile <- !is.null(nimbleOptions('generateGoldFileForUserTesting'))
outputFile <- if(generatingGoldFile) file.path(nimbleOptions('generateGoldFileForUserTesting'), goldFileName) else tempFileName

sink(outputFile)

nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)

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

# would like to test user-defined without r function, but registerDistributions()
# puts r function in calling frame and that is causing problems when called from
# within testing context

## dfoo = nimbleFunction(
##     run = function(x = double(1), log = integer(0, default = 0)) {
##         returnType(double(0))
##         return(0)
##     })

## temporarilyAssignInGlobalEnv(dfoo)

## code <- nimbleCode({
##     v[1:3] ~ dfoo()
## })

## out <- try(m <- nimbleModel(code))

## try(test_that("Test of user-supplied distribution without r function: ",
##               expect_false(is(out, 'try-error'))))

## Instead, test user-defined without r function, but with use of registerDistributions()

test_that("Test that user-defined distributions without 'r' function doesn't cause problems", {
    ## scalar density
    dfoo = nimbleFunction(
        run = function(x = double(0), log = integer(0, default = 0)) {
            returnType(double(0))
            return(0)
        })
    
    temporarilyAssignInGlobalEnv(dfoo)

    registerDistributions(list(
        dfoo = list(
            BUGSdist = "dfoo()")), userEnv = .GlobalEnv)
    
    code <- nimbleCode({
        v ~ dfoo()
    })
    m <- nimbleModel(code)
    cm <- compileNimble(m)

    ## non-scalar density
    dfoo2 = nimbleFunction(
        run = function(x = double(1), log = integer(0, default = 0)) {
            returnType(double(0))
            return(0)
        })
    
    temporarilyAssignInGlobalEnv(dfoo2)

    registerDistributions(list(
        dfoo2 = list(
            BUGSdist = "dfoo2()",
            types = c('value = double(1)'))), userEnv = .GlobalEnv)
    
    code <- nimbleCode({
        v[1:3] ~ dfoo2()
    })
    m <- nimbleModel(code)
    cm <- compileNimble(m)
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
    temporarilyAssignInGlobalEnv(dfoo3)
    temporarilyAssignInGlobalEnv(rfoo3)
    
    code =nimbleCode({
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
    temporarilyAssignInGlobalEnv(dfoo4)
    temporarilyAssignInGlobalEnv(rfoo4)

    code =nimbleCode({
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
    temporarilyAssignInGlobalEnv(dfoo5)
    temporarilyAssignInGlobalEnv(rfoo5)
    
    code =nimbleCode({
        y  ~ dfoo5(theta)
    })
    m <- nimbleModel(code)
})

sink(NULL)

if(!generatingGoldFile) {
    trialResults <- readLines(tempFileName)
    correctResults <- readLines(system.file(file.path('tests', goldFileName), package = 'nimble'))
    compareFilesByLine(trialResults, correctResults)
}

nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)
