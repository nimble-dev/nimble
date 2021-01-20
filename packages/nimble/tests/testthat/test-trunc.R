source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

context("Testing of truncation, censoring, and constraints")

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)
oldWidth <- getOption("width")
options(width = 1000)
oldMaxPrint <- getOption("max.print")
options(max.print = 1000)

goldFileName <- 'truncTestLog_Correct.Rout'
tempFileName <- 'truncTestLog.Rout'
generatingGoldFile <- !is.null(nimbleOptions('generateGoldFileForTruncTesting'))
outputFile <- if(generatingGoldFile)
                  file.path(nimbleOptions('generateGoldFileForTruncTesting'), goldFileName) else tempFileName

sink(outputFile)

nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)


# should build testing of kidney,litters,mice,lsat into test-mcmc and test-models

models <- c('kidney', 'litters', 'mice')

### test basic model building
out <- sapply(models, testBUGSmodel, useInits = TRUE)

# we haven't updated conjugacy to deal with truncation
# so litters is reporting conj post density as wrong since
# otherwise p[i,j] is a standard conjugate update

# lsat has non-explicit indexing
## changes for windows compatibility
file.copy( system.file('classic-bugs','vol1','lsat','lsat.bug', package = 'nimble'), file.path(tempdir(), "lsat.bug"))
##system(paste("cat", system.file('classic-bugs','vol1','lsat','lsat.bug', package = 'nimble'), ">>", file.path(tempdir(), "lsat.bug")))
system.in.dir("sed -i -e 's/mean(alpha\\[\\])/mean(alpha\\[1:T\\])/g' lsat.bug", dir = tempdir())
##system(paste("sed -i -e 's/mean(alpha\\[\\])/mean(alpha\\[1:T\\])/g'", file.path(tempdir(), "lsat.bug")))
system.in.dir("sed -i -e 's/p.item\\[i,\\]/p.item\\[i,1:T\\]/g' lsat.bug", dir = tempdir())
##system(paste("sed -i -e 's/p.item\\[i,\\]/p.item\\[i,1:T\\]/g'", file.path(tempdir(), "lsat.bug")))
testBUGSmodel('lsat', dir = "", model = file.path(tempdir(), "lsat.bug"), data = system.file('classic-bugs','vol1','lsat','lsat-data.R', package = 'nimble'),  inits = system.file('classic-bugs','vol1','lsat','lsat-init.R', package = 'nimble'),  useInits = TRUE)


### test MCMC
out <- sapply(models, test_mcmc, numItsC = 1000)

# this takes forever because of MCMC - rewrite to not do basic
if(FALSE) {
    test_mcmc(model = file.path(tempdir(), "lsat.bug"), inits = system.file('classic-bugs', 'vol1', 'lsat','lsat-init.R', package = 'nimble'), data = system.file('classic-bugs', 'vol1', 'lsat','lsat-data.R', package = 'nimble'), numItsC = 1000)
}

### test truncation

# test variations of T(), I() syntax

filename <- file.path(tempdir(), "mod.bug")

# basic T() test
writeLines(c("model {",
             "x ~ dnorm(mu, sd = 1)",
             "mu ~ dnorm(0, sd = 2) T(  0,5) }"),
           con = filename)
inits <- list(mu = 2)
data <- list(x = 1)
testBUGSmodel(model = filename, dir = '',
              data = data, inits = inits)

# left-truncted T() test
writeLines(c("model {",
             "x ~ dnorm(mu, sd = 1)",
             "mu ~ dnorm(0, sd = 2) T(-3, ) }"),
           con = filename)
inits <- list(mu = 2)
data <- list(x = -2)
testBUGSmodel(model = filename, dir = '',
              data = data, inits = inits)

# non-truncated T() test
writeLines(c("model {",
             "x ~ dnorm(mu, sd = 1)",
             "mu ~ dnorm(0, sd = 2) T(,) }"),
           con = filename)
inits <- list(mu = 2)
data <- list(x = -2)
testBUGSmodel(model = filename, dir = '',
              data = data, inits = inits)

# basic I() test
writeLines(c("model {",
             "x ~ dnorm(mu, sd = 1)",
             "mu ~ dnorm(0, sd = 2) I(,3) }"),
           con = filename)
inits <- list(mu = 2)
data <- list(x = 2)
testBUGSmodel(model = filename, dir = '',
              data = data, inits = inits, expectModelWarning = TRUE)

# test of T() with stochastic truncation and expressions
writeLines(c("model {",
             "x ~ dnorm(mu, sd = 1)",
             "b ~ dunif(0, 2)",
             "mu ~ dnorm(0, sd = 3) T( , exp(b + 1.2)) }"),
           con = filename)
inits <- list(mu = 2, b = 1)
data <- list(x = 2)
testBUGSmodel(model = filename, dir = '',
              data = data, inits = inits)

# test basic MCMC
test_that("Test that MCMC respects truncation bounds", {
    set.seed(0)
    
    code <- nimbleCode({
        for(i in 1:n) {
            y[i] ~ dnorm(mu, 1)
        }
        mu ~ T(dunif(-10, 10), 0.2, 2)
    })
    n <- 10
    mu <- 1.2
    constants <- list(n = n)
    data <- list(y = rnorm(n, mu, 1))
    inits <- list(mu = 1.5)
    
    Rmodel = nimbleModel(code, data = data, inits = inits,
                         constants = constants)
    mcmcConf <- configureMCMC(Rmodel)
    Rmcmc <- buildMCMC(mcmcConf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    Cmcmc$run(5000)
    smp <- as.matrix(Cmcmc$mvSamples)
    expect_gt(min(smp), 0.2, label = "minimum in MCMC", expected.label = "lower bound")
    expect_lt(max(smp), 2, label = "maximum in MCMC", expected.label = "upper bound")

    test_mcmc(model = code, data = c(data, constants), inits = inits,
              results = list(mean = list(mu = 1.5), sd = list(mu = .27 )),
              resultsTolerance = list(mean = list(mu = 0.05), sd = list(mu = .05)),
              name = 'test of MCMC with truncation')
})

### censoring (dinterval)
test_that("Test that MCMC respects right-censoring bounds", {
    set.seed(0)

    code <- nimbleCode({
        for(i in 1:n) {
            y[i] ~ dnorm(mu, sd = 10)
            is.cens[i] ~ dinterval(y[i], bnd)
        }
        mu ~ dunif(30, 60)
        bnd ~ dunif(45, 55)
    })
    n <- 20
    mu <- 45
    bnd <- 50
    y <- rnorm(n, mu, 10)
    is.cens <- y > bnd
    y[is.cens] <- NA
    y[1] <- NA
    is.cens[1] <- TRUE
    yinits <- y
    yinits[is.cens] <- 55
    yinits[!is.cens] <- NA

    constants <- list(n = n)
    data <- list(y = y, is.cens = is.cens)
    inits <- list(mu = mu, bnd = bnd, y = yinits)

    Rmodel <- nimbleModel(code, data = data, inits = inits,
                         constants = constants)
    mcmcConf <- configureMCMC(Rmodel)
    mcmcConf$addMonitors('y[1]')
    Rmcmc <- buildMCMC(mcmcConf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    Cmcmc$run(5000)
    smp <- as.matrix(Cmcmc$mvSamples)
    expect_gt(min(smp[ , 'bnd'] - max(y[!is.cens])), 0, "upper bound in MCMC exceeds observations")
    expect_gt(min(smp[ , 'y[1]'] - smp[ , 'bnd']), 0, "minimum inferred censored value in MCMC greater than bound")

    test_mcmc(model = code, data = c(data, constants), inits = inits, numItsC = 10000,
              results = list(mean = list(mu = 44.5, 'y[1]' = 56.6), sd = list(mu = 2.3, 'y[1]' = 4.95)),
              resultsTolerance = list(mean = list(mu = 0.1, 'y[1]' = 1.5), sd = list(mu = .3, 'y[1]' = .7)), name = 'test of right censoring')
})

test_that("Test that MCMC respects left-censoring bounds", {
    set.seed(0)
    code <- nimbleCode({
        for(i in 1:n) {
            y[i] ~ dnorm(mu, sd = 10)
            not.cens[i] ~ dinterval(y[i], bnd)
        }
        mu ~ dunif(30, 60)
        bnd ~ dunif(35, 45)
    })
    n <- 20
    mu <- 45
    bnd <- 40
    y <- rnorm(n, mu, 10)
    not.cens <- y > bnd
    y[!not.cens] <- NA
    y[1] <- NA
    not.cens[1] <- FALSE
    yinits <- y
    yinits[!not.cens] <- 35
    yinits[not.cens] <- NA
    
    constants <- list(n = n)
    data <- list(y = y, not.cens = not.cens)
    inits <- list(mu = mu, bnd = bnd, y = yinits)
    
    Rmodel <- nimbleModel(code, data = data, inits = inits,
                         constants = constants)
    mcmcConf <- configureMCMC(Rmodel)
    mcmcConf$addMonitors('y[1]')
    Rmcmc <- buildMCMC(mcmcConf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    Cmcmc$run(5000)
    smp <- as.matrix(Cmcmc$mvSamples)
    expect_gt(min(min(y[not.cens]) - smp[ , 'bnd']), 0, label = "lower bound in MCMC",
              expected.label = "observations")
    expect_gt(max(smp[ , 'bnd'] - smp[ , 'y[1]']), 0, label = "maximum inferred censored value in MCMC", expected.label = "bound")

    test_mcmc(model = code, data = c(data, constants), inits = inits,
              results = list(mean = list(mu = 43.4, 'y[1]' = 33), sd = list(mu = 2.4, 'y[1]' = 6.0)),
              resultsTolerance = list(mean = list(mu = .2, 'y[1]' = .5), sd = list(mu = .2, 'y[1]' = .1)), name = 'test of left censoring')
})

# interval censored
test_that("Test that MCMC respects interval censoring", {
    code <- nimbleCode({
        for(i in 1:(3*n)) {
            y[i] ~ dnorm(mu, sd = 10)
            intvl[i] ~ dinterval(y[i], bnd[1:2])
        }
        mu ~ dunif(30, 60)
    })
    n <- 5
    mu <- 45
    bnd <- c(40, 50)
    intvl <- c(rep(0, n), rep(1, n), rep(2, n))
    y <- rep(NA, 3*n)
    y[1] <- 38
    y[n+1] <- 42
    y[2*n+1] <- 53

    constants <- list(n = n, bnd = bnd)
    data <- list(y = y, intvl = intvl)
    yinits <- y
    yinits[is.na(y)] <- c(rep(35, 4), rep(45, 4), rep(55, 4))
    yinits[!is.na(y)]  <- NA
    inits <- list(mu = mu, y = yinits)

    Rmodel = nimbleModel(code, data = data, inits = inits, constants = constants)
    mcmcConf <- configureMCMC(Rmodel)
    mcmcConf$addMonitors(c('y[2]', 'y[7]', 'y[12]'))
    Rmcmc <- buildMCMC(mcmcConf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    Cmcmc$run(5000)
    smp <- as.matrix(Cmcmc$mvSamples)
    expect_lt(max(smp[ , 'y[1]']), bnd[1], label = "imputed observation in first interval")
    expect_gt(min(smp[ , 'y[7]']), bnd[1], label = "imputed observation in second interval")
    expect_lt(max(smp[ , 'y[7]']), bnd[2], label = "imputed observation in second interval")
    expect_gt(min(smp[ , 'y[12]']), bnd[2], label = "imputed observation in third interval")

    test_mcmc(model = code, data = c(data, constants), inits = inits,
              results = list(mean = list(mu = 44.77, 'y[12]' = 56.3), sd = list(mu = 2.8, 'y[12]' = 5.3)),
              resultsTolerance = list(mean = list(mu = 0.5, 'y[12]' = 0.5), sd = list(mu = .5, 'y[12]' = .5)), name = 'test of interval censoring')
})

# test of dconstraint
test_that("Test that MCMC respects constraints", {
    set.seed(0)
    code <- nimbleCode ({
        for(i in 1:n) {
            y1[i] ~ dnorm(mu1, 1)
            y2[i] ~ dnorm(mu2, 1)
        }
        mu1 ~ dunif(-10, 10)
        mu2 ~ dunif(-10, 10)
        mu.c ~ dconstraint(mu1 + mu2 > 0 & mu1 > 1)
    })
    mu1 <- 1; mu2 <- -0.5
    n <- 10
    y1 <- rnorm(n, mu1, 1)
    y2 <- rnorm(n, mu2, 1)
    data = list(y1 = y1, y2 = y2, mu.c =1)
    inits = list(mu1 = 1.5, mu2 = 0.5)
    constants = list(n = n)

    Rmodel = nimbleModel(code, data = data, inits = inits, constants = constants)
    mcmcConf <- configureMCMC(Rmodel)
    Rmcmc <- buildMCMC(mcmcConf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    Cmcmc$run(5000)
    smp <- as.matrix(Cmcmc$mvSamples)
    expect_gt(min(smp[ , 'mu1'] + smp[ , 'mu2']), 0, label = "constraint on sum of mu1 and mu2")
    expect_gt(min(smp[ , 'mu1']), 0, label = "constraint on mu1")

    test_mcmc(model = code, data = c(data, constants), inits = inits,
              results = list(mean = list(mu1 = 1.45, mu2 = -.82), sd = list(mu1 = .26, mu2 = .30)),
              resultsTolerance = list(mean = list(mu1 = 0.05, mu2 = .02), sd = list(mu1 = .03, mu2 = .03)), name = 'test of dconstraint')
})

# test dconstraint for ordering using rewrite of inhaler in our syntax

# pasting BUGS code here because inhaler has a data black that has non-standard
# use of ":" in for loop indexing, and so we can substitute our use of
# dconstraint in place of JAGS sort()

test_that("Test that MCMC respects constraints in inhaler example", {
    code <- nimbleCode({
        for (i in 1:N) {
            for (t in 1:T) {
                for (j in 1:Ncut) {
                    ##
                    ## Cumulative probability of worse response than j
                    ##
                    logit(Q[i,t,j]) <- -(a0[j] + mu[group[i],t] + b[i]);
                }
                ##
                ## Probability of response = j
                ##
                p[i,t,1] <- 1 - Q[i,t,1];
                for (j in 2:Ncut) { p[i,t,j] <- Q[i,t,j-1] - Q[i,t,j] }
                p[i,t,(Ncut+1)] <- Q[i,t,Ncut];

                response[i,t] ~ dcat(p[i,t,1:4]);
            }
            ##
            ## Subject (random) effects
            ##
            b[i] ~ dnorm(0.0, tau);
        }

        ##
        ## Fixed effects
        ##
        for (g in 1:G) {
            for(t in 1:T) {
                                        # logistic mean for group i in period t
                mu[g,t] <- beta*treat[g,t]/2 + pi*period[g,t]/2 + kappa*carry[g,t];
            }
        }
        beta ~ dnorm(0, 1.0E-06);
        pi ~ dnorm(0, 1.0E-06);
        kappa ~ dnorm(0, 1.0E-06);

        ## ordered cut points for underlying continuous latent variable
        ordered ~ dconstraint(a0[1] < a0[2] & a0[2] < a0[3]);
        for(i in 1:3) {
            a0[i] ~ dnorm(0, 1.0E-6);
        }

        tau ~ dgamma(0.001, 0.001);
        sigma <- sqrt(1/tau);
        log.sigma <- log(sigma);
    })

    constants = new.env()
    inits = new.env()
    source(system.file('classic-bugs', 'vol1', 'inhaler','inhaler-inits.R', package = 'nimble'), local = inits)
    source(system.file('classic-bugs', 'vol1', 'inhaler','inhaler-data.R', package = 'nimble'), local = constants)
    pattern <- constants$pattern

    group <- c(rep(1, 59), rep(2, 63), rep(1, 35), rep(2, 13), rep(1, 3), rep(1, 2),
               rep(1, 11), rep(2, 40), rep(1, 27), rep(2, 15), rep(1, 2), rep(1, 1),
               rep(2, 7), rep(2, 2), rep(2, 1),
               rep(1, 1), rep(2, 2), rep(1,1), rep(2, 1))
    response <- matrix(c(rep(pattern[1, ], 59+63),
                         rep(pattern[2, ], 35+13),
                         rep(pattern[3, ], 3+0),
                         rep(pattern[4, ], 2+0),
                         rep(pattern[5, ], 11+40),
                         rep(pattern[6, ], 27+15),
                         rep(pattern[7, ], 2+0),
                         rep(pattern[8, ], 1+0),
                         rep(pattern[9, ], 0+7),
                         rep(pattern[10, ], 0+2),
                         rep(pattern[11, ], 0+1),
                         rep(pattern[13, ], 1+2),
                         rep(pattern[14, ], 1+0),
                         rep(pattern[15, ], 0+1)),
                       ncol = 2, byrow = TRUE)

    constants$group = group
    data = list(ordered = 1, response = response)

    Rmodel = nimbleModel(code, data = data, inits = as.list(inits), constants = as.list(constants))
    mcmcConf <- configureMCMC(Rmodel)
    Rmcmc <- buildMCMC(mcmcConf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    Cmcmc$run(5000)
    smp <- as.matrix(Cmcmc$mvSamples)
    expect_lt(max(smp[ , 'a0[1]'] - smp[ , 'a0[2]']), 0, label = "constraint on ordering of a[1] and a[2]")
    expect_lt(max(smp[ , 'a0[2]'] - smp[ , 'a0[3]']), 0, label = "constraint on ordering of a[2] and a[3]")

    test_mcmc(model = code, inits = as.list(inits), data = c(data, as.list(constants)), numItsC = 10000,
              resampleData = TRUE, basic = FALSE,
              results = list(mean = list(a0 = c(.72, 3.93, 5.31), beta = 1.05, pi = -0.24),
                             sd = list(a0 = c(.14, .36, .51), beta = .32, pi = .20)),
              resultsTolerance = list(mean = list(a0 = c(.05, .15, .15), beta = .1, pi = .04),
                                      sd = list(a0 = c(.02, .1, .1), beta = .03, pi = .02)),
              name = 'test of ordering constraint')
    ## no basic assessment because R MCMC takes forever, even for 5 iterations
})


# test that conjugate samplers not assigned when have dependent truncated node
# also that not assigned when target node is truncated (though we could code this
# so that we have the corrected truncated density as the conjugate sampler)

test_that("Test that MCMC with truncation avoids conjugate samplers", {
    code <- nimbleCode({
        y ~ T(dnorm(mu1, 1), 0, 3)
        mu1 ~ dnorm(0, 1)
        y2 ~ dnorm(mu2, 1)
        mu2 ~ T(dnorm(theta, 1), 0, 3)
        y3 ~ dnorm(mu3, 1)
        mu3 ~ dnorm(theta, 1)
    })
    
    m <- nimbleModel(code, data = list(y = 1), inits = list(mu2 = 1))
    conf <- configureMCMC(m)

expect_identical(conf$getSamplers('mu1')[[1]]$name, 'RW', info = "incorrectly assigning conjugate sampler for mu1")
expect_identical(conf$getSamplers('mu2')[[1]]$name, 'RW', info = "incorrectly assigning conjugate sampler for mu2")
expect_identical(conf$getSamplers('mu3')[[1]]$name, 'conjugate_dnorm_dnorm', info = "incorrectly not assigning conjugate sampler for mu3")
})

# test that truncation on discrete distribution correctly uses [L, U]
# and test use with alias too

test_that("Test that truncation with discrete distribution gives correct values", {
    code <- nimbleCode({
        y[1] ~ T(dbinom(p, n), 1, 2)
        y[2] ~ T(dbinom(p, n), 0.5, 2)
        y[3] ~ T(dbinom(p, n), a, 2)
        y[4] ~ T(dbinom(p, n), -Inf, 1)
        y[5] ~ T(dbinom(p, n), n-1, Inf)
        a ~ dunif(0.5, 0.8)
    })
    
    n <- 10; p <- 0.3
    expect_warning(m <- nimbleModel(code, constants = list(n = n, p = p),
                                    inits = list(a = 0.6, y = c(1,1,1,1,n-1))),
                   "Lower bound is less than or equal to distribution lower bound")
    cm <- compileNimble(m)
    
    tmp <- dbinom(1:2, n, p)
    dens_y123 <- tmp / sum(tmp)
    tmp <- dbinom(0:1, n, p)
    dens_y4 <- tmp / sum(tmp)
    tmp <- dbinom((n-1):n, n, p)
    dens_y5 <- tmp / sum(tmp)
    
    N <- 1000
    simVals <- array(0, c(N, 5, 2))
    set.seed(0)
    for(i in 1:N) {
        simulate(m)
        simVals[i, , 1] <- m$y
    }
    set.seed(0)
    for(i in 1:N) {
        simulate(cm)
        simVals[i, , 2] <- cm$y
    }
    
    simProbs <- matrix(0, nrow = 5, ncol = 2)
    for(j in 1:5) {
        if(j < 4) val <- 1 else if(j == 4) val <- 0 else val <- n-1
        simProbs[j, 1] <- mean(simVals[ , j, 1] == val)
        simProbs[j, 2] <- mean(simVals[ , j, 2] == val)
    }
    
    diffProbsR <- abs(simProbs - c(dens_y123[1], dens_y123[1], dens_y123[1], dens_y4[1], dens_y5[1]))
    
    expect_lt(max(diffProbsR), 0.015) # "simulation does not give approximate density values"
    
    m$y[1] <- 1
    expect_equal(calculate(m, 'y[1]'), (log(dens_y123[1])), info = "incorrect R model density for y[1]=1")
    m$y[1] <- 2
    expect_equal(calculate(m, 'y[1]'), (log(dens_y123[2])), info = "incorrect R model density for y[1]=2")
    m$y[2] <- 1
    expect_equal(calculate(m, 'y[2]'), (log(dens_y123[1])), info = "incorrect R model density for y[2]=1")
    m$y[2] <- 2
    expect_equal(calculate(m, 'y[2]'), (log(dens_y123[2])), info = "incorrect R model density for y[2]=2")
    m$y[3] <- 1
    expect_equal(calculate(m, 'y[3]'), (log(dens_y123[1])), info = "incorrect R model density for y[3]=1")
    m$y[3] <- 2
    expect_equal(calculate(m, 'y[3]'), (log(dens_y123[2])), info = "incorrect R model density for y[3]=2")
    m$y[4] <- 0
    expect_equal(calculate(m, 'y[4]'), (log(dens_y4[1])), info = "incorrect R model density for y[4]=0")
    m$y[4] <- 1
    expect_equal(calculate(m, 'y[4]'), (log(dens_y4[2])), info = "incorrect R model density for y[4]=1")
    m$y[5] <- n-1
    expect_equal(calculate(m, 'y[5]'), (log(dens_y5[1])), info = "incorrect R model density for y[5]=n-1")
    m$y[5] <- n
    expect_equal(calculate(m, 'y[5]') , (log(dens_y5[2])), info = "incorrect R model density for y[5]=n")

    cm$y[1] <- 1
    expect_equal(calculate(cm, 'y[1]'), (log(dens_y123[1])), info = "incorrect C model density for y[1]=1")
    cm$y[1] <- 2
    expect_equal(calculate(cm, 'y[1]'), (log(dens_y123[2])), info = "incorrect C model density for y[1]=2")
    cm$y[2] <- 1
    expect_equal(calculate(cm, 'y[2]'), (log(dens_y123[1])), info = "incorrect C model density for y[2]=1")
    cm$y[2] <- 2
    expect_equal(calculate(cm, 'y[2]'), (log(dens_y123[2])), info = "incorrect C model density for y[2]=2")
    cm$y[3] <- 1
    expect_equal(calculate(cm, 'y[3]'), (log(dens_y123[1])), info = "incorrect C model density for y[3]=1")
    cm$y[3] <- 2
    expect_equal(calculate(cm, 'y[3]'), (log(dens_y123[2])), info = "incorrect C model density for y[3]=2")
    cm$y[4] <- 0
    expect_equal(calculate(cm, 'y[4]'), (log(dens_y4[1])), info = "incorrect C model density for y[4]=0")
    cm$y[4] <- 1
    expect_equal(calculate(cm, 'y[4]'), (log(dens_y4[2])), info = "incorrect C model density for y[4]=1")
    cm$y[5] <- n-1
    expect_equal(calculate(cm, 'y[5]'), (log(dens_y5[1])), info = "incorrect C model density for y[5]=n-1")
    cm$y[5] <- n
    expect_equal(calculate(cm, 'y[5]'), (log(dens_y5[2])), info = "incorrect C model density for y[5]=n")
})

sink(NULL)

if(!generatingGoldFile) {
    test_that("Log file matches gold file", {
        ## Not clear why "[1] TRUE" is showing up at start of trialResults (as of 2021-01)
        if(trialResults[1] == "[1] TRUE")
            trialResults <- trialResults[-1]
            trialResults <- readLines(tempFileName)
            correctResults <- readLines(system.file(file.path('tests', 'testthat', goldFileName), package = 'nimble'))
            compareFilesByLine(trialResults, correctResults)
    })
}

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)
options(width = oldWidth)
options(max.print = oldMaxPrint)
