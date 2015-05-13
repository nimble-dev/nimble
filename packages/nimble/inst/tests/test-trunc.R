source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of truncation, censoring, and constraints")

# should build testing of kidney,litters,mice,lsat into test-mcmc and test-models

models <- c('kidney', 'litters', 'mice')

### test basic model building
sapply(models, testBUGSmodel, useInits = TRUE)

# at the moment we haven't updated conjugacy to deal with truncation
# so litters is reporting conj post density as wrong since
# otherwise p[i,j] is a standard conjugate update

# lsat has non-explicit indexing
system(paste("cat", system.file('classic-bugs','vol1','lsat','lsat.bug', package = 'nimble'), ">>", file.path(tempdir(), "lsat.bug")))
system(paste("sed -i -e 's/mean(alpha\\[\\])/mean(alpha\\[1:T\\])/g'", file.path(tempdir(), "lsat.bug"))) 
system(paste("sed -i -e 's/p.item\\[i,\\]/p.item\\[i,1:T\\]/g'", file.path(tempdir(), "lsat.bug"))) 
testBUGSmodel('lsat', dir = "", model = file.path(tempdir(), "lsat.bug"), data = system.file('classic-bugs','vol1','lsat','lsat-data.R', package = 'nimble'),  inits = system.file('classic-bugs','vol1','lsat','lsat-init.R', package = 'nimble'),  useInits = TRUE)


### test MCMC
sapply(models, test_mcmc, numItsC = 1000)

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
              data = data, inits = inits)

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
mcmcspec <- configureMCMC(Rmodel)
Rmcmc <- buildMCMC(mcmcspec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
Cmcmc$run(5000)
smp <- as.matrix(Cmcmc$mvSamples)
try(test_that("Test that MCMC respects truncation bounds: ", expect_that(min(smp), is_more_than(0.2), info = "minimum in MCMC greater than lower bound")))
try(test_that("Test that MCMC respects truncation bounds: ", expect_that(min(smp), is_less_than(2), info = "maximum in MCMC less than upper bound")))

test_mcmc(model = code, data = c(data, constants), inits = inits,
          results = list(mean = list(mu = 1.5), sd = list(mu = .27 )),
          resultsTolerance = list(mean = list(mu = 0.025), sd = list(mu = .025)), name = 'test of MCMC with truncation')

### censoring (dinterval)
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

constants <- list(n = n)
data <- list(y = y, is.cens = is.cens)
inits <- list(mu = mu, bnd = bnd, y = yinits)

Rmodel = nimbleModel(code, data = data, inits = inits,
    constants = constants)
mcmcspec <- configureMCMC(Rmodel)
mcmcspec$addMonitors('y[1]')
Rmcmc <- buildMCMC(mcmcspec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
Cmcmc$run(5000)
smp <- as.matrix(Cmcmc$mvSamples)
try(test_that("Test that MCMC respects right-censoring bounds: ", expect_that(min(smp[ , 'bnd'] - max(y[!is.cens])), is_more_than(0), info = "upper bound in MCMC exceeds observations")))
try(test_that("Test that MCMC respects censoring bounds: ",
              expect_that(min(smp[ , 'y[1]'] - smp[ , 'bnd']), is_more_than(0), info = "minimum inferred censored value in MCMC greater than bound")))


test_mcmc(model = code, data = c(data, constants), inits = inits, numItsC = 10000,
          results = list(mean = list(mu = 44.5, 'y[1]' = 56.6), sd = list(mu = 2.3, 'y[1]' = 4.95)),
          resultsTolerance = list(mean = list(mu = 0.1, 'y[1]' = 1.5), sd = list(mu = .3, 'y[1]' = .7)), name = 'test of right censoring')

# left censored
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

constants <- list(n = n)
data <- list(y = y, not.cens = not.cens)
inits <- list(mu = mu, bnd = bnd, y = yinits)

Rmodel = nimbleModel(code, data = data, inits = inits,
    constants = constants)
mcmcspec <- configureMCMC(Rmodel)
mcmcspec$addMonitors('y[1]')
Rmcmc <- buildMCMC(mcmcspec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
Cmcmc$run(5000)
smp <- as.matrix(Cmcmc$mvSamples)
try(test_that("Test that MCMC respects left-censoring bounds: ", expect_that(min(min(y[not.cens]) - smp[ , 'bnd']), is_more_than(0), info = "lower bound in MCMC is less than observations")))
try(test_that("Test that MCMC respects censoring bounds: ",
              expect_that(max(smp[ , 'bnd'] - smp[ , 'y[1]']), is_more_than(0), info = "maximum inferred censored value in MCMC less than bound")))

test_mcmc(model = code, data = c(data, constants), inits = inits,
          results = list(mean = list(mu = 43.4, 'y[1]' = 33), sd = list(mu = 2.4, 'y[1]' = 5.3)),
          resultsTolerance = list(mean = list(mu = .2, 'y[1]' = .4), sd = list(mu = .2, 'y[1]' = .1)), name = 'test of left censoring')


# interval censored
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
inits <- list(mu = mu, y = yinits)

Rmodel = nimbleModel(code, data = data, inits = inits, constants = constants)
mcmcspec <- configureMCMC(Rmodel)
mcmcspec$addMonitors(c('y[2]', 'y[7]', 'y[12]'))
Rmcmc <- buildMCMC(mcmcspec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
Cmcmc$run(5000)
smp <- as.matrix(Cmcmc$mvSamples)
try(test_that("Test that MCMC respects first interval: ", expect_that(max(smp[ , 'y[1]']), is_less_than(bnd[1]), info = "imputed observation in first interval is too large")))
try(test_that("Test that MCMC respects second interval: ", expect_that(min(smp[ , 'y[7]']), is_more_than(bnd[1]), info = "imputed observation in second interval is too small")))
try(test_that("Test that MCMC respects second interval: ", expect_that(max(smp[ , 'y[7]']), is_less_than(bnd[2]), info = "imputed observation in second interval is too large")))
try(test_that("Test that MCMC respects third interval: ", expect_that(min(smp[ , 'y[12]']), is_more_than(bnd[2]), info = "imputed observation in third interval is too small")))

test_mcmc(model = code, data = c(data, constants), inits = inits,
          results = list(mean = list(mu = 44.77, 'y[12]' = 44.85), sd = list(mu = 2.8, 'y[12]' = 2.9)),
          resultsTolerance = list(mean = list(mu = 0.5, 'y[12]' = .15), sd = list(mu = .05, 'y[12]' = .1)), name = 'test of interval censoring')


# test of dconstraint
set.seed(0)
code <- nimbleCode ({
    for(i in 1:n) {
        y1[i] ~ dnorm(mu1, 1)
        y2[i] ~ dnorm(mu2, 1)
    }
    mu1 ~ dunif(-10, 10)
    mu2 ~ dunif(-10, 10)
    ind <- mu1 + mu2 > 0 & mu1 > 1
#    mu.c ~ dconstraint(mu1 + mu2 > 0 & mu1 > 1) # ind)
    mu.c ~ dconstraint(ind)
})
mu1 <- 1; mu2 <- -0.5
n <- 10
y1 <- rnorm(n, mu1, 1)
y2 <- rnorm(n, mu2, 1)
data = list(y1 = y1, y2 = y2, mu.c =1)
inits = list(mu1 = 1.5, mu2 = 0.5)
constants = list(n = n)

Rmodel = nimbleModel(code, data = data, inits = inits, constants = constants)
mcmcspec <- configureMCMC(Rmodel)
Rmcmc <- buildMCMC(mcmcspec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
Cmcmc$run(5000)
smp <- as.matrix(Cmcmc$mvSamples)
try(test_that("Test that MCMC respects constraints: ", expect_that(min(smp[ , 'mu1'] + smp[ , 'mu2']), is_more_than(0), info = "constraint on sum of mu1 and mu2 is violated")))
try(test_that("Test that MCMC respects constraints: ", expect_that(min(smp[ , 'mu1']), is_more_than(0), info = "constraint on mu1 is violated")))

test_mcmc(model = code, data = c(data, constants), inits = inits,
          results = list(mean = list(mu1 = 1.45, mu2 = -.82), sd = list(mu1 = .26, mu2 = .30)),
          resultsTolerance = list(mean = list(mu1 = 0.05, mu2 = .02), sd = list(mu1 = .03, mu2 = .03)), name = 'test of dconstraint')


# test dconstraint for ordering using rewrite of inhaler in our syntax

# pasting BUGS code here because inhaler has a data black that has non-standard
# use of ":" in for loop indexing, and so we can substitute our use of
# dconstraint in place of JAGS sort()

code <- nimbleCode({
    for (i in 1:N) {
        for (t in 1:T) {
            for (j in 1:Ncut) {
#  
# Cumulative probability of worse response than j
#
                logit(Q[i,t,j]) <- -(a0[j] + mu[group[i],t] + b[i]); 
            }
#
# Probability of response = j
#
            p[i,t,1] <- 1 - Q[i,t,1];
            for (j in 2:Ncut) { p[i,t,j] <- Q[i,t,j-1] - Q[i,t,j] }
            p[i,t,(Ncut+1)] <- Q[i,t,Ncut];
            
            response[i,t] ~ dcat(p[i,t,1:4]);
        }
#
# Subject (random) effects
#
        b[i] ~ dnorm(0.0, tau);
    }
    
#
# Fixed effects
#
    for (g in 1:G) {
        for(t in 1:T) { 
                                        # logistic mean for group i in period t
            mu[g,t] <- beta*treat[g,t]/2 + pi*period[g,t]/2 + kappa*carry[g,t]; 
        }
    }                                                             
    beta ~ dnorm(0, 1.0E-06);
    pi ~ dnorm(0, 1.0E-06);
    kappa ~ dnorm(0, 1.0E-06);
    
# ordered cut points for underlying continuous latent variable  
    ind <- a0[1] < a0[2] & a0[2] < a0[3];
    ordered ~ dconstraint(ind);
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
mcmcspec <- configureMCMC(Rmodel)
Rmcmc <- buildMCMC(mcmcspec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
Cmcmc$run(5000)
smp <- as.matrix(Cmcmc$mvSamples)
try(test_that("Test that MCMC respects constraints: ", expect_that(max(smp[ , 'a0[1]'] - smp[ , 'a0[2]']), is_less_than(0), info = "constraint on ordering of a[1] and a[2] is violated")))
try(test_that("Test that MCMC respects constraints: ", expect_that(max(smp[ , 'a0[2]'] - smp[ , 'a0[3]']), is_less_than(0), info = "constraint on ordering of a[2] and a[3] is violated")))

test_mcmc(model = code, inits = as.list(inits), data = c(data, as.list(constants)), numItsC = 10000,
          resampleData = TRUE, basic = FALSE,
        results = list(mean = list(a0 = c(.72, 3.93, 5.31), beta = 1.05, pi = -0.24),
            sd = list(a0 = c(.14, .36, .51), beta = .32, pi = .20)),
          resultsTolerance = list(mean = list(a0 = c(.05, .15, .15), beta = .1, pi = .04),
              sd = list(a0 = c(.02, .1, .1), beta = .03, pi = .02)),
          name = 'test of ordering contraint')
# no basic assessment because R MCMC takes forever, even for 5 iterations


# test that conjugate samplers not assigned when have dependent truncated node
# also that not assigned when target node is truncated (though we could code this
# so that we have the corrected truncated density as the conjugate sampler)

code <- nimbleCode({
   y ~ T(dnorm(mu1, 1), 0, 3)
   mu1 ~ dnorm(0, 1)
   y2 ~ dnorm(mu2, 1)
   mu2 ~ T(dnorm(theta, 1), 0, 3)
   y3 ~ dnorm(mu3, 1)
   mu3 ~ dnorm(theta, 1)
})

m <- nimbleModel(code)
spec <- configureMCMC(m)

try(test_that("Test that MCMC with truncation avoids conjugate samplers: ", expect_that(spec$samplerSpecs[[1]]$name, is_identical_to('RW'), info = "incorrectly assigning conjugate sampler for mu1")))
try(test_that("Test that MCMC with truncation avoids conjugate samplers: ", expect_that(spec$samplerSpecs[[3]]$name, is_identical_to('RW'), info = "incorrectly assigning conjugate sampler for mu2")))
try(test_that("Test that MCMC with truncation avoids conjugate samplers: ", expect_that(spec$samplerSpecs[[4]]$name, is_identical_to('conjugate_dnorm'), info = "incorrectly not assigning conjugate sampler for mu3")))
   
