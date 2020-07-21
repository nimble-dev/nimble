source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of default MCMC, part2 ")

### From here down is copied from the 2nd half of test-mcmc.R
## search for "windows" to see how it is split
## On Windows we cannot run all the way through test-mcmc, so instead it bails
## at a certain point.  On can then run test-mcmc2 to pick up in this file.

### Daniel's world's simplest MCMC demo

code <- nimbleCode({
    x ~ dnorm(0, 2)
    y ~ dnorm(x+1, 3)
    z ~ dnorm(y+2, 4)
})
data = list(y = 3)

test_mcmc(model = code, name = 'very simple example', data = data, resampleData = FALSE, results = list(
                                       mean = list(x = 6/5, z = 5),
                                       sd = list(x = 1/sqrt(5), z = 1/2)),
          resultsTolerance = list(mean = list(x = .1, z = .1),
            sd = list(x = .05, z = .05)))

### basic block sampler example

code <- nimbleCode({
    for(i in 1:3) {
        x[i] ~ dnorm(0, 1)
        y[i] ~ dnorm(x[i], 2)
    }
})
data = list(y = -1:1)

test_mcmc(model = code, name = 'basic no-block sampler', data = data, resampleData = FALSE, results = list(
                                       mean = list(x = c(-2/3,0,2/3)),
                                       var = list(x = rep(1/3,3))),
          resultsTolerance = list(mean = list(x = rep(.1,3)),
            var = list(x = rep(.05,3))))

test_mcmc(model = code, name = 'basic block sampler on scalars', data = data, resampleData = FALSE, results = list(
                                       mean = list(x = c(-2/3,0,2/3)),
                                       var = list(x = rep(1/3,3))),
          resultsTolerance = list(mean = list(x = rep(.1,3)),
            var = list(x = rep(.05,3))),
          samplers = list(
            list(type = 'RW_block', target = 'x[1]'),
            list(type = 'RW_block', target = 'x[2]'),
            list(type = 'RW_block', target = 'x[3]')
            ), removeAllDefaultSamplers = TRUE, numItsC = 10000)

test_mcmc(model = code, name = 'basic block sampler on vector', data = data, resampleData = FALSE, results = list(
                                       mean = list(x = c(-2/3,0,2/3)),
                                       var = list(x = rep(1/3,3))),
          resultsTolerance = list(mean = list(x = rep(.1,3)),
            var = list(x = rep(.05,3))),
          samplers = list(
            list(type = 'RW_block', target = 'x', control = list(adaptInterval = 500))
            ), numItsC = 10000)



### slice sampler example

code <- BUGScode({
    z ~ dnorm(0, 1)
    normal5_10 ~ dnorm(5, sd = 10)
    beta1_1 ~ dbeta(1, 1)
    beta3_5 ~ dbeta(3, 5)
    binom10_p5 ~ dbin(size=10, prob=0.5)
    binom20_p3 ~ dbin(size=20, prob=0.3)
})

test_mcmc(model = code, name = "slice sampler example", resampleData = FALSE, results = list(
                                       mean = list(z = 0, "beta1_1" = 0.5, "beta3_5" = 3/(3+5),
                                         "binom10_p5" = 10*.5, "binom20_p3" = 20*.3),
                                         sd = list(z = 1, "beta1_1" = sqrt(1/12),
                                           "beta3_5" = sqrt(3*5/((3+5)^2*(3+5+1))),
                                          "binom10_p5" = sqrt(10*.5*.5),
                                           "binom20_p3" = sqrt(20*.3*.7))),
          resultsTolerance = list(
                                       mean = list(z = 0.1, "beta1_1" = 0.5, "beta3_5" = .2,
                                         "binom10_p5" = .25, "binom20_p3" = .25),
                                         sd = list(z = .1, "beta1_1" = .05, "beta3_5" = .03,
                                          "binom10_p5" = .2, "binom20_p3" = .25)),
          samplers = list(list(type = 'slice', target = 'z', control = list(adaptInterval = 10)),
            list(type = 'slice', target = 'normal5_10', control = list(adaptInterval = 10)),
           list(type = 'slice', target = 'beta1_1', control = list(adaptInterval = 10)),
            list(type = 'slice', target = 'beta3_5', control = list(adaptInterval = 10)),
            list(type = 'slice', target = 'binom10_p5', control = list(adaptInterval = 10)),
            list(type = 'slice', target = 'binom20_p3', control = list(adaptInterval = 10))))




### AF_slice sampler example, default control options.


code <- nimbleCode({
  mu[1] <- 10
  mu[2] <- 20
  mu[3] <- 30
  x[1:3] ~ dmnorm(mu[1:3], prec = Q[1:3,1:3])
})

Q = matrix(c(1.0,0.2,-1.0,0.2,4.04,1.6,-1.0,1.6,10.81), nrow=3)
data = list(Q = Q)
inits = list(x = c(10, 20, 30))

test_mcmc(model = code, name = 'block sampler on multivariate node', data = data, seed = 0, numItsC = 10000,
          results = list(mean = list(x = c(10,20,30)),
                         var = list(x = diag(solve(Q)))),
          resultsTolerance = list(mean = list(x = rep(1,3)),
                                  var = list(x = c(.1, .03, .01))),
          samplers = list(list(type = 'AF_slice', target = c('x[1:3]'))))


### elliptical slice sampler 'ess'

set.seed(0)
ESScode <- quote({
    x[1:d] ~ dmnorm(mu_x[1:d], prec = prec_x[1:d, 1:d])
    y[1:d] ~ dmnorm(x[1:d], prec = prec_y[1:d, 1:d])
})
d <- 3
mu_x <- rnorm(d)
temp <- array(rnorm(d^2), c(d,d))
prec_x <- solve(temp %*% t(temp))
temp <- array(rnorm(d^2), c(d,d))
prec_y <- solve(temp %*% t(temp))
y <- rnorm(d)
ESSconstants <- list(d = d, mu_x = mu_x, prec_x = prec_x, prec_y = prec_y)
ESSdata <- list(y = y)
ESSinits <- list(x = rep(0, d))

test_mcmc(model = ESScode, data = c(ESSconstants, ESSdata), inits = ESSinits,
          name = 'exact values of elliptical slice sampler',
          seed = 0,
          exactSample = list(`x[1]` = c(-0.492880566939352, -0.214539223107114, 1.79345037297218, 1.17324496091208, 2.14095077672555, 1.60417482445964, 1.94196916651627, 2.66737323347255, 2.66744178776022, 0.253966883192744), `x[2]` = c(-0.161210109217102, -0.0726534676226932, 0.338308532423757, -0.823652445515156, -0.344130712698579, -0.132642244861469, -0.0253168895009594, 0.0701624130921676, 0.0796842215444978, -0.66369112443311), `x[3]` = c(0.278627475932455, 0.0661336950029345, 0.407055002920732, 1.98761228946318, 1.0839897275519, 1.00262648370199, 0.459841485268785, 2.59229443025387, 1.83769567435409, 1.92954706515119)),
          samplers = list(list(type = 'ess', target = 'x')))

test_mcmc(model = ESScode, data = c(ESSconstants, ESSdata), inits = ESSinits,
          name = 'results to tolerance of elliptical slice sampler',
          results = list(mean = list(x = c(1.0216463, -0.4007247, 1.1416904))),
          resultsTolerance = list(mean = list(x = c(0.01, 0.01, 0.01))),
          numItsC = 100000,
          samplers = list(list(type = 'ess', target = 'x')))



### demo2 of check conjugacy

code <- BUGScode({
    x ~ dbeta(3, 13)
    y[1] ~ dbin(x, 10)
    y[2] ~ dbin(x, 20)
})
data = list(y = c(3,4))

test_mcmc(model = code, name = 'check of beta-binom conjugacy', data = data, exactSample = list(x = c(0.195510839527966, 0.332847482503424,0.247768152764931, 0.121748195439553, 0.157842271774841, 0.197566496350904, 0.216991517500577, 0.276609942874852, 0.165733872345582, 0.144695512780252)), seed = 0)

### checkConjugacy_demo3_run.R - various conjugacies

code <- BUGScode({
    x ~ dgamma(1, 1)       # should satisfy 'gamma' conjugacy class
    a  ~ dnorm(0, x)     # should satisfy 'norm' conjugacy class
    a2 ~ dnorm(0, tau = 3*x+0)
    b  ~ dpois(0+5*x)
    b2 ~ dpois(1*x*1)
    c ~ dgamma(1, 7*x*5)
    for(i in 2:3) {
        jTau[i] <- 1
        jNorm[i] ~ dnorm(c * (a+3) - i, var = jTau[i])
        kTauSd[i] <- 2
        kLogNorm[i] ~ dlnorm(0 - a - 6*i, kTauSd[i])
    }
})

    sampleVals = list(x = c(3.950556165467749, 1.556947815895538, 1.371834934033851, 2.036442813764752, 2.247416118159410, 2.537131924778210, 2.382184991769738, 2.653737836857812, 2.934255734970981, 3.007873553270551),
                      c = c(0.010341199485849559, 0.010341199485849559, 0.003846483017887228, 0.003846483017887228, 0.003846483017887228, 0.006269117826484087, 0.009183580181658716, 0.009183580181658716, 0.006361841408434201, 0.006361841408434201))

test_mcmc(model = code, name = 'check various conjugacies', exactSample = sampleVals, seed = 0, mcmcControl = list(scale=0.01))

### Dirichlet-multinomial conjugacy

# as of v0.4, exact numerical results here have changed because
# ddirch now sometimes returns NaN rather than -Inf (when an
# alpha is proposed to be negative) -- this changes the RNG
# sequence because NaN values result in no runif() call in decide()

# single multinomial
set.seed(0)
n <- 100
alpha <- c(10, 30, 15, 60, 1)
K <- length(alpha)
p <- c(.12, .24, .09, .54, .01)
y <- rmulti(1, n, p)

code <- function() {
    y[1:K] ~ dmulti(p[1:K], n);
    p[1:K] ~ ddirch(alpha[1:K]);
    for(i in 1:K) {
        alpha[i] ~ dgamma(.001, .001);
    }
}

inits <- list(p = rep(1/K, K), alpha = rep(K, K))
data <- list(n = n, K = K, y = y)

test_mcmc(model = code, name = 'Dirichlet-multinomial example', data= data, seed = 0, numItsC = 10000,
          inits = inits,
          results = list(mean = list(p = p)),
          resultsTolerance = list(mean = list(p = rep(.06, K))))

# bad mixing for alphas; probably explains why posterior estimates for alphas changed so much as of v 0.4

# with replication

set.seed(0)
n <- 100
m <- 20
alpha <- c(10, 30, 15, 60, 1)
K <- length(alpha)
y <- p <- matrix(0, m, K)
for(i in 1:m) {
    p[i, ] <- rdirch(1, alpha)
    y[i, ] <- rmulti(1, n, p[i, ])
}

code <- function() {
    for(i in 1:m) {
        y[i, 1:K] ~ dmulti(p[i, 1:K], n);
        p[i, 1:K] ~ ddirch(alpha[1:K]);
    }
    for(i in 1:K) {
        alpha[i] ~ dgamma(.001, .001);
    }
}

inits <- list(p = matrix(1/K, m, K), alpha = rep(1/K, K))
data <- list(n = n, K = K, m = m, y = y)

test_mcmc(model = code, name = 'Dirichlet-multinomial with replication', data= data, seed = 0, numItsC = 1000,
          inits = inits, numItsC_results = 100000,
          results = list(mean = list(p = p, alpha = alpha)),
          resultsTolerance = list(mean = list(p = matrix(.05, m, K),
                                      alpha = c(5,10,10,20,.5))))

# note alphas mix poorly (and are highly correlated),
# presumably because of cross-level dependence between
# p's and alphas.  cross-level sampler would probably work well here,
# or, of course, integrating over the p's

### block sampler on MVN node

code <- nimbleCode({
    mu[1] <- 10
    mu[2] <- 20
    mu[3] <- 30
    x[1:3] ~ dmnorm(mu[1:3], prec = Q[1:3,1:3])
})

Q = matrix(c(1.0,0.2,-1.0,0.2,4.04,1.6,-1.0,1.6,10.81), nrow=3)
data = list(Q = Q)
inits = list(x = c(10, 20, 30))

test_mcmc(model = code, name = 'block sampler on multivariate node', data = data, seed = 0, numItsC = 10000,
          results = list(mean = list(x = c(10,20,30)),
          var = list(x = diag(solve(Q)))),
          resultsTolerance = list(mean = list(x = rep(1,3)),
            var = list(x = c(.1, .03, .01))),
          samplers = list(
            list(type = 'RW_block', target = 'x[1:3]')))
# caution: setting targetNodes='x' works but the initial end sampler is not removed because x[1:3] in targetNode in default sampler != 'x' in targetNodes passed in
if(FALSE) {
    Rmodel <- nimbleModel(code, constants = list(Q=Q))
    mcmcspec <- MCMCspec(Rmodel, nodes = NULL)
    mcmcspec$addSampler(type = 'RW_block', target = 'x', control = list(adaptInterval=500))
    mcmcspec$getMonitors()
    Rmcmc <- buildMCMC(mcmcspec)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    Cmcmc(200000)    ## this runs nearly instantaneously on my computer -DT
    samples <- as.matrix(nfVar(Cmcmc, 'mvSamples'))
    samples <- samples[50001:200000,]
    dim(samples)
    apply(samples, 2, mean)
    solve(Q)
    cov(samples)
    propCov <- nfVar(Cmcmc, 'samplerFunctions')[[1]]$propCov
    scale <- nfVar(Cmcmc, 'samplerFunctions')[[1]]$scale
    propCov * scale^2

    nfVar(Cmcmc, 'samplerFunctions')[[1]]$scaleHistory
    nfVar(Cmcmc, 'samplerFunctions')[[1]]$acceptanceRateHistory
    nfVar(Cmcmc, 'samplerFunctions')[[1]]$scale
    nfVar(Cmcmc, 'samplerFunctions')[[1]]$propCov
    ## why is the proposal cov w/ .99 cross-corrs?
    ## also MCMC in C takes a surprisingly long time - this might be threaded lin alg behaving badly on small matrices
}

### DT's model
mu <- c(1,2,3)
corr <- matrix(c(1,.8,0.3,.8,1,0,0.3,0,1), nrow=3)
varr <- c(1,2,3)
Sig <- diag(sqrt(varr))
Q <- Sig %*% corr %*% Sig
P <- solve(Q)

code <- nimbleCode({
#    x[1:3] ~ dmnorm(mu[1:3], cov = Q[1:3,1:3])
    x[1:3] ~ dmnorm(mu[1:3], prec = P[1:3,1:3])
})
data = list(P = P, mu = mu)

test_mcmc(model = code, name = 'second block sampler on multivariate node', data = data, seed = 0, numItsC = 100000,
          results = list(mean = list(x = mu),
          var = list(x = varr)),
          resultsTolerance = list(mean = list(x = rep(.1,3)),
            var = list(x = c(.1,.1,.1))),
          samplers = list(
            list(type = 'RW_block', target = 'x[1:3]')))



### MVN conjugate update

set.seed(0)
mu0 = 1:3
Q0 = matrix(c(1, .2, .8, .2, 2, 1, .8, 1, 2), nrow = 3)
Q = solve(matrix(c(3, 1.7, .9, 1.7, 2, .6, .9, .6, 1), nrow = 3))
a = c(-2, .5, 1)
B = matrix(rnorm(9), 3)

##### not currently working - see Perry's email of ~ 10/6/14
## code <- nimbleCode({
##   mu[1:3] ~ dmnorm(mu0[1:3], Q0[1:3, 1:3])
##   y[1:3] ~ dmnorm(asCol(a[1:3]) + B[1:3, 1:3] %*% asCol(mu[1:3]), Q[1:3, 1:3])
## })

code <- nimbleCode({
  mu[1:3] ~ dmnorm(mu0[1:3], Q0[1:3, 1:3])
  y_mean[1:3] <- asCol(a[1:3]) + B[1:3, 1:3] %*% asCol(mu[1:3])
  y[1:3] ~ dmnorm(y_mean[1:3], Q[1:3, 1:3])
})

## Simplest version of model w/o 'a' and 'B'
## a = rep(0,3)
## B = diag(rep(1,3))
## code <- nimbleCode({
##   mu[1:3] ~ dmnorm(mu0[1:3], Q0[1:3, 1:3])
##   y[1:3] ~ dmnorm(mu[1:3], Q[1:3, 1:3])
## })


mu <- mu0 + chol(solve(Q0)) %*% rnorm(3)
# make sure y is a vec not a 1-col matrix or get a dimensionality error
y <- c(a + B%*%mu + chol(solve(Q)) %*% rnorm(3))
data = list(mu0 = mu0, Q0 = Q0, Q = Q, a = a, B = B, y = y)

muQtrue = t(B) %*% Q%*%B + Q0
muMeanTrue = c(solve(muQtrue, crossprod(B, Q%*%(y-a)) + Q0%*%mu0))

test_mcmc(model = code, name = 'two-level multivariate normal', data = data, seed = 0, numItsC = 10000,
          results = list(mean = list(mu = muMeanTrue),
                           cov = list(mu = solve(muQtrue))),
          resultsTolerance = list(mean = list(mu = rep(.02,3)),
            cov = list(mu = matrix(.01, 3, 3))))


### scalar RW updates in place of conjugate mv update

test_mcmc(model = code, name = 'two-level multivariate normal with scalar updaters', data = data, seed = 0, numItsC = 100000,
          results = list(mean = list(mu = muMeanTrue),
                           cov = list(mu = solve(muQtrue))),
          resultsTolerance = list(mean = list(mu = rep(.03,3)),
            cov = list(mu = matrix(.03, 3, 3))),
          samplers = list(list(type = 'RW', target = 'mu[1]'),
            list(type = 'RW', target = 'mu[2]'),
            list(type = 'RW', target = 'mu[3]')),
          removeAllDefaultSamplers = TRUE)


## section 3 of tests
message("Continue with test-mcmc3")

