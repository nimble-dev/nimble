source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

context("Testing of default MCMC, part 3 ")

### From here down is copied from the 2nd half of test-mcmc2.R
## On Windows we cannot run all the way through test-mcmc, so instead it bails
## at a certain point.  On can then run test-mcmc2 and test-mcmc3 to pick up in this file.

## another example of MVN conjugate sampler, for test-mcmc.R
## using both cov and prec parametrizaions of MVN,
## and various linear links

set.seed(0)
prior_mean <- rep(0,5)
tmp <- array(rnorm(25), c(5,5))
tmp <- tmp + t(tmp) + 5*diag(5)
prior_cov <- tmp
a <- array(rnorm(20), c(4,5))
B <- array(NA, c(4,5,5))
for(i in c(2,4))   B[i,,] <- array(rnorm(25), c(5,5))
B[1,,] <- diag(5)
B[3,,] <- diag(5)
M_y <- array(NA, c(4,5,5))
for(i in 1:4) {
    tmp <- array(rnorm(25,i), c(5,5))
    tmp <- tmp + t(tmp) + 5*i*diag(5)
    M_y[i,,] <- tmp
}
x <- rep(0, 5)
y <- array(rnorm(20), c(4,5))

code <- nimbleCode({
    x[1:5] ~ dmnorm(mean = prior_mean[1:5], cov = prior_cov[1:5,1:5])
    for(i in 1:4)
        mu_y[i,1:5] <- asCol(a[i,1:5]) + B[i,1:5,1:5] %*% asCol(x[1:5])
    y[1,1:5] ~ dmnorm(mu_y[1,1:5], prec = M_y[1,1:5,1:5])
    y[2,1:5] ~ dmnorm(mu_y[2,1:5], cov  = M_y[2,1:5,1:5])
    y[3,1:5] ~ dmnorm(mu_y[3,1:5], prec = M_y[3,1:5,1:5])
    y[4,1:5] ~ dmnorm(mu_y[4,1:5], cov  = M_y[4,1:5,1:5])
})
constants <- list(prior_mean=prior_mean, prior_cov=prior_cov, a=a, B=B, M_y=M_y)
data <- list(y=y)
inits <- list(x=x)
Rmodel <- nimbleModel(code, constants, data, inits)
spec <- configureMCMC(Rmodel)
##spec$getSamplers()
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)
Rsamples <- as.matrix(Rmcmc$mvSamples)
set.seed(0)
Cmcmc$run(10)
Csamples <- as.matrix(Cmcmc$mvSamples)

test_that(
    'expected R sample',
    expect_equal(round(as.numeric(Rsamples), 8),
                 ##cat('c(', paste0(as.numeric(round(Rsamples,8)), collapse=', '), ')\n')
                 c(0.97473128, 0.50438666, 1.1251132, 0.83830666, 0.74077066, 0.92935482, 0.83758372, 0.98708273, 1.24199937, 0.67348127, -0.54387714, -0.60713969, -0.51392796, -0.3176801, -0.34416529, -0.08530564, -0.47160157, -0.21996584, -0.20504917, -0.77287122, 0.78462584, 0.46103509, 0.43862813, 0.49343096, 0.61020864, 0.55088287, 0.53887202, 0.49863894, 0.62691318, 0.80142839, 0.34941152, 0.06623608, 0.05624477, 0.21369178, 0.26585415, -0.1439989, -0.03133488, 0.3544062, -0.03518959, 0.27415746, 0.40977, 0.8351078, 0.25719293, 0.05663917, 0.30894028, 0.33113315, 0.47647909, 0.26143962, 0.07180759, 0.27255767)
                 ))

dif <- as.numeric(Rsamples - Csamples)
test_that('R and C equiv', expect_lt(max(abs(dif)), 1E-15))

y_prec <- array(NA, c(4,5,5))
y_prec[1,,] <-       M_y[1,,]
y_prec[2,,] <- solve(M_y[2,,])
y_prec[3,,] <-       M_y[3,,]
y_prec[4,,] <- solve(M_y[4,,])
contribution_mean <- array(NA, c(4,5))
for(i in 1:4)   contribution_mean[i,] <- t(B[i,,]) %*% y_prec[i,,] %*% (y[i,] - a[i,])
contribution_prec <- array(NA, c(4,5,5))
for(i in 1:4)   contribution_prec[i,,] <- t(B[i,,]) %*% y_prec[i,,] %*% B[i,,]
prior_prec <- solve(prior_cov)
post_prec <- prior_prec + apply(contribution_prec, c(2,3), sum)
post_cov <- solve(post_prec)
post_mean <- (post_cov %*% (prior_prec %*% prior_mean + apply(contribution_mean, 2, sum)))[,1]

Cmcmc$run(100000)
Csamples <- as.matrix(Cmcmc$mvSamples)

dif_mean <- as.numeric(apply(Csamples, 2, mean)) - post_mean
test_that('posterior mean', expect_lt(max(abs(dif_mean)), 0.001))

dif_cov <- as.numeric(cov(Csamples) - post_cov)
test_that('posterior cov', expect_lt(max(abs(dif_cov)), 0.001))



### test of conjugate Wishart

set.seed(0)

trueCor <- matrix(c(1, .3, .7, .3, 1, -0.2, .7, -0.2, 1), 3)
covs <- c(3, 2, .5)

trueCov = diag(sqrt(covs)) %*% trueCor %*% diag(sqrt(covs))
Omega = solve(trueCov)

n = 20
R = diag(rep(1,3))
mu = 1:3
Y = mu + t(chol(trueCov)) %*% matrix(rnorm(3*n), ncol = n)
M = 3
data <- list(Y = t(Y), n = n, M = M, mu = mu, R = R)

code <- nimbleCode( {
  for(i in 1:n) {
    Y[i, 1:M] ~ dmnorm(mu[1:M], Omega[1:M,1:M]);
  }
  Omega[1:M,1:M] ~ dwish(R[1:M,1:M], 4);
})

newDf = 4 + n
newR = R + tcrossprod(Y- mu)
OmegaTrueMean = newDf * solve(newR)

wishRV <- array(0, c(M, M, 10000))
for(i in 1:10000) {
  z <- t(chol(solve(newR))) %*% matrix(rnorm(3*newDf), ncol = newDf)
  wishRV[ , , i] <- tcrossprod(z)
}
OmegaSimTrueSDs = apply(wishRV, c(1,2), sd)

test_mcmc(model = code, name = 'conjugate Wishart', data = data, seed = 0, numItsC = 1000, inits = list(Omega = OmegaTrueMean),
          results = list(mean = list(Omega = OmegaTrueMean ),
            sd = list(Omega = OmegaSimTrueSDs)),
          resultsTolerance = list(mean = list(Omega = matrix(.05, M,M)),
            sd = list(Omega = matrix(0.06, M, M))))
# issue with Chol in R MCMC - probably same issue as in jaw-linear



## testing conjugate MVN updating with ragged dependencies;
## that is, dmnorm dependents of different lengths from the target node
code <- nimbleCode({
    x[1:3] ~ dmnorm(mu0[1:3], prec = ident[1:3,1:3])
    mu_y2[1:2] <- asCol(a[1:2]) + B[1:2,1:3] %*% asCol(x[1:3])
    mu_y3[1:3] <- asCol(a[1:3]) + B[1:3,1:3] %*% asCol(x[1:3])
    mu_y5[1:5] <- asCol(a[1:5]) + B[1:5,1:3] %*% asCol(x[1:3])
    y2[1:2] ~ dmnorm(mu_y2[1:2], prec = prec_y[1:2,1:2])
    y3[1:3] ~ dmnorm(mu_y3[1:3], prec = prec_y[1:3,1:3])
    y5[1:5] ~ dmnorm(mu_y5[1:5], prec = prec_y[1:5,1:5])
})

mu0 <- rep(0,3)
ident <- diag(3)
a <- 11:15
B <- matrix(1:15, nrow=5, ncol=3, byrow=TRUE)
prec_y <- diag(1:5)

constants <- list(mu0=mu0, ident=ident, a=a, B=B, prec_y=prec_y)
data <- list(y2=1:2, y3=1:3, y5=1:5)
inits <- list(x=rep(0,3))

Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
Rmcmc <- buildMCMC(spec)

Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(10)

set.seed(0)
Cmcmc$run(10)

Rsamples <- as.matrix(Rmcmc$mvSamples)
Csamples <- as.matrix(Cmcmc$mvSamples)

test_that('correct samples for ragged dmnorm conjugate update', expect_true(all(abs(as.numeric(Rsamples[,]) - c(4.96686874, 3.94112676, 4.55975130, 4.01930176, 4.47744412, 4.12927167, 4.91242131, 4.62837537, 4.54227859, 4.97237602, -1.12524733, 1.24545265, -0.13454814, 0.82755276, 0.08252775, 0.71187071, -0.31322184, -0.57462284, -0.64800963, -0.52885823, -3.92276916, -5.23904995, -4.53535941, -4.89919931, -4.66995650, -4.94181562, -4.63558011, -4.16385294, -4.03469945, -4.51128205)) < 1E-8)))

dif <- Rsamples - Csamples

test_that('R and C samples same for ragged dmnorm conjugate update', expect_true(all(abs(dif) < 2E-13)))

set.seed(0)
Cmcmc$run(200000)
Csamples <- as.matrix(Cmcmc$mvSamples)

obsmean <- apply(Csamples, 2, mean)

obsprec <- inverse(cov(Csamples))

pprec <- ident +
    t(B[1:2,1:3]) %*% prec_y[1:2,1:2] %*% B[1:2,1:3] +
    t(B[1:3,1:3]) %*% prec_y[1:3,1:3] %*% B[1:3,1:3] +
    t(B[1:5,1:3]) %*% prec_y[1:5,1:5] %*% B[1:5,1:3]


pmean <- inverse(pprec) %*% (ident %*% mu0 +
             t(B[1:2,1:3]) %*% prec_y[1:2,1:2] %*% (1:2 - a[1:2]) +
             t(B[1:3,1:3]) %*% prec_y[1:3,1:3] %*% (1:3 - a[1:3]) +
             t(B[1:5,1:3]) %*% prec_y[1:5,1:5] %*% (1:5 - a[1:5])   )


test_that('ragged dmnorm conjugate posterior mean', expect_true(all(abs(pmean - obsmean) / pmean < 0.01)))


test_that('ragged dmnorm conjugate posterior precision', expect_true(all(abs(pprec - obsprec) / pprec < 0.005)))



## testing binary sampler

code <- nimbleCode({
    a ~ dbern(0.5)
    b ~ dbern(0.6)
    c ~ dbern(0.05)
    d ~ dbin(prob=0.2, size=1)
    e ~ dbinom(prob=0.9, size=1)
    f ~ dbern(0.5)
    g ~ dbern(0.5)
    h ~ dbern(0.5)
    for(i in 1:10)
        yf[i] ~ dnorm(f, sd = 1)
    for(i in 1:10)
        yg[i] ~ dnorm(g, sd = 1)
    for(i in 1:10)
        yh[i] ~ dnorm(h, sd = 1)
})
constants <- list()
data <- list(yf = c(rep(0,2), rep(1,8)), yg = c(rep(0,8), rep(1,2)), yh = c(rep(0,5), rep(1,5)))
inits <- list(a=0, b=0, c=0, d=0, e=0, f=0, g=0, h=0)

Rmodel <- nimbleModel(code, constants, data, inits)

test_that('model$isBinary', expect_true(Rmodel$isBinary('a')))
test_that('model$isBinary', expect_true(Rmodel$isBinary('b')))
test_that('model$isBinary', expect_true(Rmodel$isBinary('c')))
test_that('model$isBinary', expect_true(Rmodel$isBinary('d')))
test_that('model$isBinary', expect_true(Rmodel$isBinary('e')))
test_that('model$isBinary', expect_true(Rmodel$isBinary('f')))
test_that('model$isBinary', expect_true(Rmodel$isBinary('g')))
test_that('model$isBinary', expect_true(Rmodel$isBinary('h')))

spec <- configureMCMC(Rmodel, nodes = NULL)
spec$addSampler('a', 'binary', print=FALSE)
spec$addSampler('b', 'binary', print=FALSE)
spec$addSampler('c', 'binary', print=FALSE)
spec$addSampler('d', 'binary', print=FALSE)
spec$addSampler('e', 'binary', print=FALSE)
spec$addSampler('f', 'binary', print=FALSE)
spec$addSampler('g', 'binary', print=FALSE)
spec$addSampler('h', 'binary', print=FALSE)
##spec$printSamplers()

Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Cmcmc$run(100000)
samples <- as.matrix(Cmcmc$mvSamples)
means <- apply(samples, 2, mean)
##means

tol <- 0.0025
test_that('binary sampler posterior', expect_lt(abs(means[['a']] - 0.5), tol))
test_that('binary sampler posterior', expect_lt(abs(means[['b']] - 0.6), tol))
test_that('binary sampler posterior', expect_lt(abs(means[['c']] - 0.05), tol))
test_that('binary sampler posterior', expect_lt(abs(means[['d']] - 0.2), tol))
test_that('binary sampler posterior', expect_lt(abs(means[['e']] - 0.9), tol))
test_that('binary sampler posterior', expect_lt(abs(means[['f']] - 0.9525), tol))
test_that('binary sampler posterior', expect_lt(abs(means[['g']] - 0.0475), tol))
test_that('binary sampler posterior', expect_lt(abs(means[['h']] - 0.5), tol))



## testing the binary sampler handles 'out of bounds' ok

code <- nimbleCode({
    px ~ dbern(0.5)
    py ~ dbern(0.5)
    x ~ dnorm(0, sd = px - 0.5)
    y ~ dnorm(0, tau = py)
})
constants <- list()
data <- list(x = 0, y = 0)
inits <- list(px = 1, py = 1)
Rmodel <- nimbleModel(code, constants, data, inits)

spec <- configureMCMC(Rmodel)
spec$printSamplers()
Rmcmc <- buildMCMC(spec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

set.seed(0)
Rmcmc$run(100)
Rsamples <- as.matrix(Rmcmc$mvSamples)
test_that('binary sampler out-of-bounds', expect_true(all(as.numeric(Rsamples) == 1)))

set.seed(0)
Cmcmc$run(100)
Csamples <- as.matrix(Cmcmc$mvSamples)
test_that('binary sampler out-of-bounds', expect_true(all(as.numeric(Csamples) == 1)))

test_that('MCMC with logProb variable being monitored builds and compiles.', {
  cat('===== Starting MCMC test of logProb variable monitoring =====')
  code <- nimbleCode({
    prob[1] <- p
    prob[2] <- 1-p
    x[1:2] ~ dmultinom(size = N, prob = prob[1:2])
    y      ~ dbinom(   size = N, prob = p)
  })
  set.seed(0)
  N <- 100
  p <- 0.3
  x1 <- rbinom(1, size=N, prob=p)
  x2 <- N - x1
  inits <- list(N = N, p = p, x = c(x1, x2), y = x1)
  Rmodel <- nimbleModel(code, constants=list(), data=list(), inits=inits)
  Cmodel <- compileNimble(Rmodel)
  conf <- configureMCMC(Rmodel, monitors = 'logProb_y')
  Rmcmc  <- buildMCMC(conf)
  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
  Cmcmc$run(10)
})


