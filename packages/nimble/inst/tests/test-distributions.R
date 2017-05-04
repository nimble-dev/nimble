## File for testing distributions provided by NIMBLE

source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
context('Testing NIMBLE distributions')

## mvt
x <- c(1, 1, 2)
mn <- c(1, 2, 3)
sc <- diag(c(1, 2, 3))
sc[1, 3] <- sc[3, 1] <- 0.1
df <- 5

## test use directly from R

truth <- mvtnorm::dmvt(x, delta = mn, sigma = sc, df = df, log = FALSE)

try(test_that("dmvt_chol calculates density correctly in R: ",
              expect_equal(dmvt_chol(x, mn, chol(sc), df, prec_param = FALSE),
                           truth,
                           info = paste0("incorrect dmvt calculation in R"))))

## test use through nimble function

nf <- nimbleFunction(
                     run = function(x = double(1), mn = double(1),
                       scale = double(2), df = double(0)) {
                       returnType(double(0))
                       ch <- chol(scale)
                       out <- dmvt_chol(x = x, mu = mn, cholesky = ch,
                                        df = df, prec_param = FALSE, log = FALSE)
                       return(out)
                     }
                     )

cnf <- compileNimble(nf)

try(test_that("Test that dmvt_chol works correctly in R nimble function: ",
              expect_equal(nf(x, mn, sc, df), (truth), 
                           info = paste0("incorrect dmvt value in R nimble function"))))

try(test_that("Test that dmvt_chol works correctly in compiled nimble function: ",
              expect_equal(cnf(x, mn, sc, df), (truth), 
                           info = paste0("incorrect dmvt value in compiled nimble function"))))

## test use in model

mvt_code <- nimbleCode({
  x[1:3] ~ dmvt(mn[1:3], scale = sc[1:3,1:3], df = df)
})

mvt_model <- nimbleModel(mvt_code, constants = list(mn = mn, sc = sc, prec = FALSE, df = df))

mvt_model$x <- x

try(test_that("Test that dmvt calculation is correct in model likelihood calculation: ",
              expect_equal(exp(mvt_model$calculate()), (truth),
                           info = paste0("incorrect likelihood value for dmvt"))))

c_mvt_model <- compileNimble(mvt_model)
c_mvt_model$x
try(test_that("Test that dmvt (compiled) calculation is correct in model likelihood calculation: ",
              expect_equal(exp(c_mvt_model$calculate()), (truth),
                           info = paste0("incorrect likelihood value for dmvt (compiled)"))))

## random sampling
# reference_samps <- mvtnorm::rmvt(n = 10000, delta = mn, sigma = sc, df = df)
r_samps <- t(replicate(10000, rmvt_chol(n = 1, mn, chol(sc), df, prec_param = FALSE)))
true_cov <- sc*df/(df-2)


try(test_that("Test that random samples (R) have correct mean: ",
              expect_equal(colMeans(r_samps), mn, 
                           tol = 0.03,
                           info = "Difference in means exceeds tolerance")))

try(test_that("Test that random samples (R) have correct covariance: ",
              expect_equal(cov(r_samps), true_cov,
                           tol = 0.1,
                           info = "Difference in covs exceeds tolerance")))


nf_sampling <- nimbleFunction(
                              run = function(mn = double(1), scale = double(2), df = double(0)) {
                                returnType(double(1))
                                ch <- chol(scale)
                                out <- rmvt_chol(n = 1, mu = mn, cholesky = ch,
                                                 df = df, prec_param = FALSE)
                                return(out)
                              }
                              )

nf_samps <- t(replicate(10000, nf_sampling(mn, sc, df)))

try(test_that("Test that random samples (nf) have correct mean: ",
              expect_equal(colMeans(nf_samps), mn, 
                           tol = 0.03,
                           info = "Difference in means exceeds tolerance")))

try(test_that("Test that random samples (nf) have correct covariance: ",
              expect_equal(cov(nf_samps), true_cov, 
                           tol = 0.1,
                           info = "Difference in covs exceeds tolerance")))

## sampling via `simulate`
simul_samp <- function(model) {
  model$simulate()
  return(model$x)
}

simul_samps <- t(replicate(10000, simul_samp(c_mvt_model)))

try(test_that("Test that random samples (simulate) have correct mean: ",
              expect_equal(colMeans(simul_samps), mn,
                           tol = 0.03,
                           info = "Difference in means exceeds tolerance")))

try(test_that("Test that random samples (simulate) have correct covariance: ",
              expect_equal(cov(simul_samps), true_cov, 
                           tol = 0.1,
                           info = "Difference in covs exceeds tolerance")))


## dmulti and dcat

set.seed(0)
normGen <- rmulti(1, 1000, prob = c(.1, .1, .8))
set.seed(0)
unnormGen <- rmulti(1, 1000, prob = c(10, 10, 80))
normResult <- dmulti(normGen, prob = c(.1, .1, .8), log = TRUE)
unnormResult <- dmulti(normGen, prob = c(10, 10, 80), log = TRUE)

try(test_that("rmulti handles 'probs' that do not sum to one: ",
              expect_identical(normGen, unnormGen,
                               info = "normalized and unnormalized probabilities give different results")))
try(test_that("dmulti handles 'probs' that do not sum to one: ",
              expect_equal(normResult, unnormResult,
                           info = "normalized and unnormalized probabilities give different results")))

set.seed(0)
normGen <- rcat(1, prob = c(.1, .1, .8))
set.seed(0)
unnormGen <- rcat(1, prob = c(10, 10, 80))
normResult <- dcat(normGen, prob = c(.1, .1, .8), log = TRUE)
unnormResult <- dcat(normGen, prob = c(10, 10, 80), log = TRUE)

try(test_that("rcat handles 'probs' that do not sum to one: ",
              expect_identical(normGen, unnormGen,
                               info = "normalized and unnormalized probabilities give different results")))
try(test_that("dcat handles 'probs' that do not sum to one: ",
              expect_equal(normResult, unnormResult,
                           info = "normalized and unnormalized probabilities give different results")))

## dinvgamma

y <- 1.1; a <- 1; c <- 2; alpha <- 3; beta <- 2; theta <- 1

manDens <- alpha*log(beta) - lgamma(alpha) - (alpha+1)*log(y) - beta/y
try(test_that("Test that dinvgamma gets correct result: ",
              expect_equal(manDens, dinvgamma(y, alpha, scale = beta, log = TRUE), tol = 1e-15)))

set.seed(0)
smp1 <- rinvgamma(100000, shape = alpha, scale = beta)
set.seed(0)
smp2 <- rinvgamma(100000, shape = alpha, rate = 1/beta)

try(test_that("Test that rinvgamma with scale gets correct result: ",
              expect_equal(beta / (alpha-1), mean(smp1), tol = 0.01,
                           info = "Difference in mean exceeds tolerance")))
try(test_that("Test that rinvgamma with rate gets correct result: ",
              expect_equal(beta / (alpha-1), mean(smp2), tol = 0.01,
                           info = "Difference in mean exceeds tolerance")))
try(test_that("Test that rinvgamma with scale gets correct result: ",
              expect_equal(beta/((alpha-1)*sqrt(alpha-2)), sd(smp1), tol = 0.1,
                           info = "Difference in sd exceeds tolerance")))
try(test_that("Test that rinvgamma with rate gets correct result: ",
              expect_equal(beta/((alpha-1)*sqrt(alpha-2)), sd(smp2), tol = 0.1,
                           info = "Difference in sd exceeds tolerance")))

quantile <- quantile(smp1, .15)
attributes(quantile) <- NULL
try(test_that("Test that pinvgamma gets correct result: ",
              expect_equal(qinvgamma(.15, alpha, scale = beta), quantile, tol = 0.005,
                           info = "Difference in quantile exceeds tolerance")))
p <- mean(smp1 < .5)
try(test_that("Test that qinvgamma gets correct result: ",
              expect_equal(pinvgamma(.5, alpha, scale = beta), p, tol = 0.005,
                           info = "Difference in probability exceeds tolerance")))


code <- nimbleCode({
  y ~ dinvgamma(a, scale = c*theta)
  theta ~ dinvgamma(alpha, beta)
})
m = nimbleModel(code, data = list(y = y),
  inits = list(theta = theta, a = a, c = c, alpha = alpha, beta = beta))
conf <- configureMCMC(m)
samplers <- conf$getSamplers()
try(test_that("dinvgamma-dinvgamma conjugacy with dependency using scale",
              expect_identical(samplers[[1]]$name, 'RW',
                               info = "conjugacy improperly detected")))

code <- nimbleCode({
  y ~ dinvgamma(a, rate = 2 + c*theta)
  theta ~ dinvgamma(alpha, beta)
})
m = nimbleModel(code, data = list(y = y),
  inits = list(theta = theta, a = a, c = c, alpha = alpha, beta = beta))
conf <- configureMCMC(m)
samplers <- conf$getSamplers()
try(test_that("dinvgamma-dinvgamma conjugacy with linear dependency",
              expect_identical(samplers[[1]]$name, 'RW',
                               info = "conjugacy improperly detected")))

code <- nimbleCode({
  y ~ dinvgamma(a, rate = c*theta)
  theta ~ dinvgamma(alpha, beta)
})
m = nimbleModel(code, data = list(y=y),
  inits = list(theta = theta, a = a, c = c, alpha = alpha, beta = beta))
conf <- configureMCMC(m)
samplers <- conf$getSamplers()
try(test_that("dinvgamma-dinvgamma conjugacy with dependency using rate",
              expect_identical(samplers[[1]]$name, 'conjugate_dinvgamma_dinvgamma',
                               info = "conjugacy not detected")))
mcmc <- buildMCMC(conf)
comp <- compileNimble(m, mcmc)
set.seed(0)
comp$mcmc$run(100)
smp <- as.matrix(comp$mcmc$mvSamples)

manualSampler <- function(n, y, a, c, alpha, beta) {
  out <- rep(0, n)
  shape = a + alpha
  scale = beta + 1/(c*y)
  set.seed(0)
  out1 <- 1/rgamma(n, shape, rate = scale)
  set.seed(0)
  out2 <- rinvgamma(n, shape, scale = scale)
  return(list(out1, out2))
}
smpMan <- manualSampler(100, y, a, c, alpha, beta) 

try(test_that("Test that invgamma conjugate sampler gets correct result: ",
              expect_identical(smp[,1], smpMan[[1]],
                               info = "NIMBLE conjugate sampler and manual sampler results differ")))

try(test_that("Test that invgamma conjugate sampler gets correct result: ",
              expect_identical(smp[,1], smpMan[[2]],
                               info = "NIMBLE conjugate sampler and manual sampler results differ")))

code <- nimbleCode({
  y ~ dinvgamma(a, scale = c*theta)
  theta ~ dgamma(alpha, beta)
})
m = nimbleModel(code, data = list(y=y),
  inits = list(theta = theta, a = a, c = c, alpha = alpha, beta = beta))
conf <- configureMCMC(m)
samplers <- conf$getSamplers()
try(test_that("dgamma-dinvgamma conjugacy with dependency using rate",
              expect_identical(samplers[[1]]$name, 'conjugate_dgamma_dinvgamma',
                               info = "conjugacy not detected")))
mcmc <- buildMCMC(conf)
comp <- compileNimble(m, mcmc)
set.seed(0)
comp$mcmc$run(10)
smp <- as.matrix(comp$mcmc$mvSamples)

manualSampler <- function(n, y, a, c, alpha, beta) {
  out <- rep(0, n)
  shape = a + alpha
  rate = beta + c/y
  set.seed(0)
  out <- rgamma(n, shape, rate = rate)
  return(out)
}
smpMan <- manualSampler(10, y, a, c, alpha, beta)

try(test_that("Test that gamma conjugate sampler with invgamma dependency gets correct result: ",
              expect_identical(smp[,1], smpMan,
                               info = "NIMBLE gamma conjugate sampler and manual sampler results differ")))

# dinvwish_chol

# tests:
# draw rinvwish_chol, check vs. mean
# check dinvwish_chol vs. direct coding
# simple conjugacy where we have y ~ dmnorm(0, W), W ~ invwish(); compare posterior mean to known posterior mean
# and test each using rate instead of scale

set.seed(1)
df <- 20
d <- 3
C <- crossprod(matrix(rnorm(d^2, 3), d))
pm <- C / (df - d - 1)

n <- 100
draws1 <- draws2 <- draws3 <- array(0, c(d, d, n))
set.seed(1)
for(i in 1:n)
  draws1[,,i] <- solve(rwish_chol(1, chol(C), df, scale_param = FALSE))
pmean1 <- apply(draws1, c(1,2), mean)

for(i in 1:n)
  draws2[,,i] <- rinvwish_chol(1, chol(C), df, scale_param = TRUE)
pmean2 <- apply(draws2, c(1,2), mean)

for(i in 1:n)
  draws3[,,i] <- rinvwish_chol(1, solve(chol(C)), df, scale_param = FALSE)
pmean3 <- apply(draws3, c(1,2), mean)

try(test_that("Test that rinvwish_chol and rwish_chol give correct results: ", {
              expect_equal(max(abs(pmean1 - pm)), 0, tol = 0.03,
                          info = "mean of inverse of rwish draws differs from truth"),
              expect_equal(max(abs(pmean2 - pm)), 0, tol = 0.03,
                          info = "mean of rinvwish with scale draws differs from truth"),
              expect_equal(max(abs(pmean3 - pm)), 0, tol = 0.03,
                          info = "mean of rinvwish with rate differs from truth")}))

dens1 <- dinvwish_chol(draw1[,,1], chol(C), df, log = TRUE)
dens2 <- dinvwish_chol(draw, chol(solve(C)), df, log = TRUE, scale_param = FALSE)

dfun <- function(W, S, nu) {
  k <- nrow(W)
  U = chol(S)
  return(-(log(2)*nu*k/2+(k*(k-1)/4)*log(pi) +sum(lgamma((nu + 1 - 1:k)/2))) + nu*sum(log(diag(U))) -
         (nu+k+1)*sum(log(diag(chol(W)))) -0.5*sum(diag(S %*% solve(W))))
}
                  
dens3 <-  dfun(draw[,,1], C, df)

try(test_that("Test that dinvwish_chol gives correct results: ", {
              expect_equal(dens1-dens3, 0, tol = 0.000001,
                          info = "dinvwish with scale differs from truth"),
              expect_equal(dens2-dens3, 0, tol = 0.000001,
                          info = "dinvwish with rate differs from truth")}))
                                                                        
set.seed(1)

trueCor <- matrix(c(1, .3, .7, .3, 1, -0.2, .7, -0.2, 1), 3)
covs <- c(3, 2, .5)

trueCov = diag(sqrt(covs)) %*% trueCor %*% diag(sqrt(covs))

n = 20
S = diag(rep(1,3))
nu = 4                                    
mu = 1:3
Y = mu + t(chol(trueCov)) %*% matrix(rnorm(3*n), ncol = n)
M = 3
data <- list(Y = t(Y), n = n, M = M)

code <- nimbleCode( {
  for(i in 1:n) {
    Y[i, 1:M] ~ dmnorm(mu[1:M], cov = Omega[1:M,1:M])
  }
  invOmega[1:M,1:M] <- inverse(Omega[1:M,1:M])
  Omega[1:M,1:M] ~ dinvwish(S[1:M,1:M], nu)
})

newDf = nu + n
newS = S + tcrossprod(Y- mu)
OmegaTrueMean = newS / (newDf - M - 1)
invOmegaTrueMean = newDf * solve(newS)
                   
invwishRV <- array(0, c(M, M, 10000))
for(i in 1:10000) {
  invwishRV[,,i] <- rinvwish_chol(1, chol(newS), df = newDf, scale_param = TRUE)
}
OmegaSimTrueSDs = apply(invwishRV, c(1,2), sd)

test_mcmc(model = code, name = 'conjugate Wishart', data = data, seed = 0, numItsC = 1000, inits = list(Omega = trueCov),
          results = list(mean = list(Omega = OmegaTrueMean, invOmega = invOmegaTrueMean, mu = mu, S = S),
            sd = list(Omega = OmegaSimTrueSDs)),
          resultsTolerance = list(mean = list(Omega = matrix(.05, M,M), invOmega = matrix(.05, M, M)),
            sd = list(Omega = matrix(0.06, M, M))))

data = list(Y=t(Y))
constants = list(n = n, M = M)
m = nimbleModel(code, data = data, inits = list(Omega = trueCov, mu = mu, S = S),
  constants = constants)
conf <- configureMCMC(m)
conf$getSamplers()
mcmc <- buildMCMC(conf)
cc = compileNimble(mcmc, m)
cc$mcmc$run(1000)
smp <- as.matrix(cc$mcmc$mvSamples)
pm <- colMeans(smp)
psd <- apply(smp, 2, sd)

try(test_that("dinvwish-dmnorm conjugacy",
              expect_identical(samplers[[1]]$name, 'conjugate_dinvwish_dmnorm',
                               info = "conjugacy not detected")))
