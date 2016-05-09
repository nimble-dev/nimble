source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of default MCMC")

### TODO: add in the special cases for dipper

if(FALSE) { # template for running JAGS for comparison
  require(R2jags)
    dir = system.file(file.path('classic-bugs', 'vol2', 'air'), package = 'nimble')
  data = new.env(); inits = new.env()
  source(file.path(dir, 'air-data.R'), data)
  source(file.path(dir, 'air-init.R'), inits)
  data = as.list(data)
  inits = list(as.list(inits))
  out1 <- jags(data = data, inits = inits,
               parameters.to.save = c('X','theta'), n.chains = 1,
               n.iter = 100000, n.burnin = 50000, n.thin = 1, model.file = file.path(dir, 'air.bug'),
               DIC = FALSE, jags.seed = 0)
  out <- as.mcmc(out1)
}

if(FALSE) {
allModels <- c(# vol1
               'blocker', 'bones', 'dyes', 'equiv', 'line', 'oxford', 'pump', 'rats', 'seeds',
               # 'bones',
               # vol2
               'dugongs')

sapply(allModels, test_mcmc, numItsC = 1000)
}

### Beginning of actual tests

test_mcmc('blocker', numItsC = 1000, resampleData = TRUE)
# 100% coverage; looks fine

test_mcmc('bones', numItsC = 10000, resampleData = TRUE)
# 100% coverage; looks fine

test_mcmc('dyes', numItsC = 1000, resampleData = TRUE)
# 100% coverage; looks fine

test_mcmc('equiv', numItsC = 1000, resampleData = TRUE)
# looks good
# testing: tau[2]=97.95, 198.8 ; tau[1]=102.2,55
# phi = -.008,.052; pi = -.1805,.052

test_mcmc('line', numItsC = 1000, resampleData = TRUE)
# 100% coverage; looks fine

test_mcmc('oxford', numItsC = 1000, resampleData = TRUE)
# probably ok; seems to overcover for 'b', but 'b' in this
# parameteriz'n is a top-level node and the multiplic'n
# by sigma seems to lead to frequentist overcoverage
# similar results in JAGS

test_mcmc('pump', numItsC = 1000, resampleData = TRUE)
# 100% coverage; looks fine

test_mcmc('rats', numItsC = 1000, resampleData = TRUE)
# 93.8% coverage; looks fine and compares well to JAGS
# however in resampleData, one of the taus wildly misses

test_mcmc('seeds', numItsC = 1000, resampleData = TRUE)
# fine

test_mcmc('dugongs', numItsC = 1000, resampleData = TRUE)
# 100% coverage; looks fine


test_mcmc('epil', model = 'epil2.bug', inits = 'epil-inits.R',
              data = 'epil-data.R', numItsC = 1000, resampleData = TRUE)
# looks ok

test_mcmc('epil', model = 'epil3.bug', inits = 'epil-inits.R',
              data = 'epil-data.R', numItsC = 1000, resampleData = TRUE)
# looks ok

test_mcmc('seeds', model = 'seedsuni.bug', inits = 'seeds-init.R',
              data = 'seeds-data.R', numItsC = 1000, resampleData = TRUE)
# looks fine - intervals for b's seem a bit large but probably ok
# particularly since default seeds.bug seems fine
# results compared to JAGS look fine
test_mcmc('seeds', model = 'seedssig.bug', inits = 'seeds-init.R',
              data = 'seeds-data.R', numItsC = 1000, resampleData = TRUE)
# looks fine - intervals for b's seem a bit large but probably ok

test_mcmc('birats', model = 'birats1.bug', inits = 'birats-inits.R',
              data = 'birats-data.R', numItsC = 1000, resampleData = TRUE)
# seems fine

test_mcmc('birats', model = 'birats3.bug', inits = 'birats-inits.R',
              data = 'birats-data.R', numItsC = 1000, resampleData = TRUE)
# seems fine

test_mcmc('birats', model = 'birats2.bug', inits = 'birats-inits.R',
            data = 'birats-data.R', numItsC = 1000, resampleData = TRUE)
# looks fine now that values() returns in order
# result changes as of v0.4 because in v0.3-1 'omega.beta' was found
# as both topNode and nontopNode and was being simulated into
# incorrectly in resampleData - this affected values further downstream

test_mcmc('ice', model = 'icear.bug', inits = 'ice-inits.R',
              data = 'ice-data.R', numItsC = 1000, resampleData = TRUE)
# resampleData gives very large magnitude betas because beta[1],beta[2] are not
# actually topNodes because of (weak) dependence on tau, and
# are simulated from their priors to have large magnitude values

# rework ice example so that beta[1] and beta[2] will be top nodes
system(paste("sed 's/tau\\*1.0E-6/1.0E-6/g'", system.file('classic-bugs','vol2','ice','icear.bug', package = 'nimble'), ">", file.path(tempdir(), "icear.bug"))) 
test_mcmc(model = file.path(tempdir(), "icear.bug"), inits = system.file('classic-bugs', 'vol2', 'ice','ice-inits.R', package = 'nimble'), data = system.file('classic-bugs', 'vol2', 'ice','ice-data.R', package = 'nimble'), numItsC = 1000, resampleData = TRUE)
# looks fine, but alpha and beta values shifted a bit (systematically) relative to JAGS results - on further inspection this is because mixing for this model is poor in both NIMBLE and JAGS - with longer runs they seem to agree (as best as one can tell given the mixing without doing a super long run)

test_mcmc('beetles', model = 'beetles-logit.bug', inits = 'beetles-inits.R',
              data = 'beetles-data.R', numItsC = 1000, resampleData = TRUE)
# getting warning; deterministic model node is NA or NaN in model initialization
# weirdness with llike.sat[8] being NaN on init (actually that makes sense), and with weird lifting of RHS of llike.sat


system(paste0("echo 'var\nY[N,T],\ndN[N,T];' >> ", file.path(tempdir(), "leuk.bug")))
system(paste("cat", system.file('classic-bugs','vol1','leuk','leuk.bug', package = 'nimble'), ">>", file.path(tempdir(), "leuk.bug")))
# need nimStep in data block as we no longer have step
system(paste("sed -i -e 's/step/nimStep/g'", file.path(tempdir(), "leuk.bug"))) 
test_mcmc(model = file.path(tempdir(), "leuk.bug"), name = 'leuk', inits = system.file('classic-bugs', 'vol1', 'leuk','leuk-init.R', package = 'nimble'), data = system.file('classic-bugs', 'vol1', 'leuk','leuk-data.R', package = 'nimble'), numItsC = 1000,
          results = list(mean = list(beta = 1.58), sd = list(beta = 0.43)),
          resultsTolerance = list(mean = list(beta = 0.02), sd = list(beta = 0.02)))

system(paste0("echo 'var\nlogx[doses];' >> ", file.path(tempdir(), "salm.bug"))) 
system(paste("cat", system.file('classic-bugs','vol1','salm','salm.bug', package = 'nimble'), ">>", file.path(tempdir(), "salm.bug")))
test_mcmc(model = file.path(tempdir(), "salm.bug"), name = 'salm', inits = system.file('classic-bugs', 'vol1', 'salm','salm-init.R', package = 'nimble'), data = system.file('classic-bugs', 'vol1', 'salm','salm-data.R', package = 'nimble'), numItsC = 1000)
# looks good compared to JAGS

system(paste("cat", system.file('classic-bugs','vol2','air','air.bug', package = 'nimble'), ">>", file.path(tempdir(), "air.bug")))
system(paste("sed -i -e 's/mean(X)/mean(X\\[\\])/g'", file.path(tempdir(), "air.bug"))) 
test_mcmc(model = file.path(tempdir(), "air.bug"), name = 'air', inits = system.file('classic-bugs', 'vol2', 'air','air-inits.R', package = 'nimble'), data = system.file('classic-bugs', 'vol2', 'air','air-data.R', package = 'nimble'), numItsC = 1000)
# theta[2] posterior is a bit off from JAGS - would be worth more investigation

system(paste("sed 's/mean(age)/mean(age\\[1:M\\])/g'", system.file('classic-bugs','vol2','jaw','jaw-linear.bug', package = 'nimble'), ">", file.path(tempdir(), "jaw-linear.bug"))) # alternative way to get size info in there
test_mcmc(model = file.path(tempdir(), "jaw-linear.bug"), name = 'jaw-linear', inits = system.file('classic-bugs', 'vol2', 'jaw','jaw-inits.R', package = 'nimble'), data = system.file('classic-bugs', 'vol2', 'jaw','jaw-data.R', package = 'nimble'), numItsC = 1000)
# C MCMC runs and seems fine; R MCMC fails as can't do Cholesky of 0 matrix in 2-point method


# vectorized version of jaw to try to deal with scalar/vec bug - not needed now that above works
if(FALSE) {
model <- function() {
  for (i in 1:N) {
     Y[i,1:M] ~ dmnorm(mu[1:M], Omega[1:M,1:M]);  # The 4 measurements for each  
  }                                   # boy are multivariate normal

  mu[1:M] <- beta0 * ones[1:M] + beta1 * (age[1:4] - mean(age[1:4]));
  beta0.uncentred <- beta0 - beta1 * mean(age[1:4]);

  beta0 ~ dnorm(0.0, 0.001); 
  beta1 ~ dnorm(0.0, 0.001); 
  Omega[1:M,1:M] ~ dwish(R[1:M,1:M], 4);	# between-child variance in length at each age	
  #Sigma2[1:M,1:M] <- inverse(Omega[1:M,1:M]);

  for (i in 1:N) {
     for  (j in 1:M) {
        resid[i,j] <- Y[i,j] - mu[j];         # residuals
        resid2[i,j] <- resid[i,j]^2;     # squared residuals
     } 
  }
  RSS <- sum(resid2[1:N,1:M]);                    # Residual Sum of Squares
}
 
inits = list(beta0 = 40, beta1 = 0)
data =list(M=4,N=20, Y = matrix(c(47.8, 46.4, 46.3, 45.1, 47.6, 52.5, 51.2, 49.8, 48.1, 
45, 51.2, 48.5, 52.1, 48.2, 49.6, 50.7, 47.2, 53.3, 46.2, 46.3, 
48.8, 47.3, 46.8, 45.3, 48.5, 53.2, 53, 50, 50.8, 47, 51.4, 49.2, 
52.8, 48.9, 50.4, 51.7, 47.7, 54.6, 47.5, 47.6, 49, 47.7, 47.8, 
46.1, 48.9, 53.3, 54.3, 50.3, 52.3, 47.3, 51.6, 53, 53.7, 49.3, 
51.2, 52.7, 48.4, 55.1, 48.1, 51.3, 49.7, 48.4, 48.5, 47.2, 49.3, 
53.7, 54.5, 52.7, 54.4, 48.3, 51.9, 55.5, 55, 49.8, 51.8, 53.3, 
49.5, 55.3, 48.4, 51.8)  , 20, 4), age = c(8, 8.5, 9, 9.5),
  R = matrix(c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1), 4 ,4),
  ones = rep(1, 4))
test_mcmc(model = model, name = 'dmnorm-dwish example', data = data, inits = inits, numItsC = 1000)
}

  
test_mcmc('pump', resampleData = TRUE, results = list(mean = list(
                    "theta[1]" = 0.06,
                    "theta[2]" = 0.10,
                    "theta[9]" = 1.58,
                    "theta[10]" = 1.97,
                    alpha = 0.73,
                    beta = 0.98)),
          resultsTolerance = list(mean = list(
            "theta[1]" = 0.01,
            "theta[2]" = 0.01,
            "theta[9]" = 0.05,
            "theta[10]" = 0.05,
            alpha = 0.1,
            beta = 0.1)))


## LogProb gap: bug fixed in after v0.3
## Problem that occurred in v0.3: because of gap in logProb_a (i.e. logProb_a[2]
## is defined but logProb_a[1] is not)
## Because logProbs get scrambled, the random walk sampler would always accept, 
## meaning the sd of proposal steps approaches Inf
gapCode <- nimbleCode({
	a[1] <- 1
	a[2] ~ dnorm(0,1)
})

test_mcmc(model = gapCode, seed = 0, numItsC = 100000,
				results = list(mean = list(`a[2]` = 0) ),
				resultsTolerance = list(mean = list(`a[2]` = 0.1)),
				samplers = list(list(type = 'RW', target = 'a[2]'))
				)



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

sampleVals = list(x = c(3.950556165467749, 1.556947815895538, 1.598959152023738, 2.223758981790340, 2.386291653164086, 3.266282048060261, 3.064019155073057, 3.229661999356182, 1.985990552839427, 2.057249437940977),
  c = c( 0.010341199485849559, 0.010341199485849559, 0.003846483017887228, 0.003846483017887228, 0.007257679932131476, 0.009680314740728335, 0.012594777095902964, 0.012594777095902964, 0.018179641351556003, 0.018179641351556003))

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
test_that('R and C equiv', expect_less_than(max(abs(dif)), 1E-15))

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
test_that('posterior mean', expect_true(all(abs(dif_mean) < 0.001)))

dif_cov <- as.numeric(cov(Csamples) - post_cov)
test_that('posterior cov', expect_true(all(abs(dif_cov) < 0.001)))



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

test_that('R and C samples same for ragged dmnorm conjugate update', expect_true(all(abs(dif) < 1E-13)))

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
test_that('binary sampler posterior', expect_less_than(abs(means[['a']] - 0.5), tol))
test_that('binary sampler posterior', expect_less_than(abs(means[['b']] - 0.6), tol))
test_that('binary sampler posterior', expect_less_than(abs(means[['c']] - 0.05), tol))
test_that('binary sampler posterior', expect_less_than(abs(means[['d']] - 0.2), tol))
test_that('binary sampler posterior', expect_less_than(abs(means[['e']] - 0.9), tol))
test_that('binary sampler posterior', expect_less_than(abs(means[['f']] - 0.9525), tol))
test_that('binary sampler posterior', expect_less_than(abs(means[['g']] - 0.0475), tol))
test_that('binary sampler posterior', expect_less_than(abs(means[['h']] - 0.5), tol))



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

