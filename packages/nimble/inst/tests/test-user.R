source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of user-supplied distributions and functions in BUGS code")

## User-supplied functions

dbl <- nimbleFunction(
    run = function(x = double(0)) {
        returnType(double(0))
        return(2*x)
    }
    )
# if not in Global, nimble's scoping can't find it in the
# environment created by testthat
assign('dbl', dbl, envir = .GlobalEnv)
# not working at moment as enclosing env of fun is getting messed up


## two arguments to the nimbleFunction
dblSum <- nimbleFunction(
    run = function(x = double(0), y = double(0)) {
        returnType(double(0))
        return(2*(x+y))
    }
    )
assign('dblSum', dblSum, envir = .GlobalEnv)

## vector input and output
vecdbl <- nimbleFunction(
    run = function(x = double(1)) {
        returnType(double(1))
        return(2*x)
    }
    )
assign('vecdbl', vecdbl, envir = .GlobalEnv)


code <- nimbleCode({
    x ~ dnorm(0, 1)
    dx ~ dnorm(dbl(x), sd = .01)
    y[1:K] ~ dmnorm(vecdbl(mu[1:K]), cov = .0001*I[1:K, 1:K])
    mu[1:K] ~ dmnorm(zeros[1:K], cov = I[1:K, 1:K])
    z ~ dnorm(0, 1)
    dz ~ dnorm(dblSum(x, z), sd = .01)
    # vectorized fun applied to scalar nodes-based variable

    for(i in 1:K) {
        theta[i] ~ dnorm(0, 1)
    }
    w[1:K] ~ dmnorm(vecdbl(theta[1:K]), cov = .0001*I[1:K, 1:K])
})

K <- 3
m <- nimbleModel(code, inits = list(x = 0.25, y = 1:K, mu = 1:K,
                           z = 0.5, theta = rep(.5, K), w = rep(1, K)),
                 constants = list(K = K, zeros = rep(0, K), I = diag(K)))

cm <- compileNimble(m)

set.seed(0)
simulate(m)
set.seed(0)
simulate(cm)

for(var in c('dx', 'y', 'dz', 'w')) {
    try(test_that("Test that R and C models agree with user-supplied functions: ",
                  expect_that(get(var, m), equals(get(var, cm),
                                             info = paste0(var, " values differ")))))
}
try(test_that("Test that values based on user-supplied functions are correct: ",
              expect_that(abs(2*cm$x - cm$dx), is_less_than(.03),
                          info = paste0("x and dx are not consistent"))))
try(test_that("Test that values based on user-supplied functions are correct: ",
              expect_that(max(abs(2*cm$mu - cm$y)), is_less_than(.03),
                          info = paste0("mu and y are not consistent"))))
try(test_that("Test that values based on user-supplied functions are correct: ",
              expect_that(abs(2*(cm$x + cm$z) - cm$dz), is_less_than(.03),
                          info = paste0("x plus z and dz are not consistent"))))
try(test_that("Test that values based on user-supplied functions are correct: ",
              expect_that(max(abs(2*cm$theta - cm$w)), is_less_than(.03),
                          info = paste0("theta and w are not consistent"))))


## User-supplied distributions

dmyexp <- nimbleFunction(
    run = function(x = double(0), rate = double(0), log = integer(0)) {
        returnType(double(0))
        logProb <- log(rate) - x*rate
        if(log) {
            return(logProb)
        } else {
            return(exp(logProb))
        }
    })
assign('dmyexp', dmyexp, envir = .GlobalEnv)


rmyexp <- nimbleFunction(
    run = function(n = integer(0), rate = double(0)) {
        returnType(double(0))
        if(n != 1) nimPrint("rmyexp only allows n = 1; using n = 1.")
        dev <- runif(1, 0, 1)
        return(-log(1-dev) / rate)
    }
    )
assign('rmyexp', rmyexp, envir = .GlobalEnv)

pmyexp <- nimbleFunction(
    run = function(q = double(0), rate = double(0), lower.tail = integer(0), log.p = integer(0)) {
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
assign('pmyexp', pmyexp, envir = .GlobalEnv)

qmyexp <- nimbleFunction(
    run = function(p = double(0), rate = double(0), lower.tail = integer(0), log.p = integer(0)) {
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
assign('qmyexp', qmyexp, envir = .GlobalEnv)


ddirchmulti <- nimbleFunction(
    run = function(x = double(1), alpha = double(1), size = double(0), log = integer(0)) {
        returnType(double(0))
        logProb <- lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum(lgamma(alpha + x)) - lgamma(sum(alpha) + size)
        if(log) {
            return(logProb)
        } else {
            return(exp(logProb))
        }
    })
assign('ddirchmulti', ddirchmulti, envir = .GlobalEnv)

rdirchmulti <- nimbleFunction(
    run = function(n = integer(0), alpha = double(1), size = double(0)) {
        returnType(double(1))
        if(n != 1) nimPrint("rdirchmulti only allows n = 1; using n = 1.")
        p <- rdirch(1, alpha)
        return(rmulti(1, size = size, prob = p))
    })
assign('rdirchmulti', rdirchmulti, envir = .GlobalEnv)

registerDistributions(list(
    dmyexp = list(
        BUGSdist = "dmyexp(rate, scale)",
        Rdist = "dmyexp(rate = 1/scale)",
        altParams = "scale = 1/rate",
        pqAvail = TRUE),
    ddirchmulti = list(
        BUGSdist = "ddirchmulti(alpha, size)",
        types = c('value = double(1)', 'alpha = double(1)', 'size = double(0)'))
        )
    )
                      

code1 <- BUGScode({
    for(i in 1:n1) {
        y1[i] ~ dmyexp(rate = r1)
        y2[i] ~ dmyexp(scale = s2)
    }
    r1 <- 1 / s1
    s1 ~ dunif(0, 100)
    s2 ~ dunif(0, 100)
})

code2 <- BUGScode({
    for(i in 1:n2) {
        y3[i] ~ dpois(lambda)
    }
    lambda ~ T(dmyexp(scale = 5), 0, upper)
})

code3 <- BUGScode({
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

try(test_that("Test that truncation works with user-supplied distribution: ",
              expect_that(max(smp[ , 'lambda']), is_less_than(upper),
                          info = paste0("parameter exceeds upper bound"))))


deregisterDistributions('ddirchmulti')
try(test_that("Test that deregistration of user-supplied distributions works: ",
              expect_that(is.null(nimble:::nimbleUserNamespace$distributions[['ddirchmulti']]), equals(TRUE),
                          info = paste0("ddirchmulti has not been deregistered"))))

