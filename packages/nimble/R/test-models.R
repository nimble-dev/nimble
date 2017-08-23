source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

allModels <- c(# vol1
               'blocker', 'bones', 'dyes', 'equiv', 'line', 'oxford', 'pump', 'rats',
               # vol2
               'dugongs')

## Until NIMBLE 0.6-6, providing an unnecessary variable in data
## generated an error.  Now it generates a warning.  This test
## checks that unwanted data does not generate an error.
test_that('unnecessary data do not break model building', {
    data <- list(a = 1, unwantedVariable = 2)

    ## This use of try() within a test is done carefully.
    ## The following expectations determine if an error was thrown.
    ## Otherwise we would need additional infrastructure.
    toyModel <-
        try(
            nimbleModel(
                nimbleCode({
                    a ~ dnorm(0,1)
                }),
                data = 123
            )
        )

    expect_false(inherits(toyModel, 'try-error'),
                 'nimbleModel stopped due to unnecessary data.')
    expect_true(inherits(toyModel, 'modelBaseClass'),
                'nimbleModel turned out wrong.')
})

sapply(allModels, testBUGSmodel, useInits = TRUE)

# special cases in vol1: 'epil', 'leuk', 'salm', 'seeds'

# Problem cases
# dinterval, I(), T(), random indexing:
# kidney, litters, lsat, mice

# data preparation issue in data block: inhaler (has (k+):(k) style indexing

# various cases where we need to refer to a differently-named .bug file:

testBUGSmodel('epil', model = 'epil2.bug', inits = 'epil-inits.R',
              data = 'epil-data.R', useInits = TRUE)
testBUGSmodel('epil', model = 'epil3.bug', inits = 'epil-inits.R',
              data = 'epil-data.R', useInits = TRUE)
testBUGSmodel('seeds', model = 'seedsuni.bug', inits = 'seeds-init.R',
              data = 'seeds-data.R', useInits = FALSE)
testBUGSmodel('seeds', model = 'seedssig.bug', inits = 'seeds-init.R',
              data = 'seeds-data.R', useInits = FALSE)
testBUGSmodel('birats', model = 'birats1.bug', inits = 'birats-inits.R',
              data = 'birats-data.R', useInits = TRUE)
testBUGSmodel('birats', model = 'birats3.bug', inits = 'birats-inits.R',
              data = 'birats-data.R', useInits = TRUE)
testBUGSmodel('ice', model = 'icear.bug', inits = 'ice-inits.R',
              data = 'ice-data.R', useInits = TRUE)
testBUGSmodel('beetles', model = 'beetles-logit.bug', inits = 'beetles-inits.R',
              data = 'beetles-data.R', useInits = TRUE)
testBUGSmodel('birats', model = 'birats2.bug', inits = 'birats-inits.R',
              data = 'birats-data.R', useInits = TRUE)

# various cases where we need to modify the BUGS code, generally the indexing

# test of leuk; needs var info on Y and dN since these are used in data block
writeLines(c("var", "Y[N,T],", "dN[N,T];"), con = file.path(tempdir(), "leuk.bug"))
##system(paste0("echo 'var\nY[N,T],\ndN[N,T];' >> ", file.path(tempdir(), "leuk.bug")))
system.in.dir(paste("cat leuk.bug >>", file.path(tempdir(), "leuk.bug")), dir = system.file('classic-bugs','vol1','leuk', package = 'nimble'))
##system(paste("cat", system.file('classic-bugs','vol1','leuk','leuk.bug', package = 'nimble'), ">>", file.path(tempdir(), "leuk.bug")))
# need nimStep in data block as we no longer have step
system.in.dir(paste("sed -i -e 's/step/nimStep/g'", file.path(tempdir(), "leuk.bug")))
##system(paste("sed -i -e 's/step/nimStep/g'", file.path(tempdir(), "leuk.bug")))
testBUGSmodel('leuk', dir = "", model = file.path(tempdir(), "leuk.bug"), data = system.file('classic-bugs','vol1','leuk','leuk-data.R', package = 'nimble'),  inits = system.file('classic-bugs','vol1','leuk','leuk-init.R', package = 'nimble'),  useInits = TRUE)

# salm: need dimensionality of logx
writeLines(c("var","logx[doses];"), con = file.path(tempdir(), "salm.bug"))
##system(paste0("echo 'var\nlogx[doses];' >> ", file.path(tempdir(), "salm.bug")))
system.in.dir(paste("cat salm.bug >>", file.path(tempdir(), "salm.bug")), dir = system.file('classic-bugs','vol1','salm', package = 'nimble'))
##system(paste("cat", system.file('classic-bugs','vol1','salm','salm.bug', package = 'nimble'), ">>", file.path(tempdir(), "salm.bug")))
#system(paste("sed -i -e 's/logx\\[\\]/logx\\[1:doses\\]/g'", file.path(tempdir(), "salm.bug"))) # alternative way to get size info in there
testBUGSmodel('salm', dir = "", model = file.path(tempdir(), "salm.bug"), data = system.file('classic-bugs','vol1','salm','salm-data.R', package = 'nimble'),  inits = system.file('classic-bugs','vol1','salm','salm-init.R', package = 'nimble'),  useInits = TRUE)

file.copy(system.file('classic-bugs','vol2','air','air.bug', package = 'nimble'), file.path(tempdir(), "air.bug"), overwrite=TRUE)
system.in.dir("sed -i -e 's/mean(X)/mean(X\\[\\])/g' air.bug", dir = tempdir())
##system(paste("cat", system.file('classic-bugs','vol2','air','air.bug', package = 'nimble'), ">>", file.path(tempdir(), "air.bug")))
##system(paste("sed -i -e 's/mean(X)/mean(X\\[\\])/g'", file.path(tempdir(), "air.bug")))
testBUGSmodel('air', dir = "", model = file.path(tempdir(), "air.bug"), data = system.file('classic-bugs','vol2','air','air-data.R', package = 'nimble'),  inits = system.file('classic-bugs','vol2','air','air-inits.R', package = 'nimble'),  useInits = TRUE)


# need the next line or gives error:
#1: In replaceConstantsRecurse(x, constEnv, constNames) :
#  Code age was given as known but evaluates to a non-scalar.  This is probably # not what you want.
system.in.dir(paste("sed 's/mean(age)/mean(age\\[1:M\\])/g' jaw-linear.bug > ", file.path(tempdir(), "jaw-linear.bug")), dir = system.file('classic-bugs','vol2','jaw', package = 'nimble')) # alternative way to get size info in there
##system(paste("sed 's/mean(age)/mean(age\\[1:M\\])/g'", system.file('classic-bugs','vol2','jaw','jaw-linear.bug', package = 'nimble'), ">", file.path(tempdir(), "jaw-linear.bug"))) # alternative way to get size info in there
testBUGSmodel('jaw', dir = "", model = file.path(tempdir(), "jaw-linear.bug"), inits = system.file('classic-bugs', 'vol2', 'jaw','jaw-inits.R', package = 'nimble'), data = system.file('classic-bugs', 'vol2', 'jaw','jaw-data.R', package = 'nimble'), useInits = TRUE)


testBUGSmodel('dipper', dir = system.file('classic-bugs', 'other', 'dipper', package = 'nimble'), useInits = FALSE)

# various simple tests of multivariate nodes

# simple test of dmnorm/wishart

K <- 2
y <- c(.1, .3)
model <- function() {
  y[1:K] ~ dmnorm(mu[1:K], prec[1:K,1:K]);
  for(i in 1:K) {
    mu0[i] <- 0
  }
  R[1,1] <- .01
  R[2,2] <- .01
  R[1,2] <- 0
  R[2,1] <- 0
  Omega[1,1] <- .01
  Omega[2,2] <- .01
  Omega[1,2] <- 0
  Omega[2,1] <- 0

  mu[1:K] ~ dmnorm(mu0[1:K], Omega[1:K,1:K])
  prec[1:K,1:K] ~ dwish(R[1:K,1:K], 5)
}

inits <- list(mu = c(0,0), prec = matrix(c(.005,.001,.001,.005), 2))
data <- list(K = K, y = y)

testBUGSmodel(example = 'testN', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)

# test multi/Dirichlet

set.seed(0)
n <- 100
alpha <- c(10, 30, 15, 60, 1)
K <- length(alpha)
if(require(nimble)) {
  p <- rdirch(1, alpha)
  y <- rmulti(1, n, p)
} else {
  p <- c(.12, .24, .10, .53, .01)
  y <- c(rmultinom(1, n, p))
}

model <- function() {
  y[1:K] ~ dmulti(p[1:K], n);
  p[1:K] ~ ddirch(alpha[1:K]);
  for(i in 1:K) {
    log(alpha[i]) ~ dnorm(0, sd = 100);
  }
  # log(alpha) ~ dmnorm(0, .001)
}

inits <- list(p = rep(1/K, K), alpha = rep(1/K, K))
data <- list(n = n, K = K, y = y)

testBUGSmodel(example = 'test', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)

# simple test of dt()

model <- function() {
  y ~ dt(1, mu, 1)
  mu ~ dt(3, 0, 1)
}

inits <- list(mu = 0)
data <- list(y = 3)

testBUGSmodel(example = 'testt', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)

# simple test of negbin

model <- function() {
  y ~ dnegbin(p, n)
  p ~ dbeta(3,3)
}

data = list(n = 10, y = 3)
inits <- list(p = 0.5)

testBUGSmodel(example = 'testnb', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)

# test of multivariate data nodes as rows of 3-d array
# and of multivariate latent nodes as rows of matrix

set.seed(0)
n <- 100
m <- 10
g <- 3
alpha <- c(10, 30, 15, 60, 1)
K <- length(alpha)
p <- matrix(0, g, K)
y <- array(0, c(m, g, K))

if(require(nimble)) {
  for(i in seq_len(g))
    p[i, ] <- rdirch(1, alpha)
} else {
  p[1,]  <- c(.12, .24, .10, .53, .01)
  p[2,]  <- c(.2, .3, .05, .2, .25)
  p[3,]  <- c(.05, .05, .10, .3, .5)
  rmulti <- rmultinom
}

for(i in seq_len(g))
  for(j in seq_len(m))
    y[j, i, ] <- rmulti(1, n, p[i, ])


model <- function() {
  for(i in 1:g)
    for(j in 1:m)
      y[j, i, 1:K] ~ dmulti(p[i, 1:K], n);
  for(i in 1:g)
    p[i, 1:K] ~ ddirch(alpha[1:K]);
  for(i in 1:K) {
    log(alpha[i]) ~ dnorm(0, sd = 100);
  }
}

inits <- list(p = matrix(1/K, g, K), alpha = rep(1/K, K))
data <- list(g =g, m = m, n = n, K = K, y = y)

testBUGSmodel(example = 'test', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)

# test handling of lumped data and constants, and overwriting of
# data by inits

test_that("test of distinguishing lumped data and constants:", {

code <- nimbleCode({
    x ~ dnorm(mu,sig)
    mu ~ dnorm(0, 1)
    y ~ dunif(0,1)
})
xVal <- 0.5
m <- nimbleModel(code, constants = list(sig = 1, x = xVal, y = 0.1))
expect_equal(m$isData('x'), TRUE, info = "'x' not set as data in first test")
expect_equal(c(m$x), xVal, info = "value of 'x' not correctly set in first test")
expect_equal('sig' %in% m$getVarNames(), FALSE, info = "sig not set as constant in first test")

code <- nimbleCode({
    x ~ dnorm(mu,1)
    mu ~ dnorm(0, 1)
})
m <- nimbleModel(code, data = list(x = xVal))
expect_equal(m$isData('x'), TRUE, info = "'x' not set as data in second test")
expect_equal(c(m$x), xVal, info = "value of 'x' not correctly set in second test")

code <- nimbleCode({
    x ~ dnorm(mu,1)
    mu ~ dnorm(0, 1)
})
m <- nimbleModel(code, constants = list(x = xVal))
expect_equal(m$isData('x'), TRUE, info = "'x' not set as data in third test")
expect_equal(c(m$x), xVal, info = "value of 'x' not correctly set in third test")

code <- nimbleCode({
    y[1] ~ dnorm(beta*x[1], 1)
    y[2] ~ dnorm(beta*x[2], 1)
    beta ~ dnorm(0, 1)
 })
m <- nimbleModel(code, data = list(y = c(1, 2), x = c(1, 5)))
expect_equal(m$isData('x'), c(TRUE, TRUE), info = "'x' not set as data in fourth test")
expect_equal('x' %in% m$getVarNames(), TRUE, info = "'x' is not set as variable in fourth test")
expect_equal(sum(c("x[1]", "x[2]") %in% m$getNodeNames()), 0, info = "'x[1]' appears incorrectly in nodes in fourth test")


})

test_that("test of preventing overwriting of data values by inits:", {

code <- nimbleCode({
    x[1] ~ dnorm(mu,1)
    x[2] ~ dnorm(mu,1)
    mu ~ dnorm(0, 1)
})
xVal <- c(3, NA)
xInit <- c(4, 4)
m <- nimbleModel(code, constants = list(x = xVal), inits = list(x = xInit))
try(expect_equal(m$isData('x'), c(TRUE, FALSE), info = "'x' data flag is not set correctly in fourth test"))
try(expect_equal(m$x, c(xVal[1], xInit[2]), info = "value of 'x' not correctly set in fourth test"))
try(expect_equal(c('x[1]','x[2]') %in% m$getNodeNames(), c(TRUE, TRUE), info = "'x' nodes note correctly set in fourth test"))

code <- nimbleCode({
    x[1] ~ dnorm(mu,1)
    x[2] ~ dnorm(mu,1)
    mu ~ dnorm(0, 1)
})
m <- nimbleModel(code, data = list(x = xVal), inits = list(x = xInit))
expect_equal(m$isData('x'), c(TRUE, FALSE), info = "'x' data flag is not set correctly in fifth test")
expect_equal(m$x, c(xVal[1], xInit[2]), info = "value of 'x' not correctly set in fifth test")
expect_equal(c('x[1]','x[2]') %in% m$getNodeNames(), c(TRUE, TRUE), info = "'x' nodes note correctly set in fifth test")

})

test_that("test of the handling of missing covariates:", {

    code <- nimbleCode({
    y[1] ~ dnorm(beta*x[1], 1)
    y[2] ~ dnorm(beta*x[2], 1)
    beta ~ dnorm(0, 1)
})
m <- nimbleModel(code, data = list(y = c(1, 2)), constants = list(x = c(1, NA)))
expect_equal('x' %in% m$getVarNames(), FALSE, info = "'x' is not set as constant in first test")

code <- nimbleCode({
    y[1] ~ dnorm(beta*x[1], 1)
    y[2] ~ dnorm(beta*x[2], 1)
    beta ~ dnorm(0, 1)
})
m <- nimbleModel(code, data = list(y = c(1, 2), x = c(1, NA)))
expect_equal('x' %in% m$getVarNames(), TRUE, info = "'x' is not set as variable in second test")
expect_equal(m$isData('x'), c(TRUE, FALSE), info = "'x' data flags are not set correctly in second test")
expect_equal(sum(c("x[1]", "x[2]") %in% m$getNodeNames()), 0, info = "'x' appears incorrectly in nodes in second test")

code <- nimbleCode({
    y[1] ~ dnorm(beta*x[1], 1)
    y[2] ~ dnorm(beta*x[2], 1)
    beta ~ dnorm(0, 1)
    x[2] ~ dnorm(0, 1)
})
m <- nimbleModel(code, data = list(y = c(1, 2)), constants = list(x = c(1, NA)))
expect_equal('x' %in% m$getVarNames(), TRUE, info = "'x' is not set as variable in second test")
expect_equal(m$isData('x'), c(TRUE, FALSE), info = "'x' data flags are not set correctly in second test")
expect_equal(sum("x[1]" %in% m$getNodeNames()), 0, info = "'x[1]' appears incorrectly in nodes in second test")
expect_equal(sum("x[2]" %in% m$getNodeNames()), 1, info = "'x[2]' does not appears in nodes in second test")

code <- nimbleCode({
    y[1] ~ dnorm(beta*x[1], 1)
    y[2] ~ dnorm(beta*x[2], 1)
    beta ~ dnorm(0, 1)
    x[2] ~ dnorm(0, 1)
})
m <- nimbleModel(code, data = list(y = c(1, 2), x = c(1, NA)))
expect_equal('x' %in% m$getVarNames(), TRUE, info = "'x' is not set as variable in second test")
expect_equal(m$isData('x'), c(TRUE, FALSE), info = "'x' data flags are not set correctly in second test")
expect_equal(sum("x[1]" %in% m$getNodeNames()), 0, info = "'x[1]' appears incorrectly in nodes in second test")
expect_equal(sum("x[2]" %in% m$getNodeNames()), 1, info = "'x[2]' does not appears in nodes in second test")

})


# test of use of alias names for distributions

model <- function() {
  y ~ dnbinom(p, n)
  p ~ dbeta(3,3)
  yalt ~ dnegbin(p, n)

  y2 ~ dbin(p2, n)
  p2 ~ dbeta(3,3)
  yalt2 ~ dbinom(p2, n)
}

data = list(y = 3, yalt = 3, y2 = 3, yalt2 = 3)
constants = list(n = 10)
inits <- list(p = 0.5, p2 = 0.5)

testBUGSmodel(example = 'test_dist_aliases', dir = "",
              model = model, data = c(constants, data), inits = inits,
              useInits = TRUE)

m <- nimbleModel(body(model), constants = constants, inits = inits, check = FALSE)
cm <- compileNimble(m)

set.seed(0)
simulate(m, c('y','y2'))
set.seed(0)
simulate(m, c('yalt','yalt2'))
set.seed(0)
simulate(cm, c('y','y2'))
set.seed(0)
simulate(cm, c('yalt','yalt2'))

try(test_that("dnegbin and dnbinom give same results in R model", expect_equal(
    m$y, m$yalt, info = "simulate gives different results for dnegbin and dnbinom")))
try(test_that("dnegbin and dnbinom give same results in R model", expect_equal(
    getLogProb(m, 'y'), getLogProb(m, 'yalt'), info = "calculate gives different results for dnegbin and dnbinom")))
try(test_that("dbin and dbinom give same results in R model", expect_equal(
    m$y2, m$yalt2, info = "simulate gives different results for dbin and dbinom")))
try(test_that("dbin and dbinom give same results in R model", expect_equal(
    getLogProb(m, 'y2'), getLogProb(m, 'yalt2'), info = "calculate gives different results for dbin and dbinom")))
try(test_that("dnegbin and dnbinom give same results in C model", expect_equal(
    cm$y, cm$yalt, info = "simulate gives different results for dnegbin and dnbinom")))
try(test_that("dnegbin and dnbinom give same results in C model", expect_equal(
    getLogProb(cm, 'y'), getLogProb(cm, 'yalt'), info = "calculate gives different results for dnegbin and dnbinom")))
try(test_that("dbin and dbinom give same results in C model", expect_equal(
    cm$y2, cm$yalt2, info = "simulate gives different results for dbin and dbinom")))
try(test_that("dbin and dbinom give same results in C model", expect_equal(
    getLogProb(cm, 'y2'), getLogProb(cm, 'yalt2'), info = "calculate gives different results for dbin and dbinom")))
