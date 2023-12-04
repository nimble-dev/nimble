source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = TRUE)

context("Testing of NIMBLE model building and operation")

## If you do *not* want to write to results files
##    comment out the sink() call below.  And consider setting verbose = FALSE 
## To record a new gold file, nimbleOptions('generateGoldFileForMCMCtesting') should contain the path to the directory where you want to put it
## e.g. nimbleOptions(generateGoldFileForMCMCtesting = getwd())
## Comparison to the gold file won't work until it is installed with the package.

allModels <- c(# vol1
    'blocker', 'bones', 'dyes', 'equiv', 'line', 'pump', 'rats',
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
                data = data
            )
        )

    expect_false(inherits(toyModel, 'try-error'),
                 'nimbleModel stopped due to unnecessary data.')
    expect_true(inherits(toyModel, 'modelBaseClass'),
                'nimbleModel turned out wrong.')
})

out <- sapply(allModels, testBUGSmodel, useInits = TRUE)

testBUGSmodel('oxford', useInits = TRUE, expectModelWarning = "'tau' has initial values but is not")

## special cases in vol1: 'epil', 'leuk', 'salm', 'seeds'

## Problem cases
## dinterval, I(), T(), random indexing:
## kidney, litters, lsat, mice

## data preparation issue in data block: inhaler (has (k+):(k) style indexing

## various cases where we need to refer to a differently-named .bug file:

testBUGSmodel('epil', model = 'epil2.bug', inits = 'epil-inits.R',
              data = 'epil-data.R', useInits = TRUE, expectModelWarning = "'tau.b' has initial values but is not")
testBUGSmodel('epil', model = 'epil3.bug', inits = 'epil-inits.R',
              data = 'epil-data.R', useInits = TRUE)
testBUGSmodel('seeds', model = 'seedsuni.bug', inits = 'seeds-init.R',
              data = 'seeds-data.R', useInits = FALSE)
testBUGSmodel('seeds', model = 'seedssig.bug', inits = 'seeds-init.R',
              data = 'seeds-data.R', useInits = FALSE)
testBUGSmodel('birats', model = 'birats1.bug', inits = 'birats-inits.R',
              data = 'birats-data.R', useInits = TRUE, expectModelWarning = "'Omega.beta' has initial values but is not")
testBUGSmodel('birats', model = 'birats3.bug', inits = 'birats-inits.R',
              data = 'birats-data.R', useInits = TRUE, expectModelWarning = "'Omega.beta' has initial values but is not")
testBUGSmodel('ice', model = 'icear.bug', inits = 'ice-inits.R',
              data = 'ice-data.R', useInits = TRUE)
testBUGSmodel('beetles', model = 'beetles-logit.bug', inits = 'beetles-inits.R',
              data = 'beetles-data.R', useInits = TRUE)
testBUGSmodel('birats', model = 'birats2.bug', inits = 'birats-inits.R',
              data = 'birats-data.R', useInits = TRUE, expectModelWarning = "'tau.beta' has initial values but is not")

## various cases where we need to modify the BUGS code, generally the indexing

## test of leuk; needs var info on Y and dN since these are used in data block
writeLines(c("var", "Y[N,T],", "dN[N,T];"), con = file.path(tempdir(), "leuk.bug"))
##system(paste0("echo 'var\nY[N,T],\ndN[N,T];' >> ", file.path(tempdir(), "leuk.bug")))
system.in.dir(paste("cat leuk.bug >>", file.path(tempdir(), "leuk.bug")), dir = system.file('classic-bugs','vol1','leuk', package = 'nimble'))
##system(paste("cat", system.file('classic-bugs','vol1','leuk','leuk.bug', package = 'nimble'), ">>", file.path(tempdir(), "leuk.bug")))
## need nimStep in data block as we no longer have step
system.in.dir(paste("sed -i -e 's/step/nimStep/g'", file.path(tempdir(), "leuk.bug")))
##system(paste("sed -i -e 's/step/nimStep/g'", file.path(tempdir(), "leuk.bug")))
testBUGSmodel('leuk', dir = "", model = file.path(tempdir(), "leuk.bug"), data = system.file('classic-bugs','vol1','leuk','leuk-data.R', package = 'nimble'),  inits = system.file('classic-bugs','vol1','leuk','leuk-init.R', package = 'nimble'), useInits = TRUE, expectModelWarning = "'tau' has initial values but is not")

## salm: need dimensionality of logx
writeLines(c("var","logx[doses];"), con = file.path(tempdir(), "salm.bug"))
##system(paste0("echo 'var\nlogx[doses];' >> ", file.path(tempdir(), "salm.bug")))
system.in.dir(paste("cat salm.bug >>", file.path(tempdir(), "salm.bug")), dir = system.file('classic-bugs','vol1','salm', package = 'nimble'))
##system(paste("cat", system.file('classic-bugs','vol1','salm','salm.bug', package = 'nimble'), ">>", file.path(tempdir(), "salm.bug")))
##system(paste("sed -i -e 's/logx\\[\\]/logx\\[1:doses\\]/g'", file.path(tempdir(), "salm.bug"))) ## alternative way to get size info in there
testBUGSmodel('salm', dir = "", model = file.path(tempdir(), "salm.bug"), data = system.file('classic-bugs','vol1','salm','salm-data.R', package = 'nimble'),  inits = system.file('classic-bugs','vol1','salm','salm-init.R', package = 'nimble'),  useInits = TRUE)

out <- file.copy(system.file('classic-bugs','vol2','air','air.bug', package = 'nimble'), file.path(tempdir(), "air.bug"), overwrite=TRUE)
system.in.dir("sed -i -e 's/mean(X)/mean(X\\[\\])/g' air.bug", dir = tempdir())
##system(paste("cat", system.file('classic-bugs','vol2','air','air.bug', package = 'nimble'), ">>", file.path(tempdir(), "air.bug")))
##system(paste("sed -i -e 's/mean(X)/mean(X\\[\\])/g'", file.path(tempdir(), "air.bug")))
testBUGSmodel('air', dir = "", model = file.path(tempdir(), "air.bug"), data = system.file('classic-bugs','vol2','air','air-data.R', package = 'nimble'),  inits = system.file('classic-bugs','vol2','air','air-inits.R', package = 'nimble'),  useInits = TRUE)


## need the next line or gives error:
##1: In replaceConstantsRecurse(x, constEnv, constNames) :
##  Code age was given as known but evaluates to a non-scalar.  This is probably ## not what you want.
system.in.dir(paste("sed 's/mean(age)/mean(age\\[1:M\\])/g' jaw-linear.bug > ", file.path(tempdir(), "jaw-linear.bug")), dir = system.file('classic-bugs','vol2','jaw', package = 'nimble')) ## alternative way to get size info in there
##system(paste("sed 's/mean(age)/mean(age\\[1:M\\])/g'", system.file('classic-bugs','vol2','jaw','jaw-linear.bug', package = 'nimble'), ">", file.path(tempdir(), "jaw-linear.bug"))) ## alternative way to get size info in there
testBUGSmodel('jaw', dir = "", model = file.path(tempdir(), "jaw-linear.bug"), inits = system.file('classic-bugs', 'vol2', 'jaw','jaw-inits.R', package = 'nimble'), data = system.file('classic-bugs', 'vol2', 'jaw','jaw-data.R', package = 'nimble'), useInits = TRUE, expectModelWarning = "'beta2' has initial values but is not")



testBUGSmodel('dipper', dir = system.file('classic-bugs', 'other', 'dipper', package = 'nimble'), useInits = FALSE)

## various simple tests of multivariate nodes

## simple test of dmnorm/wishart

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

### Repeat simple test of dmnorm/wishart with a tweak
## to test that copyIfNeeded operates correctly for a non-complete matrix

K <- 2
y <- c(.1, .3)
model <- function() {
    y[1:K] ~ dmnorm(mu[1:K], prec[1:K,1:K]);
    for(i in 1:(K+1)) {
        mu0[i] <- 0
    }
    R[1:3, 1:3] <- 0.01 * diag(3)
    Omega[1:3, 1:3] <- 0.01 * diag(3)
    cholR[1:3, 1:3] <- chol(R[1:3, 1:3])
    mu[1:K] ~ dmnorm(mu0[1:K], Omega[1:K,1:K])
    ## one shouldn't chop up a cholesky matrix like this,
    ## but it is diagonal and the only purpose here is to test
    ## the copyIfNeeded mechanism
    prec[1:K,1:K] ~ dwish(cholesky = cholR[1:K,1:K], df = 5, scale_param = 1)
}

inits <- list(mu = c(0,0), prec = matrix(c(.005,.001,.001,.005), 2))
data <- list(K = K, y = y)

testBUGSmodel(example = 'testN', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)

## test multi/Dirichlet

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
    ## log(alpha) ~ dmnorm(0, .001)
}

inits <- list(p = rep(1/K, K), alpha = rep(1/K, K))
data <- list(n = n, K = K, y = y)

testBUGSmodel(example = 'test', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)

## Repeat test of multi/Dirichlet using Km1 for length of p (and y).
## This forces additional check in copyIfNeeded
set.seed(0)
n <- 100
alpha <- c(10, 30, 15, 60, 1)
K <- length(alpha)
Km1 <- K-1
if(require(nimble)) {
    p <- rdirch(1, alpha[1:Km1])
    y <- rmulti(1, n, p)
} else {
    p <- c(.12, .24, .10, .53)
    p <- p/sum(p)
    y <- c(rmultinom(1, n, p))
}
model <- function() {
    y[1:Km1] ~ dmulti(p[1:Km1], n);
    p[1:Km1] ~ ddirch(alpha[1:Km1]);
    for(i in 1:K) {
        log(alpha[i]) ~ dnorm(0, sd = 100);
    }
    ## log(alpha) ~ dmnorm(0, .001)
}
inits <- list(p = rep(1/Km1, Km1), alpha = rep(1/K, K))
data <- list(n = n, K = K, y = y, Km1 = Km1)

testBUGSmodel(example = 'test', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)

## simple test of dt()

model <- function() {
    y ~ dt(1, mu, 1)
    mu ~ dt(3, 0, 1)
}

inits <- list(mu = 0)
data <- list(y = 3)

testBUGSmodel(example = 'testt', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)

## simple test of negbin

model <- function() {
    y ~ dnegbin(p, n)
    p ~ dbeta(3,3)
}

data = list(n = 10, y = 3)
inits <- list(p = 0.5)

testBUGSmodel(example = 'testnb', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)

## test of multivariate data nodes as rows of 3-d array
## and of multivariate latent nodes as rows of matrix

set.seed(0)
n <- 100
m <- 10
g <- 3
alpha <- c(10, 30, 15, 60, 1)
K <- length(alpha)
p <- matrix(0, g, K)
y <- array(0, c(m, g, K))

for(i in seq_len(g))
    p[i, ] <- rdirch(1, alpha)
## We had this for when NIMBLE not available but not clear why that would ever be the case.
    ## p[1,]  <- c(.12, .24, .10, .53, .01)
    ## p[2,]  <- c(.2, .3, .05, .2, .25)
    ## p[3,]  <- c(.05, .05, .10, .3, .5)
    ## rmulti <- rmultinom


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

## test handling of lumped data and constants, and overwriting of
## data by inits

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

    expect_message(m <- nimbleModel(code, constants = list(x = xVal), inits = list(x = xInit)), "Ignoring non-NA values in `inits` for data nodes")
    expect_equal(m$isData('x'), c(TRUE, FALSE), info = "'x' data flag is not set correctly in fourth test")
    expect_equal(m$x, c(xVal[1], xInit[2]), info = "value of 'x' not correctly set in fourth test")
    expect_equal(c('x[1]','x[2]') %in% m$getNodeNames(), c(TRUE, TRUE), info = "'x' nodes note correctly set in fourth test")

    code <- nimbleCode({
        x[1] ~ dnorm(mu,1)
        x[2] ~ dnorm(mu,1)
        mu ~ dnorm(0, 1)
    })
    expect_message(m <- nimbleModel(code, data = list(x = xVal), inits = list(x = xInit)), "Ignoring non-NA values in `inits` for data nodes")
    expect_equal(m$isData('x'), c(TRUE, FALSE), info = "'x' data flag is not set correctly in fifth test")
    expect_equal(m$x, c(xVal[1], xInit[2]), info = "value of 'x' not correctly set in fifth test")
    expect_equal(c('x[1]','x[2]') %in% m$getNodeNames(), c(TRUE, TRUE), info = "'x' nodes note correctly set in fifth test")

})

test_that("test of using dimensions of inits when dimension information not available:", {
    code <- nimbleCode({
        for(i in 1:3) {
            y[i] ~ dnorm(mu[k[i]], 1)
            k[i] ~ dcat(p[1:5])
        }
    })
    expect_error(m <- nimbleModel(code, data = list(y = rep(1, 3))), info = "expected error because dimension of mu is unknown")
    m <- nimbleModel(code, data = list(y = rep(1, 3)), inits = list(k = rep(1, 3), mu = 1:5))
    expect_equal(m$modelDef$dimensionsList$mu, 5, info = "dimension for mu not equal to that given in inits")
    expect_message(m <- nimbleModel(code, data = list(y = rep(1, 3)), inits = list(k = rep(1, 3), mu = 1:8), dimensions = list(mu = 5)), "Inconsistent dimensions between inits and dimensions")
})

test_that("test of using dimensions of data when dimension information not available:", {
    code <- nimbleCode({
        for(i in 1:3) {
            y[i] ~ dnorm(mu[k[i]], 1)
            k[i] ~ dcat(p[1:5])
        }
    })
    expect_error(m <- nimbleModel(code, data = list(y = rep(1, 3))), info = "expected error because dimension of mu is unknown")
    m <- nimbleModel(code, data = list(y = rep(1, 3), mu = 1:5), inits = list(k = rep(1, 3)))
    expect_equal(m$modelDef$dimensionsList$mu, 5, info = "dimension for mu not equal to that given in data")
    nimbleOptions(verbose = FALSE)
    expect_error(m <- nimbleModel(code, data = list(y = rep(1, 3), mu = 1:8), inits = list(k = rep(1, 3)), dimensions = list(mu = 5)), info = "expected error because of dimension mismatch")  # error emitted by setData() and warning by assignDimensions()
    nimbleOptions(verbose = TRUE)
    
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

test_that("test of error trapping for indexes that are zero or less:", {
    code <- nimbleCode( {
        for(i in 1:n)
            y[i] ~ dnorm(mu[n-i], 1)})
    expect_error(m <- nimbleModel(code, constants = list(n=3)), 'index value of zero or less')
})

## test of use of alias names for distributions

test_that("use of dbin/dbinom and dnegbin/dnbinom are identical", {
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
    
    expect_equal(m$y, m$yalt, info = "simulate gives different results for dnegbin and dnbinom")
    expect_equal(getLogProb(m, 'y'), getLogProb(m, 'yalt'),
                 info = "calculate gives different results for dnegbin and dnbinom")
    expect_equal(m$y2, m$yalt2, info = "simulate gives different results for dbin and dbinom")
    expect_equal(getLogProb(m, 'y2'), getLogProb(m, 'yalt2'),
                 info = "calculate gives different results for dbin and dbinom")
    expect_equal(cm$y, cm$yalt, info = "compiled simulate gives different results for dnegbin and dnbinom")
    expect_equal(getLogProb(cm, 'y'), getLogProb(cm, 'yalt'),
                 info = "compiled calculate gives different results for dnegbin and dnbinom")
    expect_equal(cm$y2, cm$yalt2, info = "compiled simulate gives different results for dbin and dbinom")
    expect_equal(getLogProb(cm, 'y2'), getLogProb(cm, 'yalt2'),
                 info = "compiled calculate gives different results for dbin and dbinom")
})


test_that("test of using data frame as 'data' in model:", {
    code <- nimbleCode({
        for(i in 1:3) 
            for(j in 1:2) 
                y[i,j] ~ dnorm(0, 1)
        mu ~dnorm(0,1)
    })
    expect_error(m <- nimbleModel(code, data = list(y = data.frame(a = 1:3, b = c('a','b','c')))), info = "expected error because data frame entry to data is non-numeric")
    y <- data.frame(a = rnorm(3), b = rnorm(3))
    m <- nimbleModel(code, data = list(y = y))
    cm <- compileNimble(m)
    y <- as.matrix(y); dimnames(y) <- NULL
    expect_identical(m$y, y, info = "input data frame as data not handled correctly")
    expect_identical(cm$y, y, info = "input data frame as data not handled correctly")
    expect_error(m$setData(list(y = data.frame(a = 1:3, b = c('a','b','c')))))
})

test_that("test of using ragged arrays in a model:", {
    mc <- nimbleCode({
        for(i in 1:2) {
            Z[i, 1:n[i]] <- 2*X[i, 1:n[i]]
        }
    })
    
    n <- c(2, 3)
    X <- matrix(1:6, nrow = 2)
    constants <- list(n = n, X = X)
    nimbleOptions(verbose = FALSE)
    expect_silent(m <- nimbleModel(mc, constants = constants))
    nimbleOptions(verbose = TRUE)
})

test_that("warnings for multiply-defined model nodes:", {
    code <- nimbleCode({
        tmp ~ dnorm(0,1)
        for(i in 1:3) {
            y[i] ~ dnorm(0,1)
            mu ~ dnorm(mu0[i],1)
        }
    })
    expect_warning(m <- nimbleModel(code), "'i' on the left-hand side of 'mu ~ ", fixed = TRUE)
    code <- nimbleCode({
        tmp ~ dnorm(0,1)
        for(i in 1:3) {
            y[i] ~ dnorm(0,1)
            mu ~ dnorm(0,1)
        }
    })
    expect_warning(m <- nimbleModel(code), "'i' on the left-hand side of 'mu ~ ", fixed = TRUE)
    code <- nimbleCode({
        tmp ~ dnorm(0,1)
        for(i in 1:3) {
            for(j in 1:3)
                mu[i+2,1] ~ dnorm(0,1)
        }
    })
    expect_warning(m <- nimbleModel(code), "'j' on the left-hand side of 'mu[i + 2, 1] ~ ", fixed = TRUE)
    code <- nimbleCode({
        tmp ~ dnorm(0,1)
        for(i in 1:3) {
            for(j in 1:3)
                for(k in 1:3)
                mu[i+2,1,1] ~ dnorm(0,1)
        }
    })
    expect_warning(m <- nimbleModel(code), "'j,k' on the left-hand side of 'mu[i + 2, 1, 1] ~ ", fixed = TRUE)
})

test_that("handling of missing indexes of expressions:", {
    code = nimbleCode({
        mn[1:2] <- (X[1:2,1:2] %*% beta[1:2,1:2])[,1]
        y[1:2] ~ dmnorm( mn[1:2], pr[1:2,1:2])
    })
    m = nimbleModel(code, data = list(y = rnorm(2)),
                    inits = list(X = matrix(1, 2,2), beta = matrix(2,2,2), pr = diag(2)))
    cm <- compileNimble(m)
    expect_true(is.numeric(cm$calculate('y')), "incorrectly not dealing with missing index in ()[] expression")

    code = nimbleCode({
        mn[1:2] <- (X[1:2,] %*% beta[1:2,1:2])[,1]
        y[1:2] ~ dmnorm( mn[1:2], pr[1:2,1:2])
    })
    m = nimbleModel(code, data = list(y = rnorm(2)),
                    inits = list(X = matrix(1, 2,2), beta = matrix(2,2,2), pr = diag(2)))
    cm <- compileNimble(m)
    expect_true(is.numeric(cm$calculate('y')), "incorrectly not dealing with missing index in ()[] expression")
    
    code = nimbleCode({
        mn[1:2] <- (X[1:2,] %*% beta[1:2,1:2])[,1]
        y[1:2] ~ dmnorm( mn[1:2], pr[1:2,1:2])
    })
    expect_error(m <- nimbleModel(code, data = list(y = rnorm(2)),
                                  inits = list(beta = matrix(2,2,2), pr = diag(2))),
                 "missing indices", info = "not catching missing indices")

    code = nimbleCode({
    mn[1:2] <- (X[1:2,] %*% beta[1:2,1:2])[,1]
    y[1:2] ~ dmnorm( mn[1:2], pr[1:2,1:2])
    })
    ## Having trouble with consistency in whether output or message is produced,
    ## so just run nimbleModel and test_that should fail if error occurs.
    m <- nimbleModel(code, data = list(y = rnorm(2)),
                    inits = list(beta = matrix(2,2,2), pr = diag(2)),
                    dimensions = list(X = c(2,2)))
                  ## "model building finished",
                  ## info = "incorrectly handling missing indices with dims present")

    code = nimbleCode({
        mn[1:2] <- (X[1:2,1:2] %*% beta[1:2,1:2])[k[,1],1]
        y[1:2] ~ dmnorm( mn[1:2], pr[1:2,1:2])
    })
    expect_error(m <- nimbleModel(code, data = list(y = rnorm(2)),
                                  inits = list(X = matrix(1, 2, 2), beta = matrix(2,2,2), pr = diag(2))),
                 "missing indices", info = "not catching missing indices in model variable within indexing of ()")
    
    code = nimbleCode({
        mn[1:2] <- (X[1:2,1:2] %*% beta[1:2,1:2])[k[,1],1]
        y[1:2] ~ dmnorm( mn[1:2], pr[1:2,1:2])
    })
    m = nimbleModel(code, data = list(y = rnorm(2)),
                    inits = list(X = matrix(1, 2, 2), beta = matrix(2,2,2), pr = diag(2)),
                    dimensions = list(k = c(2,2)))
    cm <- compileNimble(m)  # if compilation fails, test_that should catch this; having trouble using expect_message as behavior of whether a message is detected seems to differ when running tests locally versus Travis.
})

test_that("handling of missing indexes of expressions, part 2:", {
    ## Testing that case like `myfun()[,1]` handled similarly to the above case.
    myfun0 <- nimbleFunction(
    run = function() {
        returnType(double(2))
        out = matrix(3.1, 3, 3)
        return(out)
    })

    myfun1 <- nimbleFunction(
        run = function(x = double(0)) {
            returnType(double(2))
            out = matrix(x, 3, 3)
            return(out)
        })
    
    myfun2 <- nimbleFunction(
        run = function(x = double(1), y = double(0)) {
            returnType(double(2))
            out = matrix(x[1], 3, 3)
            return(out)
        })

    temporarilyAssignInGlobalEnv(myfun0)
    temporarilyAssignInGlobalEnv(myfun1)
    temporarilyAssignInGlobalEnv(myfun2)

    code <- nimbleCode({
        a[1:3] <- myfun0()[,1]      
    })
    m <- nimbleModel(code)
    cm <- compileNimble(m)
    expect_true(is.numeric(cm$calculate('a')))

    code <- nimbleCode({
        a[1:3] <- (myfun0())[,1]      
    })
    m <- nimbleModel(code)
    cm <- compileNimble(m)
    expect_true(is.numeric(cm$calculate('a')))
    
    code <- nimbleCode({
        a[1:3] <- myfun1(b)[,1]      
    })
    m <- nimbleModel(code, inits = list(b = 3.1))
    cm <- compileNimble(m)
    expect_true(is.numeric(cm$calculate('a')))

    code <- nimbleCode({
        a[1:3] <- myfun2(b[1:3, 1], 7)[,1]      
    })
    m <- nimbleModel(code)
    cm <- compileNimble(m)
    expect_true(is.numeric(cm$calculate('a')))

    code <- nimbleCode({
        a[1:3] <- myfun2(b[ , 1], 7)[,1]      
    })
    expect_error(m <- nimbleModel(code), "missing indices")

    code <- nimbleCode({
        a[1:3] <- myfun2(b[1:3, 1], 7)[,1]      
    })
    m <- nimbleModel(code, inits = list(b = matrix(rnorm(9), 3, 3)))
    cm <- compileNimble(m)
    expect_true(is.numeric(cm$calculate('a')))

    code <- nimbleCode({
        a[1:3] <- myfun0()[k[,1],1]      
    })
    m <- nimbleModel(code, inits = list(k=matrix(2, 3,3)))
    cm <- compileNimble(m)
    expect_true(is.numeric(cm$calculate('a')))

    code <- nimbleCode({
        a[1:3] <- myfun0()[k[,1],1]      
    })
    expect_error(m <- nimbleModel(code), "missing indices")
})

test_that("warning when RHS only nodes used as dynamic indexes", {
    code <- nimbleCode({
        for(i in 1:3) 
            y[i] ~ dnorm(mu[k[i]+1], 1)
        for(i in 1:3)
            mu[i] ~ dnorm(0,1)
    })
    
    expect_message(m <- nimbleModel(code, inits = list(k = rep(1,3))),
                   "Detected use of non-constant indexes")
    expect_message(m <- nimbleModel(code, data = list(k = rep(1,3))),
                   "Detected use of non-constant indexes")

    nimbleOptions(verbose = FALSE)
    expect_silent(m <- nimbleModel(code, constants = list(k = rep(1,3))))
    nimbleOptions(verbose = TRUE)
    
    myfun <- nimbleFunction(run = function(x = double()) {
        returnType(double())
        return(1)
    })
    temporarilyAssignInGlobalEnv(myfun)
    
    code <- nimbleCode({
        for(i in 1:3) 
            y[i] ~ dnorm(mu[myfun(k[i])], 1)
        for(i in 1:3)
            mu[i] ~ dnorm(0,1)
    })
    expect_message(m <- nimbleModel(code, inits = list(k = rep(1,3))),
                   "Detected use of non-constant indexes")

    code <- nimbleCode({
        for(i in 1:3)
            y[i] ~ dnorm(mu[k[j[i]]+1], 1)
        for(i in 1:3)
            mu[i] ~ dnorm(0,1)
    })
    expect_message(m <- nimbleModel(code, inits = list(k = rep(1,3)), constants = list(j=1:3)), "Detected use of non-constant indexes")
    expect_message(m <- nimbleModel(code, inits = list(k = rep(1,3), j=1:3)),
                   "Detected use of non-constant indexes")

    ## We don't detect when deterministic node intervenes.
    code <- nimbleCode({
        for(i in 1:3) {
            y[i] ~ dnorm(mu[k[i]+1], 1)
            k[i] <- kk[i] + 1
        }
        for(i in 1:5)
            mu[i] ~ dnorm(0,1)
    })

    ## Checking that no warning; if this were to warn, it would cause error
    nimbleOptions(verbose = FALSE)
    expect_silent(m <- nimbleModel(code, inits = list(k = rep(1,3)), constants = list(kk = 1:3)))
    expect_silent(m <- nimbleModel(code, inits = list(k = rep(1,3), kk = 1:3)))
    nimbleOptions(verbose = TRUE)
    
    code <- nimbleCode({
        for(i in 1:3) 
            y[i] ~ dnorm(mu[2, k[i]+j[i]],1)
        for(i in 1:3)
            for(ii in 1:4)
                mu[i, ii] ~ dnorm(0,1)
    })
    expect_message(m <- nimbleModel(code, inits = list(k = rep(1,3), j = rep(1,3))),
                   "Detected use of non-constant indexes")
    expect_message(m <- nimbleModel(code, inits = list(k = rep(1,3)), constants = list(j = rep(1,3))),
                   "Detected use of non-constant indexes")

    code <- nimbleCode({
        for(i in 1:3) 
            y[i] ~ dnorm(mu[k[i]+1, k[i]+j[i]],1)
        for(i in 1:3)
            for(ii in 1:4)
                mu[i, ii] ~ dnorm(0,1)
    })
    expect_message(m <- nimbleModel(code, inits = list(k = rep(1,3), j = rep(1,3))),
                   "Detected use of non-constant indexes")
    expect_message(m <- nimbleModel(code, inits = list(k = rep(1,3)), constants = list(j = rep(1,3))),
                   "Detected use of non-constant indexes")


    code <- nimbleCode({
        for(i in 1:3) 
            y[i,1:2] ~ dmnorm(mu[1:2, k[i]], prec[1:2,1:2])
        for(i in 1:3)
            mu[1:2, i] ~ dmnorm(z[1:2], prec[1:2,1:2])
    })
    expect_message(m <- nimbleModel(code, inits = list(k = rep(1,3))),
                   "Detected use of non-constant indexes")

    ## Checking that no warning; if this were to warn, it would cause error
    nimbleOptions(verbose = FALSE)
    expect_silent(m <- nimbleModel(code, constants = list(k = rep(1,3))))
    nimbleOptions(verbose = TRUE)
 

    ## To test for problem raised in issue #996
    code <- nimbleCode({
        for(i in 1:5)
            y[i] ~ dnorm(X[idx[i, 1],idx[i, 2]]*X[idx[i, 1], idx[i, 4]], 1)
        for(i in 1:6)
            for(j in 1:6)
                X[i,j] ~ dnorm(0,1)
    })
    
    expect_message(m <- nimbleModel(code, data = list(y=rnorm(5), idx = matrix(rep(1:6), 6, 6))),
                   "Detected use of non-constant indexes")
    expect_message(m <- nimbleModel(code, data = list(y=rnorm(5)), inits = list(idx = matrix(rep(1:6), 6, 6))),
                   "Detected use of non-constant indexes")

})

test_that("handling of contiguous blocks", {
    ## Before issue 1015 fixed, this would error out.
    code <- nimbleCode({
        R[1:n,1:n] <- exp(-dist[1:n, 1:n])
        for(i in 1:n) 
            y[i] ~ dnorm(0, sd = R[i,i])
    })
    n <- 4
    dd <- matrix(1, n, n); dd[1,4] <- dd[4,1] <- dd[2,3] <- dd[3,2] <- sqrt(2)
    diag(dd) <- 0
    model <- nimbleModel(code, constants = list(n = n), inits = list(dist = dd))

    indArr <- matrix(1, 3, 3)
    diag(indArr) <- c(2, 3, 1)
    out <- nimble:::makeVertexNamesFromIndexArray2(indArr, varName = 'y')
    expect_identical(out$names, c('y[1%.s%3, 1%.s%3]', 'y[1, 1]', 'y[2, 2]'))

    indArr <- matrix(1, 3, 3)
    diag(indArr) <- c(2, 3, 4)
    out <- nimble:::makeVertexNamesFromIndexArray2(indArr, varName = 'y')
    expect_identical(out$names, c('y[1%.s%3, 1%.s%3]', 'y[1, 1]', 'y[2, 2]', 'y[3, 3]'))

    indArr <- matrix(1, 4, 4)
    indArr[2:3, 1:3] <- 2
    out <- nimble:::makeVertexNamesFromIndexArray2(indArr, varName = 'y')
    expect_identical(out$names, c('y[1%.s%4, 1%.s%4]', 'y[2:3, 1:3]'))

    indArr <- matrix(1, 4, 4)
    indArr[2:3, c(1,3)] <- indArr[2:3, c(1,3)] <- 2
    out <- nimble:::makeVertexNamesFromIndexArray2(indArr, varName = 'y')
    expect_identical(out$names, c('y[1%.s%4, 1%.s%4]', 'y[2:3, 1%.s%3]'))

    indArr <- array(1, c(3, 3, 3))
    indArr[1:2, 1:2, 2] <- 2
    out <- nimble:::makeVertexNamesFromIndexArray2(indArr, varName = 'y')
    expect_identical(out$names, c('y[1%.s%3, 1%.s%3, 1%.s%3]', 'y[1:2, 1:2, 2]'))

    indArr <- array(1, c(3, 3, 3))
    indArr[1:2, 1:2, 1:2] <- 2
    out <- nimble:::makeVertexNamesFromIndexArray2(indArr, varName = 'y')
    expect_identical(out$names, c('y[1%.s%3, 1%.s%3, 1%.s%3]', 'y[1:2, 1:2, 1:2]'))

    indArr <- array(1, c(3, 3, 3))
    indArr[1, 1, 1] <- indArr[2, 2, 1] <- 2
    out <- nimble:::makeVertexNamesFromIndexArray2(indArr, varName = 'y')
    expect_identical(out$names, c('y[1%.s%3, 1%.s%3, 1%.s%3]', 'y[1%.s%2, 1%.s%2, 1]'))

    indArr <- array(1, c(3, 3, 3))
    indArr[1, 1, 1] <- indArr[2, 2, 2] <- 2
    out <- nimble:::makeVertexNamesFromIndexArray2(indArr, varName = 'y')
    expect_identical(out$names, c('y[1%.s%3, 1%.s%3, 1%.s%3]', 'y[1%.s%2, 1%.s%2, 1%.s%2]'))

    ## Another case that would formerly error out
    code <- nimbleCode({
        R[1:3,1:3] <- exp(-dist[1:3, 1:3])
	y[1] ~ dnorm(0, sd = R[1,1])
	y[2] ~ dnorm(0, sd = R[1,2])
	y[3] ~ dnorm(0, sd = R[2,1])
	y[4] ~ dnorm(0, sd = R[2,3])
	y[5] ~ dnorm(0, sd = R[3,2])
	y[6] ~ dnorm(0, sd = R[3,3])
    })
    n <- 4
    dd <- matrix(1, n, n); dd[1,4] <- dd[4,1] <- dd[2,3] <- dd[3,2] <- sqrt(2)
    diag(dd) <- 0
    model <- nimbleModel(code, inits = list(dist = dd[1:3,1:3]))

    
})

test_that("error produced when variable used in index", {
    code <- nimbleCode({
        for(i in 1:2)
            y[i] ~ dnorm(0,1)
    })
    m <- nimbleModel(code)
    idx <- 1
    ## Check introduced in v0.10.1 disabled because not robust; see NCT issue 293
    expect_failure(expect_error(m$expandNodeNames("y[idx]"), "parseEvalNumericMany: a variable was found"))
})

test_that("dmvt usage", {
    code <- nimbleCode({
    	 y1[1:n] ~ dmvt(z[1:n], pr[1:n, 1:n], 4)
	 y2[1:n] ~ dmvt(z[1:n], scale = pr[1:n, 1:n], df = 4)
    })
    n <- 3
    m <- nimbleModel(code, inits = list(pr = diag(rep(2,n))), constants = list(n = n))	 
    expect_identical(m$getParam('y1[1:3]', 'prec'), diag(rep(2, n)))
    expect_equal(m$getParam('y2[1:3]', 'prec'), diag(rep(0.5, n)))
})

test_that("bad size or dimension of initial values", {

    code <- nimbleCode({
        for(i in 1:3)
            z[i] ~ dnorm(0,1)
        a[1:3,1:2] <- b[1:3,1:2]
    })
    
    ## Bad dim for nimbleModel
    expect_error(m <- nimbleModel(code, inits = list(z = matrix(rnorm(9), 3))),
                 "inconsistent dimensionality")
    
    m <- nimbleModel(code)
    expect_message(m$setInits(list(z = matrix(rnorm(9), 3))), "Incorrect size or dimension")
    expect_output(cm <- compileNimble(m), "Incorrect number of dimensions")
    expect_identical(cm$z, rep(0, 3))    
    expect_message(cm$setInits(list(z = matrix(rnorm(9), 3))), "Incorrect size or dimension")
    
    expect_error(m <- nimbleModel(code, inits = list(b = rnorm(2))),
          "inconsistent dimensionality")         
    
    m <- nimbleModel(code)
    expect_message(m$setInits(list(b = rnorm(2))), "Incorrect size or dimension")
    expect_output(cm <- compileNimble(m), "R object of different size")
    expect_identical(cm$b, matrix(0, 3, 2))
    expect_message(cm$setInits(list(b = rnorm(2))), "Incorrect size or dimension")

    ## Too many values

    ## For better or worse, length of 5 gets baked in based on inits.
    ## TODO: do we want to reconsider whether this should error out?
    m <- nimbleModel(code, inits = list(z = rnorm(5)))
    cm <- compileNimble(m)
    expect_identical(m$z, cm$z)

    m <- nimbleModel(code)
    expect_message(m$setInits(list(z=rnorm(5))), "Incorrect size or dimension")
    expect_output(cm <- compileNimble(m), "R object of different size")
    expect_identical(cm$z, rep(0, 3))
    expect_message(cm$setInits(list(z=rnorm(5))), "Incorrect size or dimension")
    
    ### matrix

    ## For better or worse, 4x2 gets baked in based on inits.
    ## TODO: do we want to reconsider whether this should error out?
    m <- nimbleModel(code, inits = list(b = matrix(rnorm(8),4,2)))
    cm <- compileNimble(m)
    expect_identical(m$b, cm$b)

    m <- nimbleModel(code)
    expect_message(m$setInits(list(b = matrix(rnorm(8),4,2))), "Incorrect size or dimension")
    expect_output(cm <- compileNimble(m), "R object of different size")
    expect_identical(cm$b, matrix(0, 3, 2))
    expect_message(cm$setInits(list(b = matrix(rnorm(8),4,2))), "Incorrect size or dimension")
    

    ## Too few values

    expect_error(m <- nimbleModel(code, inits = list(z = rnorm(2))),
                 "dimensions specified are smaller")

    m <- nimbleModel(code)
    expect_message(m$setInits(list(z=rnorm(2))), "Incorrect size or dimension")
    expect_output(cm <- compileNimble(m), "R object of different size")
    expect_identical(cm$z, rep(0, 3))
    expect_message(cm$setInits(list(z=rnorm(2))), "Incorrect size or dimension")
    
    ### matrix

    expect_error(m <- nimbleModel(code, inits = list(b = matrix(rnorm(4),2,2))),
                 "dimensions specified are smaller")

    m <- nimbleModel(code)
    expect_message(m$setInits(list(b = matrix(rnorm(4),2,2))), "Incorrect size or dimension")
    expect_output(cm <- compileNimble(m), "R object of different size")
    expect_identical(cm$b, matrix(0, 3, 2))
    expect_message(cm$setInits(list(b = matrix(rnorm(4),2,2))), "Incorrect size or dimension")

})

test_that("Example of splitVertices bug from Issue 1268 works.", {
  code <- nimbleCode({
    for (i in 1:2) {
      for (j in 1:2){
        b[i,j] <- X[i, index[j]]
      }
    }
    for (i in 1:2) {
      y[i] <- sum(a[1:2, index[i]])
    }
  })
  expect_no_error(Rmodel <- nimbleModel(code, constants = list(index=c(1,1))))
})

test_that("Warning printed when indexing info in user environment.", {
    code <- nimbleCode({
        for(i in 1:N)
            y[i] ~ dnorm(0,1)
    })
    N <- 3
    expect_message(m <- nimbleModel(code, constants = list(foo=3)),
                   "Information has been found in the user's environment")
})


options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
