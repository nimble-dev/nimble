source(system.file(file.path('tests', 'testthat', 'AD_test_utils.R'), package = 'nimble'))
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)

nimbleOptions(useADcholAtomic = TRUE)
nimbleOptions(useADsolveAtomic  = TRUE)
nimbleOptions(useADmatMultAtomic = TRUE)
nimbleOptions(useADmatInverseAtomic  = TRUE)

RwarnLevel <- options('warn')$warn
options(warn = 1)

relTol <- eval(formals(test_ADModelCalculate)$relTol)
relTol[3] <- 1e-6
relTol[4] <- 1e-4

# Standalone test
test_that("basics of R calls to PDinverse_logdet and dmnormAD work", {
  set.seed(1)
  cov <- crossprod(matrix(rnorm(9), nrow=3))
  PDild <- PDinverse_logdet(cov)
  prec <- solve(cov)
  expect_equal(PDild[1:9][upper.tri(prec, diag=TRUE)], prec[upper.tri(prec, diag=TRUE)])
  expect_equal(PDild[10], determinant(cov, TRUE)$modulus[1])
  mu <- c(.2, .3, .1)
  x <- c(.1, .2, .3)
  chol_cov <- chol(cov)
  expect_equal(dmnorm_inv_ld(x, mu, mat = cov, inv_ld = PDild, prec_param=FALSE),
               dmnorm_chol(x, mu, chol_cov, prec_param=FALSE))
  expect_equal(dmnorm_inv_ld(x, mu, mat = cov, inv_ld = PDild, prec_param=FALSE, log=TRUE),
               dmnorm_chol(x, mu, chol_cov, prec_param=FALSE, log=TRUE))

  PDild2 <- PDinverse_logdet(prec)
  expect_equal(PDild2[1:9][upper.tri(cov, diag=TRUE)], cov[upper.tri(cov, diag=TRUE)])
  expect_equal(PDild2[10], determinant(prec, TRUE)$modulus[1])
  chol_prec <- chol(prec)
  expect_equal(dmnorm_inv_ld(x, mu, mat = prec, inv_ld = PDild2, prec_param=TRUE),
               dmnorm_chol(x, mu, chol_prec, prec_param=TRUE))
  expect_equal(dmnorm_inv_ld(x, mu, mat = prec, inv_ld = PDild2, prec_param=TRUE, log=TRUE),
               dmnorm_chol(x, mu, chol_prec, prec_param=TRUE, log=TRUE))
})

# Test in a model with dmnormAD

set.seed(1)
n <- 3
locs <- cbind(runif(n),runif(n))
dd <- matrix(nrow = n, ncol = n)
for(i in 1:n) for(j in 1:n) dd[i,j] <- sqrt(sum((locs[i,]-locs[j,])^2))
code <- nimbleCode({
  sigma ~ dunif(0, 20)
  rho ~ dunif(0,3)
  cov[1:3,1:3] <- sigma*sigma*exp(-dd[1:3,1:3]/rho)
  inv_ld[1:10] <- PDinverse_logdet(cov[1:3, 1:3])
  x[1:3] ~ dmnormAD(mean = mu[1:3], cov = cov[1:3, 1:3])
})
inits <- list(sigma = 0.8, rho = 1.2, mu = c(.2, .3, .1))
constants = list(dd=dd)
data = list(x = c(.1, .2, .3))
model <- nimbleModel(code, data = data, constants = constants, inits = inits, buildDerivs=TRUE)

relTolTmp <- relTol
relTolTmp[2] <- 1e-6
relTolTmp[3] <- 1e-2
relTolTmp[4] <- 1e-2
relTolTmp[5] <- 1e-13

test_ADModelCalculate(model, useParamTransform = TRUE, relTol = relTolTmp,
                      newConstantNodes = list(dd = dd*1.1, x=c(.15, .25, .35)),
                      newUpdateNodes = list(mu = c(.22, .32, .12)),
                      checkCompiledValuesIdentical = FALSE, check01vs012jacIdentical = FALSE,
                      checkDoubleUncHessian = TRUE,
                      useFasterRderivs = TRUE, verbose = FALSE, name = 'dmnormAD with inv_ld')

## Depending on the order in which the test above appears in the sequence of testing,
## we have seen the following crash:
# ============================================
# testing HMC/MAP-based scenario
# --------------------------------------------
# Testing initial wrt values with initial constantNodes
#  Using wrt:  sigma rho
# corrupted size vs. prev_size while consolidating
# Aborted

cov <- crossprod(matrix(rnorm(25), 5))
mu <- c(0.2,0.1,0.3,0.15,0.25)
x <- c(0.1,0.2,0.15,0.3,0.12)

PDinverse_logdet_test <- make_AD_test2(
  op = list(
    name = "PDinverse_logdet with prec_param FALSE: positive definite inverse and log determinant test",
    opParam = list(name = "PDinverse_logdet"),
    expr = quote({
      pdl <- PDinverse_logdet(mat)
      out <- pdl
    }),
    args = list(
      mat = quote(double(2))
    ),
    outputType = quote(double(1))
  ),
  argTypes = c(mat='double(2)'),
  wrt = c('mat'),
  inputs = list(record = list(mat = cov),
                test   = list(mat = cov+0.1))
)
PDinverse_logdet_test_out <- test_AD2(PDinverse_logdet_test)

dmnormAD_test1 <- make_AD_test2(
  op = list(
    name = "dmnormAD test with cov input and log fixed",
    opParam = list(name = "dmormAD test 1"),
    expr = quote({
      pld <- PDinverse_logdet(cov)
      logProb <- dmnorm_inv_ld(x, mu, mat=cov, inv_ld=pld, prec_param=0, log=TRUE)
      out <- logProb
    }),
    args = list(
      x = quote(double(1)),
      mu = quote(double(1)),
      cov = quote(double(2))
    ),
    outputType = quote(double(0))
  ),
  argTypes = c(x='double(1)', mu='double(1)', cov='double(2)'),
  wrt = c('x','mu','mat'),
  inputs = list(record = list(cov = cov, mu=mu, x=x),
                test   = list(cov = cov+0.1, mu=mu-0.07, x=x+0.03))
)
dmnormAD_test1_out <- test_AD2(dmnormAD_test1)

dmnormAD_test1b <- make_AD_test2(
    op = list(
        name = "dmnormAD test with cov input and log dynamics",
        opParam = list(name = "dmormAD test 1"),
        expr = quote({
            pld <- PDinverse_logdet(cov)
            logProb <- dmnorm_inv_ld(x, mu, mat=cov, inv_ld=pld, prec_param=0, log=log)
            out <- logProb
        }),
        args = list(
            x = quote(double(1)),
            mu = quote(double(1)),
            cov = quote(double(2)),
            log = quote(double(0))
        ),
        outputType = quote(double(0))
    ),
    argTypes = c(x='double(1)', mu='double(1)', cov='double(2)', log='double(0)'),
    wrt = c('x','mu','mat'),
    inputs = list(record = list(cov = cov, mu=mu, x=x, log=1),
                  test   = list(cov = cov+0.1, mu=mu-0.07, x=x+0.03, log=0))
)
dmnormAD_test1b_out <- test_AD2(dmnormAD_test1b)

prec <- solve(cov)
dmnormAD_test2 <- make_AD_test2(
  op = list(
    name = "dmnormAD test with prec input and log fixed",
    opParam = list(name = "dmormAD test 1"),
    expr = quote({
      pld <- PDinverse_logdet(prec)
      logProb <- dmnorm_inv_ld(x, mu, mat=prec, inv_ld=pld, prec_param=1, log=TRUE)
      out <- logProb
    }),
    args = list(
      x = quote(double(1)),
      mu = quote(double(1)),
      prec = quote(double(2))
    ),
    outputType = quote(double(0))
  ),
  argTypes = c(x='double(1)', mu='double(1)', prec='double(2)'),
  wrt = c('x','mu','mat'),
  inputs = list(record = list(prec = prec, mu=mu, x=x),
                test   = list(prec = prec+0.1, mu=mu-0.07, x=x+0.03))
)
dmnormAD_test2_out <- test_AD2(dmnormAD_test2)

dmnormAD_test2b <- make_AD_test2(
    op = list(
        name = "dmnormAD test with prec input and log dynamic",
        opParam = list(name = "dmormAD test 1"),
        expr = quote({
            pld <- PDinverse_logdet(prec)
            logProb <- dmnorm_inv_ld(x, mu, mat=prec, inv_ld=pld, prec_param=1, log=log)
            out <- logProb
        }),
        args = list(
            x = quote(double(1)),
            mu = quote(double(1)),
            prec = quote(double(2)),
            log = quote(double(0))
        ),
        outputType = quote(double(0))
    ),
    argTypes = c(x='double(1)', mu='double(1)', prec='double(2)', log='double(0)'),
    wrt = c('x','mu','mat'),
    inputs = list(record = list(prec = prec, mu=mu, x=x, log=0),
                  test   = list(prec = prec+0.1, mu=mu-0.07, x=x+0.03, log=1))
)
dmnormAD_test2b_out <- test_AD2(dmnormAD_test2b)

## Check that rmnorm_inv_ld works and draws the same numbers as from rmnorm_chol
test_that("rmnorm_inv_ld works and matches rmnorm_chol", {
  set.seed(1)
  cov <- crossprod(matrix(rnorm(25), 5))
  PDild_cov <- PDinverse_logdet(cov)
  chol_cov <- chol(cov)

  prec <- solve(cov)
  PDild_prec <- PDinverse_logdet(prec)
  chol_prec <- chol(prec)

  mu <- c(0.2,0.1,0.3,0.15,0.25)

  # PD_ild (PDinverse_logdet) will be ignored for method 2, rmnorm_chol
  nf <- nimbleFunction(
    run = function(mu=double(1), mat=double(2), PDild=double(1),
                   prec_param=double(0), method = integer(0)) {
      if(method == 1)
        ans <- rmnorm_inv_ld(1, mu, mat=mat, inv_ld = PDild, prec_param=prec_param)
      if(method == 2)
        ans <- rmnorm_chol(1, mu, cholesky = mat, prec_param=prec_param)
      return(ans)
      returnType(double(1))
    }
  )
  cnf <- compileNimble(nf)

  # uncompiled
  set.seed(1); x_PLD_0 <- nf(mu, cov, PDild_cov, 0, 1)
  set.seed(1); x_chol_0 <- nf(mu, chol_cov, c(0), 0, 2)
  set.seed(1); x_PLD_1 <- nf(mu, prec, PDild_prec, 1, 1)
  set.seed(1); x_chol_1 <- nf(mu, chol_prec, c(0), 1, 2)

  # compiled
  set.seed(1); Cx_PLD_0 <- cnf(mu, cov, PDild_cov, 0, 1)
  set.seed(1); Cx_chol_0 <- cnf(mu, chol_cov, c(0), 0, 2)
  set.seed(1); Cx_PLD_1 <- cnf(mu, prec, PDild_prec, 1, 1)
  set.seed(1); Cx_chol_1 <- cnf(mu, chol_prec, c(0), 1, 2)

  # look all together
  #rbind(x_PLD_0, x_chol_0, x_PLD_1, x_chol_1,
  #      Cx_PLD_0, Cx_chol_0, Cx_PLD_1, Cx_chol_1)

  expect_equal(x_PLD_0, x_chol_0)
  expect_equal(x_PLD_1, x_chol_1)
  expect_equal(Cx_PLD_0, Cx_chol_0)
  expect_equal(Cx_PLD_1, Cx_chol_1)
  expect_equal(x_PLD_0, Cx_PLD_0)
  expect_equal(x_PLD_1, Cx_PLD_1)

  # Check use of recycling rule for the mean
  mu <- c(0, 100)

  set.seed(1); x_PLD_0 <- nf(mu, cov, PDild_cov, 0, 1)
  set.seed(1); x_chol_0 <- nf(mu, chol_cov, c(0), 0, 2)
  set.seed(1); x_PLD_1 <- nf(mu, prec, PDild_prec, 1, 1)
  set.seed(1); x_chol_1 <- nf(mu, chol_prec, c(0), 1, 2)

  # compiled
  set.seed(1); Cx_PLD_0 <- cnf(mu, cov, PDild_cov, 0, 1)
  set.seed(1); Cx_chol_0 <- cnf(mu, chol_cov, c(0), 0, 2)
  set.seed(1); Cx_PLD_1 <- cnf(mu, prec, PDild_prec, 1, 1)
  set.seed(1); Cx_chol_1 <- cnf(mu, chol_prec, c(0), 1, 2)

  expect_equal(x_PLD_0, x_chol_0)
  expect_equal(x_PLD_1, x_chol_1)
  expect_equal(Cx_PLD_0, Cx_chol_0)
  expect_equal(Cx_PLD_1, Cx_chol_1)
  expect_equal(x_PLD_0, Cx_PLD_0)
  expect_equal(x_PLD_1, Cx_PLD_1)

  mu <- c(0, 10, 20, 30, 40, 50, 60) # too long, only first 5 should be used

  set.seed(1); x_PLD_0 <- nf(mu, cov, PDild_cov, 0, 1)
  set.seed(1); x_chol_0 <- nf(mu, chol_cov, c(0), 0, 2)
  set.seed(1); x_PLD_1 <- nf(mu, prec, PDild_prec, 1, 1)
  set.seed(1); x_chol_1 <- nf(mu, chol_prec, c(0), 1, 2)

  # compiled
  set.seed(1); Cx_PLD_0 <- cnf(mu, cov, PDild_cov, 0, 1)
  set.seed(1); Cx_chol_0 <- cnf(mu, chol_cov, c(0), 0, 2)
  set.seed(1); Cx_PLD_1 <- cnf(mu, prec, PDild_prec, 1, 1)
  set.seed(1); Cx_chol_1 <- cnf(mu, chol_prec, c(0), 1, 2)

  expect_equal(x_PLD_0, x_chol_0)
  expect_equal(x_PLD_1, x_chol_1)
  expect_equal(Cx_PLD_0, Cx_chol_0)
  expect_equal(Cx_PLD_1, Cx_chol_1)
  expect_equal(x_PLD_0, Cx_PLD_0)
  expect_equal(x_PLD_1, Cx_PLD_1)
})

test_that("non-assignment of dmnormAD if cholesky param used", {
    code <- nimbleCode({
        y[1:3] ~ dmnorm(mu[1:3], cholesky = L[1:3,1:3], prec_param = 0)
    })
    expect_message(m <- nimbleModel(code, inits=list(mu=rep(0,3),L=diag(3))),
                   "Detected use of `cholesky` parameterization")
})

test_that("lifting of PDinverse_logdet", {
    code <- nimbleCode({
        y[1:3] ~ dmnorm(mu[1:3], Q[2, 1:3, 5:7, 1])
    })
    m <- nimbleModel(code, data = list(y=rnorm(3)))
    which_PDinverse <- grep("PDinverse_logdet", m$getNodeNames())
    expect_identical(length(m[[m$getNodeNames()[which_PDinverse]]]), 10L)
})

# getParam

test_that('getParam', {

    code = nimbleCode({
        a[1:4] ~ dmnorm(mu[1:4],pr[1:4,1:4])
    })
    pr1 = diag(4)
    pr1[1,2]=pr1[2,1]=.3
    pr1[1,3]=pr1[3,1]=-.1
    pr1[2,3]=pr1[3,2] = 0.4
    pr2 <- pr1
    pr1[1,2]=pr1[2,1]=.5
    pr1[1,3]=pr1[3,1]=-.2
    pr1[2,3]=pr1[3,2] = 0.3

    m = nimbleModel(code, inits =list(mu=rep(1,4), pr = pr1))
    cm = compileNimble(m)

    cm$pr <- pr2
    cm$calculate(cm$getDependencies('pr'))

    expect_identical(pr1, m$getParam('a', 'prec'))
    expect_identical(pr2, cm$getParam('a', 'prec'))
    expect_equal(solve(pr1), m$getParam('a', 'cov'))
    expect_equal(solve(pr2), cm$getParam('a', 'cov'))

    code = nimbleCode({
        a[1:3] ~ dmnorm(mu[1:3],cov=pr[1:3,1:3])
    })
    pr1 = diag(3)
    pr1[1,2]=pr1[2,1]=.3
    pr2 <- pr1
    pr1[1,2]=pr1[2,1]=.5
    m = nimbleModel(code, inits =list(mu=rep(1,3), pr = pr1))
    cm = compileNimble(m)

    cm$pr <- pr2
    cm$calculate(cm$getDependencies('pr'))

    expect_identical(pr1, m$getParam('a', 'cov'))
    expect_identical(pr2, cm$getParam('a', 'cov'))
    expect_equal(solve(pr1), m$getParam('a', 'prec'))
    expect_equal(solve(pr2), cm$getParam('a', 'prec'))
})

# test-models

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

# Polya-gamma

test_that('polyagamma validity checks', {
    ## Various valid basic model structures
    code <- nimbleCode({
        for(i in 1:n) {
            p0[i] <- expit(b0+b1*x[i])
            y0[i] ~ dbern(p0[i])

            logit(p1[i]) <- a0 + a1*x[i]
            y1[i]~dbern(p1[i])

            logit(p2[i]) <- inprod(x2[i,1:p], b2[1:p])
            y2[i]~dbin(p2[i], size = m)

            logit(p3[i])  <- x2[i,1:p] %*% b3[1:p]
            y3[i]~dbin(p3[i], 1)
        }
        m ~ dpois(5)

        b0~dnorm(0, sd=10)
        b1~dnorm(0, sd=10)
        a0~dnorm(0, sd=10)
        a1~dnorm(0, sd=10)

        b2[1:p] ~ dmnorm(mu[1:p], Q[1:p,1:p])
        for(i in 1:p)
            b3[i] ~ dnorm(0, sd = 10)
    })

    n <- 10
    p <- 3
    constants <- list(n=n, p=p, x = runif(n), x2 = matrix(runif(n*p),n))
    ys <- rep(1,n)
    data <- list(y0 = ys, y1 = ys, y2 = ys, y3 = ys)

    m <- nimbleModel(code, constants = constants, data = data)
    conf <- configureMCMC(m, nodes = 'm')
    expect_silent(conf$addSampler(type='polyagamma', target=c('b0','b1')))
    expect_silent(conf$addSampler(type='polyagamma', target=c('a0','a1')))
    expect_silent(conf$addSampler(type='polyagamma', target=c('b2')))
    expect_silent(conf$addSampler(type='polyagamma', target=c('b3')))
    expect_silent(mcmc <- buildMCMC(conf))


    ## dnorm + dmnorm
    code <- nimbleCode({
        for(i in 1:n) {
            p0[i] <- expit(b[1]+u[k[i]]+b[2]*x[i])
            y0[i] ~ dbern(p0[i])
        }

        b[1:2] ~ dmnorm(z[1:2], pr[1:2,1:2])

        for(i in 1:p)
            u[i] ~ dnorm(0,1)
    })


    n <- 9
    p <- 3
    constants <- list(n=n, p=p, x = runif(n), k = rep(1:3, each = 3), z=rep(0,2), pr=diag(2))
    ys <- rep(1,n)
    data <- list(y0 = ys)

    m <- nimbleModel(code, constants = constants, data = data)
    conf <- configureMCMC(m, nodes = NULL)
    expect_silent(conf$addSampler(type='polyagamma', target=c('b','u')))
    expect_silent(mcmc <- buildMCMC(conf))
    mcmc$run(1)
    expect_equal(mcmc$samplerFunctions[[1]]$X,
                     cbind(rep(1,n), constants$x, c(rep(1,3),rep(0,6)),
                           c(rep(0,3),rep(1,3),rep(0,3)),
                           c(rep(0,6),rep(1,3))))
    lens <- c(2,1,1,1); names(lens) <- c('b[1:2]','u[1]','u[2]','u[3]')
    expect_equal(mcmc$samplerFunctions[[1]]$nodeLengths, lens)

})

test_that('Wishart-dmnorm conjugacy', {
   n <- 3
   R <- diag(rep(1,3))
   mu <- 1:3
   Y <- matrix(rnorm(9), 3)
   data <- list(Y = Y, mu = mu, R = R)
   code <- nimbleCode( {
       for(i in 1:3) {
           Y[i, 1:3] ~ dmnorm(mu[1:3], Omega[1:3,1:3]);
       }
       Omega[1:3,1:3] ~ dwish(R[1:3,1:3], 4);
   })
   m <- nimbleModel(code, data = data)
   conf <- configureMCMC(m)
   expect_identical(conf$getSamplers()[[1]]$name, "conjugate_dwish_dmnormAD_identity")

   code <- nimbleCode( {
       for(i in 1:3) {
           Y[i, 1:3] ~ dmnorm(mu[1:3], cov = Omega[1:3,1:3]);
       }
       Omega[1:3,1:3] ~ dwish(R[1:3,1:3], 4);
   })
   m <- nimbleModel(code, data = data)
   conf <- configureMCMC(m)
   expect_identical(conf$getSamplers()[[1]]$name, "RW_wishart")

   code <- nimbleCode( {
       for(i in 1:3) {
           Y[i, 1:3] ~ dmnorm(mu[1:3], Omega[1:3,1:3]);
       }
       Omega[1:3,1:3] ~ dinvwish(R[1:3,1:3], 4);
   })
   m <- nimbleModel(code, data = data)
   conf <- configureMCMC(m)
   expect_identical(conf$getSamplers()[[1]]$name, "RW_wishart")

   code <- nimbleCode( {
       for(i in 1:3) {
           Y[i, 1:3] ~ dmnorm(mu[1:3], cov = Omega[1:3,1:3]);
       }
       Omega[1:3,1:3] ~ dinvwish(R[1:3,1:3], 4);
   })
   m <- nimbleModel(code, data = data)
   conf <- configureMCMC(m)
   expect_identical(conf$getSamplers()[[1]]$name, "conjugate_dinvwish_dmnormAD_identity")
})

# Run various NIMBLE tests that use dmnorm with dmnormAD instead.

# From test-mcmc.R

## verbose: set to FALSE
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

## MCMC progress bar: set to FALSE
nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)

## MCMC orderSamplersPosteriorPredictiveLast - save current setting
nimbleReorderPPsamplersSetting <- getNimbleOption('MCMCorderPosteriorPredictiveSamplersLast')

## MCMC use usePosteriorPredictiveSampler - save current setting
nimbleUsePosteriorPredictiveSamplerSetting <- getNimbleOption('MCMCusePosteriorPredictiveSampler')

## MCMC calculation include predictive dependencies - save current setting
nimbleUsePredictiveDependenciesSetting <- nimbleOptions('MCMCusePredictiveDependenciesInCalculations')

## MCMC warn about unsampled nodes - save current setting
nimbleWarnUnsampledNodesSetting <- nimbleOptions('MCMCwarnUnsampledStochasticNodes')

test_that('elliptical slice sampler setup', {
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
              samplers = list(list(type = 'ess', target = 'x')), avoidNestedTest = TRUE)

    })

test_that('block sampler on MVN node setup', {
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
                  list(type = 'RW_block', target = 'x[1:3]')), avoidNestedTest = TRUE)
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
})

test_that('second block sampler on multivariate node', {
### DT's model
    mu <- c(1,2,3)
    corr <- matrix(c(1,.8,0.3,.8,1,0,0.3,0,1), nrow=3)
    varr <- c(1,2,3)
    Sig <- diag(sqrt(varr))
    Q <- Sig %*% corr %*% Sig
    P <- solve(Q)

    code <- nimbleCode({
        x[1:3] ~ dmnorm(mu[1:3], prec = P[1:3,1:3])
    })
    data = list(P = P, mu = mu)

    test_mcmc(model = code, name = 'second block sampler on multivariate node', data = data, seed = 0, numItsC = 100000,
              results = list(mean = list(x = mu),
                             var = list(x = varr)),
              resultsTolerance = list(mean = list(x = rep(.1,3)),
                                      var = list(x = c(.1,.1,.1))),
              samplers = list(
                  list(type = 'RW_block', target = 'x[1:3]')), avoidNestedTest = TRUE)
    })

test_that('MVN conjugate setup', {
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
                                      cov = list(mu = matrix(.01, 3, 3))), avoidNestedTest = TRUE)


### scalar RW updates in place of conjugate mv update

    test_mcmc(model = code, name = 'two-level multivariate normal with scalar updaters', data = data, seed = 0, numItsC = 100000,
              results = list(mean = list(mu = muMeanTrue),
                             cov = list(mu = solve(muQtrue))),
              resultsTolerance = list(mean = list(mu = rep(.03,3)),
                                      cov = list(mu = matrix(.03, 3, 3))),
              samplers = list(list(type = 'RW', target = 'mu[1]'),
                              list(type = 'RW', target = 'mu[2]'),
                              list(type = 'RW', target = 'mu[3]')),
              removeAllDefaultSamplers = TRUE, avoidNestedTest = TRUE)

})


test_that('another MVN conjugate sampler setup', {
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
    Rmcmc <- buildMCMC(spec)

    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

    set.seed(0)
    Rmcmc$run(10)
    Rsamples <- as.matrix(Rmcmc$mvSamples)
    set.seed(0)
    Cmcmc$run(10)
    Csamples <- as.matrix(Cmcmc$mvSamples)

    expect_equal(round(as.numeric(Rsamples), 8),
                 c(0.97473128, 0.50438666, 1.1251132, 0.83830666, 0.74077066, 0.92935482, 0.83758372, 0.98708273, 1.24199937, 0.67348127, -0.54387714, -0.60713969, -0.51392796, -0.3176801, -0.34416529, -0.08530564, -0.47160157, -0.21996584, -0.20504917, -0.77287122, 0.78462584, 0.46103509, 0.43862813, 0.49343096, 0.61020864, 0.55088287, 0.53887202, 0.49863894, 0.62691318, 0.80142839, 0.34941152, 0.06623608, 0.05624477, 0.21369178, 0.26585415, -0.1439989, -0.03133488, 0.3544062, -0.03518959, 0.27415746, 0.40977, 0.8351078, 0.25719293, 0.05663917, 0.30894028, 0.33113315, 0.47647909, 0.26143962, 0.07180759, 0.27255767),
                 info = 'R sample not correct compared to known result')

    dif <- as.numeric(Rsamples - Csamples)
    expect_true(max(abs(dif)) < 1E-15, info = 'R and C equiv')

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
    expect_true(max(abs(dif_mean)) < 0.001, info = 'posterior mean')

    dif_cov <- as.numeric(cov(Csamples) - post_cov)
    expect_true(max(abs(dif_cov)) < 0.001, info = 'posterior cov')
})


test_that('conjugate Wishart setup', {
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
                                      sd = list(Omega = matrix(0.06, M, M))), avoidNestedTest = TRUE)
                                        # issue with Chol in R MCMC - probably same issue as in jaw-linear

})

test_that('conjugate Wishart setup with scaling', {
    set.seed(0)

    trueCor <- matrix(c(1, .3, .7, .3, 1, -0.2, .7, -0.2, 1), 3)
    covs <- c(3, 2, .5)
    tau <- 4

    trueCov = diag(sqrt(covs)) %*% trueCor %*% diag(sqrt(covs))
    Omega = solve(trueCov) / tau

    n = 20
    R = diag(rep(1,3))
    mu = 1:3
    Y = mu + t(chol(trueCov)) %*% matrix(rnorm(3*n), ncol = n)
    M = 3
    data <- list(Y = t(Y), n = n, M = M, mu = mu, R = R)

    code <- nimbleCode( {
        for(i in 1:n) {
            Y[i, 1:M] ~ dmnorm(mu[1:M], tmp[1:M,1:M])
        }
        tmp[1:M,1:M] <- tau * Omega[1:M,1:M]
        Omega[1:M,1:M] ~ dwish(R[1:M,1:M], 4);
    })

    newDf = 4 + n
    newR = R + tcrossprod(Y - mu)*tau
    OmegaTrueMean = newDf * solve(newR)

    wishRV <- array(0, c(M, M, 10000))
    for(i in 1:10000) {
        z <- t(chol(solve(newR))) %*% matrix(rnorm(3*newDf), ncol = newDf)
        wishRV[ , , i] <- tcrossprod(z)
    }
    OmegaSimTrueSDs = apply(wishRV, c(1,2), sd)

    m <- nimbleModel(code, data = data[1],constants=data[2:5],inits = list(Omega = OmegaTrueMean, tau = tau))
    conf <- configureMCMC(m)
    expect_equal(conf$getSamplers()[[1]]$name, "conjugate_dwish_dmnormAD_multiplicativeScalar",
                 info = "conjugate dmnormAD-dwish with scaling not detected")

    test_mcmc(model = code, name = 'conjugate Wishart, scaled', data = data, seed = 0, numItsC = 1000, inits = list(Omega = OmegaTrueMean, tau = tau),
              results = list(mean = list(Omega = OmegaTrueMean ),
                             sd = list(Omega = OmegaSimTrueSDs)),
              resultsTolerance = list(mean = list(Omega = matrix(.05, M,M)),
                                      sd = list(Omega = matrix(0.06, M, M))), avoidNestedTest = TRUE)
                                        # issue with Chol in R MCMC - probably same issue as in jaw-linear

})

test_that('using RW_wishart sampler on non-conjugate Wishart node', {
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
            Y[i, 1:M] ~ dmnorm(mu[1:M], Omega[1:M,1:M])
        }
        Omega[1:M,1:M] ~ dwish(R[1:M,1:M], 4)
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
    allData <- data
    constants <- list(n = allData$n, M = allData$M, mu = allData$mu, R = allData$R)
    data <- list(Y = allData$Y)
    inits <- list(Omega = OmegaTrueMean)

    Rmodel <- nimbleModel(code, constants, data, inits)
    Rmodel$calculate()
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('Omega', 'RW_wishart')
    Rmcmc <- buildMCMC(conf)
    compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
    Cmodel <- compiledList$model; Cmcmc <- compiledList$mcmc
    set.seed(0)
    samples <- runMCMC(Cmcmc, 10000)
    d1 <- as.numeric(apply(samples, 2, mean)) - as.numeric(OmegaTrueMean)
    difference <- sum(round(d1,9) - c(0.024469145, -0.011872571, -0.045035297, -0.011872571, -0.003443918, 0.009363410, -0.045035297,  0.009363410,  0.049971420))
    expect_true(difference == 0)
})


test_that('using RW_wishart sampler on inverse-Wishart distribution', {
    code <- nimbleCode( {
        for(i in 1:n) {
            Y[i, 1:M] ~ dmnorm(mu[1:M], cov = C[1:M,1:M])
        }
        C[1:M,1:M] ~ dinvwish(R[1:M,1:M], 4)
    })

    set.seed(0)
    trueCor <- matrix(c(1, .3, .7, .3, 1, -0.2, .7, -0.2, 1), 3)
    covs <- c(3, 2, .5)
    trueCov = diag(sqrt(covs)) %*% trueCor %*% diag(sqrt(covs))
    n = 20
    R = diag(rep(1,3))
    mu = 1:3
    Y = mu + t(chol(trueCov)) %*% matrix(rnorm(3*n), ncol = n)
    M = 3
    data <- list(Y = t(Y), n = n, M = M, mu = mu, R = R)

    newDf = 4 + n
    newR = R + tcrossprod(Y- mu)
    CTrueMean <- newR / newDf
    CTrueMeanVec <- as.numeric(CTrueMean)
    allData <- data
    constants <- list(n = allData$n, M = allData$M, mu = allData$mu, R = allData$R)
    data <- list(Y = allData$Y)
    inits <- list(C = CTrueMean)
    niter <- 50000

    Rmodel <- nimbleModel(code, constants, data, inits)
    Rmodel$calculate()
    conf <- configureMCMC(Rmodel)
    Rmcmc <- buildMCMC(conf)
    compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
    Cmodel <- compiledList$model; Cmcmc <- compiledList$mcmc
    set.seed(0)
    samples <- runMCMC(Cmcmc, niter)
    conjMean <- as.numeric(apply(samples, 2, mean))
    conjSD <- as.numeric(apply(samples, 2, sd))

    Rmodel <- nimbleModel(code, constants, data, inits)
    Rmodel$calculate()
    conf <- configureMCMC(Rmodel, nodes = NULL)
    conf$addSampler('C', 'RW_wishart')
    Rmcmc <- buildMCMC(conf)
    compiledList <- compileNimble(list(model=Rmodel, mcmc=Rmcmc))
    Cmodel <- compiledList$model; Cmcmc <- compiledList$mcmc
    set.seed(0)
    samples <- runMCMC(Cmcmc, niter)
    RWMean <- as.numeric(apply(samples, 2, mean))
    RWSD <- as.numeric(apply(samples, 2, sd))

    expect_true(all(round(as.numeric(RWMean - conjMean), 9) == c(-0.001651758, -0.009675571, 0.004894809, -0.009675571, 0.015533882, -0.008256095, 0.004894809, -0.008256095, 0.002119615)))
    expect_true(all(round(as.numeric(RWSD - conjSD), 9) == c(0.022803503, -0.010107015, 0.012342044, -0.010107015, 0.006191412, -0.000091101, 0.012342044, -0.000091101, 0.001340032)))
})

test_that('detect conjugacy when scaling Wishart, inverse Wishart cases', {
    code <- nimbleCode({
        mycov[1:p, 1:p]  <- lambda * Sigma[1:p,1:p] / eta
        y[1:p] ~ dmnorm(z[1:p], cov = mycov[1:p, 1:p])
        Sigma[1:p,1:p] ~ dinvwish(S[1:p,1:p], nu)
    })
    m  <- nimbleModel(code, constants = list(p = 3))
    expect_identical(length(m$checkConjugacy('Sigma')), 1L, 'inverse-Wishart case')

    code <- nimbleCode({
        mycov[1:p, 1:p]  <- lambda * Sigma[1:p,1:p] / eta
        y[1:p] ~ dmnorm(z[1:p], prec = mycov[1:p, 1:p])
        Sigma[1:p,1:p] ~ dwish(S[1:p,1:p], nu)
    })
    m  <- nimbleModel(code, constants = list(p = 3))
    expect_identical(length(m$checkConjugacy('Sigma')), 1L, 'Wishart case')

    code <- nimbleCode({
        mycov[1:p, 1:p]  <- lambda * Sigma[1:p,1:p] / eta
        y[1:p] ~ dmnorm(z[1:p], cov = mycov[1:p, 1:p])
        Sigma[1:p,1:p] ~ dwish(S[1:p,1:p], nu)
    })
    m  <- nimbleModel(code, constants = list(p = 3))
    expect_identical(length(m$checkConjugacy('Sigma')), 0L, 'inverse-Wishart not-conj case')

    code <- nimbleCode({
        mycov[1:p, 1:p]  <- lambda * Sigma[1:p,1:p] / eta
        y[1:p] ~ dmnorm(z[1:p], prec = mycov[1:p, 1:p])
        Sigma[1:p,1:p] ~ dinvwish(S[1:p,1:p], nu)
    })
    m  <- nimbleModel(code, constants = list(p = 3))
    expect_identical(length(m$checkConjugacy('Sigma')), 0L, 'Wishart not-conj case')

    code <- nimbleCode({
        mycov[1:p, 1:p]  <- lambda[1:p,1:p] * Sigma[1:p,1:p]
        y[1:p] ~ dmnorm(z[1:p], cov = mycov[1:p, 1:p])
        Sigma[1:p,1:p] ~ dinvwish(S[1:p,1:p], nu)
    })
    m  <- nimbleModel(code, constants = list(p = 3))
    expect_identical(length(m$checkConjugacy('Sigma')), 0L, 'Wishart case')

    code <- nimbleCode({
        mycov[1:p, 1:p]  <- lambda[1:p,1:p] %*% Sigma[1:p,1:p]
        y[1:p] ~ dmnorm(z[1:p], cov = mycov[1:p, 1:p])
        Sigma[1:p,1:p] ~ dinvwish(S[1:p,1:p], nu)
    })
    m  <- nimbleModel(code, constants = list(p = 3))
    expect_identical(length(m$checkConjugacy('Sigma')), 0L, 'Wishart case')
})


test_that("realized conjugacy links are working", {

    ## dmnorm cases
    code <- nimbleCode({
        mu[1:3] ~ dmnorm(z[1:3], pr[1:3,1:3])
        for(i in 1:2)
            y1[i, 1:3] ~ dmnorm(mu[1:3], pr[1:3,1:3])
        mn2[1:3] <- b0[1:3] + mu[1:3]
        for(i in 1:2) {
            y2[i, 1:3] ~ dmnorm(mn2[1:3], pr[1:3,1:3])
        }
        mn3[1:3] <- A[1:3,1:3]%*%mu[1:3]
        for(i in 1:2) {
            y3[i, 1:3] ~ dmnorm(mn3[1:3], pr[1:3,1:3])
        }
        mn4[1:3] <- b0[1:3] + A[1:3,1:3]%*%mu[1:3]
        for(i in 1:2) {
            y4[i, 1:3] ~ dmnorm(mn4[1:3], pr[1:3,1:3])
        }
    })
    m <- nimbleModel(code, data = list (y1 = matrix(rnorm(6),2),
                               y2 = matrix(rnorm(6),2),
                               y3 = matrix(rnorm(6),2),
                               y4 = matrix(rnorm(6),2)),
                     inits = list(b0 = rnorm(3), A=matrix(1:9, 3), pr = diag(3)))
    conf <- configureMCMC(m)
    mcmc <- buildMCMC(conf)

    expect_identical(conf$getSamplers()[[1]]$name, "conjugate_dmnormAD_dmnormAD_additive_dmnormAD_identity_dmnormAD_linear_dmnormAD_multiplicative")
    expect_identical(mcmc$samplerFunctions[[1]]$N_dep_dmnormAD_identity, 2L)
    expect_identical(mcmc$samplerFunctions[[1]]$N_dep_dmnormAD_additive, 2L)
    expect_identical(mcmc$samplerFunctions[[1]]$N_dep_dmnormAD_multiplicative, 2L)
    expect_identical(mcmc$samplerFunctions[[1]]$N_dep_dmnormAD_linear, 2L)

    expect_identical(c('dep_dmnormAD_identity_coeff', 'dep_dmnormAD_additive_coeff') %in%
                ls(mcmc$samplerFunctions[[1]]), rep(FALSE, 2))
    expect_identical(c('dep_dmnormAD_multiplicative_coeff', 'dep_dmnormAD_linear_coeff') %in%
                ls(mcmc$samplerFunctions[[1]]), rep(TRUE, 2))
    expect_identical(c('dep_dmnormAD_identity_offset', 'dep_dmnormAD_multiplicative_offset') %in%
                ls(mcmc$samplerFunctions[[1]]), rep(FALSE, 2))
    expect_identical(c('dep_dmnormAD_additive_offset', 'dep_dmnormAD_linear_offset') %in%
                ls(mcmc$samplerFunctions[[1]]), rep(TRUE, 2))

    expect_identical(mcmc$samplerFunctions[[1]]$dep_dmnormAD_identity_nodeNames, c('y1[1, 1:3]', 'y1[2, 1:3]'))
    expect_identical(mcmc$samplerFunctions[[1]]$dep_dmnormAD_additive_nodeNames, c('y2[1, 1:3]', 'y2[2, 1:3]'))
    expect_identical(mcmc$samplerFunctions[[1]]$dep_dmnormAD_multiplicative_nodeNames, c('y3[1, 1:3]', 'y3[2, 1:3]'))
    expect_identical(mcmc$samplerFunctions[[1]]$dep_dmnormAD_linear_nodeNames, c('y4[1, 1:3]', 'y4[2, 1:3]'))

    ## dwish case
    code <- nimbleCode({
        for(i in 1:2) {
            y1[i, 1:3] ~ dmnorm(mu[1:3], pr[1:3,1:3])
        }
        pr2[1:3,1:3] <- d*pr[1:3,1:3]
        for(i in 1:2) {
            y2[i, 1:3] ~ dmnorm(mu[1:3], pr2[1:3,1:3])
        }
        pr[1:3,1:3] ~ dwish(R[1:3,1:3], 8)
    })
    m <- nimbleModel(code, data = list (y1 = matrix(rnorm(6),2),
                               y2 = matrix(rnorm(6),2)),
                     inits = list(pr = diag(3), R = diag(3)))
    conf <- configureMCMC(m)
    mcmc <- buildMCMC(conf)

    expect_identical(conf$getSamplers()[[1]]$name, "conjugate_dwish_dmnormAD_identity_dmnormAD_multiplicativeScalar")
    expect_identical(mcmc$samplerFunctions[[1]]$N_dep_dmnormAD_identity, 2L)
    expect_identical(mcmc$samplerFunctions[[1]]$N_dep_dmnormAD_multiplicativeScalar, 2L)

    expect_identical('dep_dmnormAD_identity_coeff' %in%
                ls(mcmc$samplerFunctions[[1]]), FALSE)
    expect_identical('dep_dmnormAD_multiplicativeScalar_coeff' %in%
                ls(mcmc$samplerFunctions[[1]]), TRUE)
    expect_identical(c('dep_dmnormAD_identity_offset', 'dep_dmnormAD_multiplicativeScalar_offset') %in%
                ls(mcmc$samplerFunctions[[1]]), rep(FALSE, 2))

    expect_identical(mcmc$samplerFunctions[[1]]$dep_dmnormAD_identity_nodeNames, c('y1[1, 1:3]', 'y1[2, 1:3]'))
    expect_identical(mcmc$samplerFunctions[[1]]$dep_dmnormAD_multiplicativeScalar_nodeNames, c('y2[1, 1:3]', 'y2[2, 1:3]'))
})


options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)
nimbleOptions(buildModelDerivs = BMDopt)
