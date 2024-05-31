source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)

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

test_that('noncentered sampler works', {

    ## error-trapping
    code <- nimbleCode({
        for(j in 1:g) {
            for(i in 1:n) {
                y[i,j] ~ dpois(exp(lambda[i,j]))
                lambda[i,j] <- b[j]
            }
            b[j] ~ dnorm(b0, sd = sigma)
        }
        d[1:2] ~ dmnorm(b0v[1:2], pr[1:2,1:2])
        for(i in 1:2)
            b0v[i] <- b0
        b0 ~ dflat()
        sigma ~ dunif(0, 50)
    })

    set.seed(5)
    n <- 30
    g <- 10
    b <- rnorm(g) + 1
    lambda <- matrix(rep(b, each = n), n, g)
    y <- matrix(rpois(n*g, exp(lambda)), n, g)

    m <- nimbleModel(code, data = list(y = y, d = rnorm(2)), constants = list(n = n, g = g),
                     inits = list(b0 = 0, sigma = 1, pr = diag(2)))
    conf <- configureMCMC(m, monitors = c('b0','b','sigma'))
    conf$addSampler('b0', 'noncentered')
    expect_error(mcmc <- buildMCMC(conf), "all dependent nodes must be univariate")

    code <- nimbleCode({
        for(j in 1:g) {
            for(i in 1:n) {
                y[i,j] ~ dpois(exp(lambda[i,j]))
                lambda[i,j] <- b[j]
            }
            b[j] ~ dexp(b0)
        }
        b0 ~ dflat()
    })

    set.seed(5)
    n <- 30
    g <- 10
    b <- rnorm(g) + 1
    lambda <- matrix(rep(b, each = n), n, g)
    y <- matrix(rpois(n*g, exp(lambda)), n, g)

    m <- nimbleModel(code, data = list(y = y), constants = list(n = n, g = g),
                     inits = list(b0 = 0, sigma = 1, pr = diag(2)))
    conf <- configureMCMC(m, monitors = c('b0','b'))
    conf$addSampler('b0', 'noncentered')
    expect_error(mcmc <- buildMCMC(conf), "the distribution `dexp` does not have")
    
    colSDs <- function(x) apply(x,2,sd)
    
    code <- nimbleCode({
        for(j in 1:g) {
            for(i in 1:n) {
                y[i,j] ~ dpois(exp(lambda[i,j]))
                lambda[i,j] <- b[j]
            }
            b[j] ~ dnorm(b0, sd = sigma)
        }
        b0 ~ dflat()
        sigma ~ dunif(0, 50)
    })

    set.seed(5)
    n <- 30
    g <- 10
    b <- rnorm(g) + 1
    lambda <- matrix(rep(b, each = n), n, g)
    y <- matrix(rpois(n*g, exp(lambda)), n, g)

    m <- nimbleModel(code, data = list(y = y), constants = list(n = n, g = g),
                     inits = list(b0 = 0, sigma = 1))
    conf <- configureMCMC(m, monitors = c('b0','b','sigma'))
    conf$addSampler('b0', 'noncentered')
    conf$addSampler('sigma', 'noncentered', control=list(logScale = TRUE, samplerParam = 'scale'))
    mcmc <- buildMCMC(conf)
    cm <- compileNimble(m)
    cmcmc <- compileNimble(mcmc, project = m)
    resultsNoncLog <- runMCMC(cmcmc, 26000, nburnin=1000)

    m <- nimbleModel(code, data = list(y = y), constants = list(n = n, g = g),
                     inits = list(b0 = 0, sigma = 1))
    conf <- configureMCMC(m, monitors = c('b0','b','sigma'))
    conf$addSampler('b0', 'noncentered', control=list(samplerType = 'slice'))
    conf$addSampler('sigma', 'noncentered', control=list(samplerType = 'slice', samplerParam = 'scale'))
    mcmc <- buildMCMC(conf)
    cm <- compileNimble(m)
    cmcmc <- compileNimble(mcmc, project = m)
    resultsNoncSlice <- runMCMC(cmcmc, 26000, nburnin=1000)

    m <- nimbleModel(code, data = list(y = y), constants = list(n = n, g = g),
                     inits = list(b0 = 0, sigma = 1))
    conf <- configureMCMC(m, monitors = c('b0','b','sigma'))
    mcmc <- buildMCMC(conf)
    cm <- compileNimble(m)
    cmcmc <- compileNimble(mcmc, project = m)
    resultsDefault <- runMCMC(cmcmc, 26000, nburnin = 1000)

    expect_equal(colMeans(resultsDefault), colMeans(resultsNoncLog), tolerance = .015)
    expect_equal(colSDs(resultsDefault), colSDs(resultsNoncLog), tolerance = .015)

    expect_equal(colMeans(resultsDefault), colMeans(resultsNoncSlice), tolerance = .015)
    expect_equal(colSDs(resultsDefault), colSDs(resultsNoncSlice), tolerance = .015)

    ## Check correctness of using 'location' vs 'scale' in nonsensical ways.
    m <- nimbleModel(code, data = list(y = y), constants = list(n = n, g = g),
                     inits = list(b0 = 0, sigma = 1))
    conf <- configureMCMC(m, monitors = c('b0','b','sigma'))
    conf$addSampler('b0', 'noncentered', control = list(samplerParam = 'scale'))
    conf$addSampler('sigma', 'noncentered', control=list(samplerParam = 'location')) 
    mcmc <- buildMCMC(conf)
    cm <- compileNimble(m)
    cmcmc <- compileNimble(mcmc, project = m)
    resultsNoncWeird <- runMCMC(cmcmc, 26000, nburnin=1000)

    m <- nimbleModel(code, data = list(y = y), constants = list(n = n, g = g),
                     inits = list(b0 = 0, sigma = 1))
    conf <- configureMCMC(m, monitors = c('b0','b','sigma'))
    conf$addSampler('b0', 'noncentered', control = list(samplerType = 'slice', samplerParam = 'scale'))
    conf$addSampler('sigma', 'noncentered', control=list(samplerType = 'slice', samplerParam = 'location')) 
    mcmc <- buildMCMC(conf)
    cm <- compileNimble(m)
    cmcmc <- compileNimble(mcmc, project = m)
    resultsNoncWeirdSlice <- runMCMC(cmcmc, 26000, nburnin=1000)

    expect_equal(colMeans(resultsDefault), colMeans(resultsNoncWeird), tolerance = .01)
    expect_equal(colSDs(resultsDefault), colSDs(resultsNoncWeird), tolerance = .01)

    expect_equal(colMeans(resultsDefault), colMeans(resultsNoncWeirdSlice), tolerance = .01)
    expect_equal(colSDs(resultsDefault), colSDs(resultsNoncWeirdSlice), tolerance = .015)

    
    ## nonstandard model (non-normal effects)
    code <- nimbleCode({
        for(j in 1:g) {
            for(i in 1:n) {
                y[i,j] ~ dpois(lambda[i,j])
                lambda[i,j] <- b[j]
            }
            b[j] ~ dgamma(mean = b0, sd = sigma)
        }
        b0 ~ dhalfflat()
        sigma ~ T(dnorm(0, sd = 10), 0, Inf) # dunif(0, 50)
    })

    set.seed(5)
    n <- 30
    g <- 10
    b <- rnorm(g) + 1
    lambda <- matrix(rep(b, each = n), n, g)
    y <- matrix(rpois(n*g, exp(lambda)), n, g)

    m <- nimbleModel(code, data = list(y = y), constants = list(n = n, g = g),
                     inits = list(b0 = 1, sigma = 1))
    conf <- configureMCMC(m, monitors = c('b0','b','sigma'))
    conf$addSampler('b0', 'noncentered')
    conf$addSampler('sigma', 'noncentered', control=list(logScale = TRUE, samplerParam = 'scale'))
    mcmc <- buildMCMC(conf)
    cm <- compileNimble(m)
    cmcmc <- compileNimble(mcmc, project = m)
    resultsNoncLog <- runMCMC(cmcmc, 101000, nburnin=1000)

    m <- nimbleModel(code, data = list(y = y), constants = list(n = n, g = g),
                     inits = list(b0 = 1, sigma = 1))
    conf <- configureMCMC(m, monitors = c('b0','b','sigma'))
    conf$addSampler('b0', 'noncentered', control = list(samplerType = 'slice'))
    conf$addSampler('sigma', 'noncentered', control = list(samplerType = 'slice', samplerParam = 'scale'))
    mcmc <- buildMCMC(conf)
    cm <- compileNimble(m)
    cmcmc <- compileNimble(mcmc, project = m)
    resultsNoncSlice <- runMCMC(cmcmc, 101000, nburnin=1000)

    m <- nimbleModel(code, data = list(y = y), constants = list(n = n, g = g),
                     inits = list(b0 = 1, sigma = 1))
    conf <- configureMCMC(m, 'b', monitors = c('b0','b','sigma'))
    conf$addSampler('b0','slice')  # better mixing than default
    conf$addSampler('sigma','slice')
    mcmc <- buildMCMC(conf)
    cm <- compileNimble(m)
    cmcmc <- compileNimble(mcmc, project = m)
    resultsSlice <- runMCMC(cmcmc, 101000, nburnin = 1000)

    m1 <- colMeans(resultsSlice)
    s1 <- colSDs(resultsSlice)
    m2 <- colMeans(resultsNoncLog)
    s2 <- colSDs(resultsNoncLog)
    expect_equal(m1[1:g], m2[1:g], tolerance = .01)
    expect_equal(m1[(g+1):(g+2)], m2[(g+1):(g+2)], tolerance = .05)
    expect_equal(s1[1:g], s2[1:g], tolerance = .01)
    expect_equal(s1[(g+1):(g+2)], s2[(g+1):(g+2)], tolerance = .1)

    m3 <- colMeans(resultsNoncSlice)
    s3 <- colSDs(resultsNoncSlice)
    expect_equal(m1[1:g], m3[1:g], tolerance = .01)
    expect_equal(m1[(g+1):(g+2)], m3[(g+1):(g+2)], tolerance = .05)
    expect_equal(s1[1:g], s3[1:g], tolerance = .01)
    expect_equal(s1[(g+1):(g+2)], s3[(g+1):(g+2)], tolerance = .1)
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)

