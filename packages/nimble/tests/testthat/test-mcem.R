source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)
nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)

context("Testing of MCEM")

nimbleOptions(verbose = TRUE)
goldFileName <- 'mcemTestLog_Correct.Rout'
tempFileName <- 'mcemTestLog.Rout'
generatingGoldFile <- !is.null(nimbleOptions('generateGoldFileForMCEMtesting'))
outputFile <- if(generatingGoldFile) file.path(nimbleOptions('generateGoldFileForMCEMtesting'), goldFileName) else tempFileName

sink(outputFile)

test_that("Test that MCEM finds the MLE", {
    pumpCode <- nimbleCode({
        for (i in 1:N){
            theta[i] ~ dgamma(alpha,beta)
            lambda[i] <- theta[i]*t[i]
            x[i] ~ dpois(lambda[i])
        }
        alpha ~ dexp(1.0)
        beta ~ dgamma(0.1,1.0)
    })
    
    pumpConsts <- list(N = 10,
                       t = c(94.3, 15.7, 62.9, 126, 5.24,
                             31.4, 1.05, 1.05, 2.1, 10.5))
    
    pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))
    
    pumpInits <- list(alpha = 1, beta = 1,
                      theta = rep(0.1, pumpConsts$N))
    
    pump <- nimbleModel(code = pumpCode, name = 'pump',
                        constants = pumpConsts,
                        data = pumpData,
                        inits = pumpInits,
                        check = FALSE)
    
    compileNimble(pump)
    
    ## build an MCEM algorithm with Ascent-based convergence criterion
    pumpMCEM <- buildMCEM(model = pump,
                          latentNodes = 'theta', burnIn = 300,
                          mcmcControl = list(adaptInterval = 100),
                          boxConstraints = list( list( c('alpha', 'beta'),
                                                      limits = c(0, Inf) ) ),
                          C = 0.01, alpha = .01, beta = .01, gamma = .01, buffer = 1e-6)
    ## C changed from .001 to .01 to make the test run faster
    ## Correspondingly changed tolerance below from 0.01 to 0.04
    set.seed(0)
    out <- pumpMCEM$run(initM = 1000)
    names(out) <- NULL
    
    mle <- c(0.82, 1.26)
    
    expect_equal(out, mle, tolerance = 0.04)
})

## below example obtained from https://people.ucsc.edu/~abrsvn/bayes_winbugs_jags_4.r
## covariance matrix from buildMCEM compared to cov. matrix obtained from lme4 package

test_that("Test that asymptotic covariance is correct", {
    mixedEffCode <- nimbleCode({
        ## Priors
        mu.int~dnorm(0, 0.0001) # mean for random intercepts
        mu.slope~dnorm(0, 0.0001) # mean for random slopes
        sigma.int~dunif(0, 100) # SD of intercepts
        sigma.slope~dunif(0, 100) # SD of slopes
        rho~dunif(-1, 1) # correlation between intercepts and slopes
        Sigma.B[1, 1] <- pow(sigma.int, 2) # We start assembling the var-covar matrix for the random effects
        Sigma.B[2, 2] <- pow(sigma.slope, 2)
        Sigma.B[1, 2] <- rho*sigma.int*sigma.slope
        Sigma.B[2, 1] <- Sigma.B[1, 2]
        covariance <- Sigma.B[1, 2]
        Tau.B[1:2, 1:2] <- inverse(Sigma.B[1:2, 1:2])
        for (i in 1:ngroups) {
            B.hat[i, 1] <- mu.int
            B.hat[i, 2] <- mu.slope
            B[i, 1:2]~dmnorm(B.hat[i, 1:2], Tau.B[1:2, 1:2]) # the pairs of correlated random effects
            alpha[i] <- B[i, 1] # random intercept
            beta[i] <- B[i, 2] # random slope
        }
        sigma~dunif(0, 100) # Residual standard deviation
        tau <- 1/(sigma*sigma)
                                        # Likelihood
        for (i in 1:n) {
            mu[i] <- alpha[groups[i]]+beta[groups[i]]*x[i]
            y[i]~dnorm(mu[i], tau)
        }
    })

    mixedEffGroups <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4,
                        4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 
                        8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10)
    mixedEffData <- list(y = c(248.4, 254.8, 380.4, 332.7, 329.1, -31.7, 2.8, 
                               181.6, 228.6, 105.5, 226.3, 229.7, 208, 246.2, 
                               215.8, 113.4, 275.5, 194.2, 252.8, 207.4, 290.3, 
                               229, 220.9, 252.1, 196.1, 262.7, 244.3, 233.5,
                               289.3, 207.3, 68.7, 149.6, 157.3, 243.5, 174.2, 
                               133.5, 199.8, 284.8, 272, 240.8, 328, 309, 372.4,
                               32.6, 282.9, 167.5, 217.1, 175.5, 183.6, 145.1),
                         x =  c(-0.4, -0.5, 1.4, 0.8, 0.5, -1.9, -2, -0.5, -0.1,
                                -1.3, 1.5, 1.2, -0.6, 0.6, -0.7, -0.8, 1.1, -0.2,
                                0.4, 0.1, 1.4, -0.5, -0.9, -0.6, -0.3, 0.4, 0.4, 
                                0.3, 1.3, -0.5, -1.8, -0.3, -0.5, 0.6, -0.3, -1.3, 
                                1.3, 1.4, 1.2, 0, 1.4, 1.2, 1.3, -1.9, 0.1, 0.4, 0.5, -0.5, -0.9, -1.5))
    mixedEffConsts = list(ngroups=max(mixedEffGroups), n=length(mixedEffGroups), groups = mixedEffGroups)
    mixedEffMod <- nimbleModel(mixedEffCode, data = mixedEffData, constants = mixedEffConsts)
    compileNimble(mixedEffMod)
    mixedEffMCEM <- buildMCEM(model = mixedEffMod,
                              latentNodes = 'B', burnIn = 300,
                              mcmcControl = list(adaptInterval = 100),
                              boxConstraints = list(),
                              C = 0.05, alpha = .001, beta = .01, gamma = .01, buffer = 1e-6)
    set.seed(0)
    mixedEffMLEs <- c( mu.int = 221.79, mu.slope = 60.6,  sigma.int = 24.69, sigma.slope = 31.87, rho = 0.24,  sigma = 25.54) 
    mixedEffCov <- mixedEffMCEM$estimateCov(mixedEffMLEs)
    lme4Cov <- matrix(c(87.16, 22.2, 22.2, 135.3), nrow = 2)
    ## test that st.devs and sqrt of covs are close enough to truth.  we only examine cov of mu.int and mu.slope
    expect_equal(sqrt(lme4Cov[1,1]), sqrt(mixedEffCov[1,1]), tolerance =  2)
    expect_equal(sqrt(lme4Cov[2,2]), sqrt(mixedEffCov[2,2]), tolerance =  2)
    expect_equal(sign(lme4Cov[2,1]), sign(mixedEffCov[2,1]))
    expect_equal(sqrt(abs(lme4Cov[2,1])), sqrt(abs(mixedEffCov[2,1])), tolerance = 1)
})

sink(NULL)

if(!generatingGoldFile) {
    test_that("Log file matches gold file", {
        trialResults <- readLines(tempFileName)
        correctResults <- readLines(system.file(file.path('tests', 'testthat', goldFileName), package = 'nimble'))
        compareFilesByLine(trialResults, correctResults)
    })
}

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)
