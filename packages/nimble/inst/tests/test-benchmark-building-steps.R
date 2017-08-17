## This runs benchmarks of the steps of building a model, compiling
## the model, configuring the MCMC, building the MCMC, and compiling
## the MCMC.

source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context('Benchmarking model and MCMC building and compiling steps')
cat('\n')

timeSteps <- function(code, data = list(), constants = list(), inits = list()) {
    times <- numeric()
    times['Build model'] <- system.time(
        Rmodel <-
            nimbleModel(code = code,
                        data = data,
                        constants = constants,
                        inits = inits,
                        check = FALSE,
                        calculate = FALSE)
    )[3]
    times['Compile model'] <- system.time(
        Cmodel <-
            compileNimble(Rmodel)
    )[3]
    times['Configure MCMC'] <- system.time(
        MCMCconf <-
            configureMCMC(Rmodel)
    )[3]
    times['Configure MCMC (No Conj)'] <- system.time(
        MCMCconfNoConj <-
            configureMCMC(Rmodel, useConjugacy = FALSE)
    )[3]
    times['Build MCMC (No Conj)'] <- system.time(
        RMCMC <-
            buildMCMC(MCMCconfNoCon)
    )[3]
    times['Compile MCMC'] <- system.time(
        CMCMC <-
            compileNimble(RMCMC, project = Rmodel)
    )[3]
    times
}

test_that('Benchmarking model and MCMC building and compiling steps',
{
    caseNames <- character()

    caseNames[1] <- 'theta->mu[1:1000]->y[1:1000]'
    code1 <- nimbleCode({
        theta ~ dnorm(0,1)
        for(i in 1:1000) mu[i] ~ dnorm(theta, sd = 1)
        for(i in 1:1000) y[i] ~ dnorm(mu[i], sd = 1)
    })
    y1 <- rnorm(1000, 0, 2)
    
    profile1 <- timeSteps(code = code1, data = list(y = y1))

    caseNames[2] <- 'theta->mu[1:100]->y[1:100, 1:20]'
    code2 <- nimbleCode({
        theta ~ dnorm(0,1)
        for(i in 1:100) mu[i] ~ dnorm(theta, sd = 1)
        for(i in 1:100)
            for(j in 1:20) y[i, j] ~ dnorm(mu[i], sd = 1)
    })
    y2 <- matrix(rnorm(2000, 0, 2), nrow = 100)
    profile2 <- timeSteps(code = code2, data = list(y = y2))

    results <- rbind(profile1, profile2)
    rownames(results) <- caseNames

    print(results)
}
)
