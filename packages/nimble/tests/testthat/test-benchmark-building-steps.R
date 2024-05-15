## This runs benchmarks of the steps of building a model, compiling
## the model, configuring the MCMC, building the MCMC, and compiling
## the MCMC.

source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

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
            buildMCMC(MCMCconfNoConj)
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

    ## 1000 is a good size for full benchmarking
    ## 10 is a good size for routine testing
    Benchmark1length <- 10
    ## following will be like 'theta->mu[1:10]->y[1:10]'
    caseNames[1] <- paste0('theta->mu[1:',
                           Benchmark1length,
                           ']->y[1:',
                           Benchmark1length,
                           ']')
    code1 <- nimbleCode({
        theta ~ dnorm(0,1)
        for(i in 1:Benchmark1length) mu[i] ~ dnorm(theta, sd = 1)
        for(i in 1:Benchmark1length) y[i] ~ dnorm(mu[i], sd = 1)
    })
    y1 <- rnorm(Benchmark1length, 0, 2)
    
    profile1 <- timeSteps(code = code1,
                          data = list(y = y1),
                          constants = list(
                              Benchmark1length = Benchmark1length)
                          )

    ## 100x20 is a good size for full benchmarking
    ## 5 x 2 is a good size for routing testing 
    Benchmark2dims <- c(5, 2)
    ## following will be like 'theta->mu[1:5]->y[1:5, 1:2]'
    caseNames[2] <- paste0('theta->mu[1:',
                           Benchmark2dims[1],
                           ']->y[1:',
                           Benchmark2dims[1],
                           '1:',
                           Benchmark2dims[2],
                           ']')
    code2 <- nimbleCode({
        theta ~ dnorm(0,1)
        for(i in 1:Benchmark2dims[1]) mu[i] ~ dnorm(theta, sd = 1)
        for(i in 1:Benchmark2dims[1])
            for(j in 1:Benchmark2dims[2]) y[i, j] ~ dnorm(mu[i], sd = 1)
    })
    y2 <- matrix(rnorm(prod(Benchmark2dims), 0, 2), nrow = Benchmark2dims[1])
    profile2 <- timeSteps(code = code2,
                          data = list(y = y2),
                          constants = list(
                              Benchmark2dims = Benchmark2dims)
                          )

    results <- rbind(profile1, profile2)
    rownames(results) <- caseNames

    print(results)
}
)

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
