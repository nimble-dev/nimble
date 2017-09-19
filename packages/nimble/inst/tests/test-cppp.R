rm(list=ls())
library(nimble)
library(parallel)
library(testthat)
setwd("~/Dropbox")
source('occupancy/analysis/cppp/src/calcCPPP.R')


runExample <- function(pumpData, ## data
                       expected.accepted, ## expected cppp value
                       pumpInits, ## initial values
                       pumpConsts, ## constants
                       pumpCode, ## nimble code
                       nthin=1,
                       nrun=1000) { ## samples to thin
    ## build model
    R.model <- nimbleModel(code=pumpCode,
                           constants=pumpConsts,
                           data=pumpData,
                           inits=pumpInits,
                           check=FALSE,
                           calculate=FALSE)
    message('R model created')

    ## configure and build mcmc
    mcmc.spec <- configureMCMC(R.model,
                               print=FALSE,
                               monitors = c("alpha", "beta", "theta"),
                               thin=nthin)
    mcmc <- buildMCMC(mcmc.spec)
    message('MCMC built')

    ## compile model in C++
    D.model <- compileNimble(R.model)
    D.mcmc <- compileNimble(mcmc, project = R.model)
    D.mcmc$run(nrun)
    output <- as.matrix(D.mcmc$mvSamples)
    message('NIMBLE model compiled')

    pumpDiscMeasure <- nimbleFunction(
        setup = function(model, discFunctionArgs){
            dataNames <- discFunctionArgs[[1]]
            lambdaNames <- 'lambda'
        },
        run = function(){
            dataVals <- values(model, dataNames)
            lambdaVals <- values(model, lambdaNames)
            disc <- 0
            for(i in 1:dim(dataVals)[1]){
                disc <- disc + ((dataVals[i] - lambdaVals[i])/sqrt(lambdaVals[i]))
            }
            returnType(double(0))
            return(disc)
        },
        contains = virtualDiscFunction
    )

    mcmcGenFunc <- function(model){
        mcmc.spec <- configureMCMC(model,
                                   print=FALSE,
                                   monitors = c("alpha", "beta", "theta"),
                                   thin=nthin)
        mcmc <- buildMCMC(mcmc.spec)
    }

    pumpCPPP <- runCPPP(R.model,
                             origMCMCOutput=output,
                             mcmcCreator = mcmcGenFunc,
                             dataNames = 'x',
                             nCPPPMCMCIter = 100,
                             nPPPCalcIters = 10,
                             nSimPPPVals = 10,
                             burnInProp = 0.1,
                             thin = nthin,
                             nBootstrapSDReps=10,
                             discFuncGenerator=pumpDiscMeasure,
                             discFunctionArgs = list('x'),
                             nCores = 1)
    if(pumpCPPP$CPPP <= 0.05){
        passCPPP <- FALSE
    } else {
        passCPPP = TRUE
    }
}


makePumpModelInput <- function(N){
    ## N must be between 2 and 10 to match the original example with
    ## real data
    pumpCode <- nimbleCode({
        for (i in 1:N){
            theta[i] ~ dgamma(alpha, beta)
            lambda[i] <- theta[i]*t[i]
            x[i] ~ dpois(lambda[i])
        }
        alpha ~ dexp(1.0)
        beta ~ dgamma(0.1,1.0)
    })

    pumpConsts <- list(N = N,
                       t = c(94.3, 15.7, 62.9, 126, 5.24,
                             31.4, 1.05, 1.05, 2.1, 10.5)[1:N])

    pumpInits <- list(alpha = 1, beta = 1,
                      theta = rep(0.1, pumpConsts$N))

    return(list(pumpCode=pumpCode,
                pumpInit=pumpInits,
                pumpConsts=pumpConsts))

}


test.examples <- list(
    list(name = 'too little data',
         data = list(x = c(5, 1, 5, 14)),
         model.input =  makePumpModelInput(4),
         expected.accepted = FALSE),
    list(name = 'rejected',
         data = list(x = rep(1, 10)),
         model.input =  makePumpModelInput(10),
         expected.accepted = FALSE),
    list(name = 'accepted',
         data = list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)),
         model.input =  makePumpModelInput(10),
         expected.accepted = TRUE))


for (example in test.examples) {
    test_that(example$name, {
        pumpData <- example$data
        runExample(pumpData,
                   expected.accepted = example$expected.accepted,
                   example$model.input$pumpInits,
                   example$model.input$pumpConsts,
                   example$model.input$pumpCode)
    })
}

