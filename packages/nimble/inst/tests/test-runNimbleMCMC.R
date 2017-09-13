
if(FALSE) {
runExample <- function(pumpData, ## data
                       expected.accepted, ## expected output
                       pumpInits, ## initial values
                       pumpConsts, ## constants
                       pumpCode, ## nimble code
                       nthin=1,
                       nrun=1000) { ## samples to thin
    samples <- runNimbleMCMC(code=pumpCode,
                             inits=pumpInits,
                             constants=pumpConsts,
                             data=pumpData)

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

    pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))

    pumpConsts <- list(N = N,
                       t = c(94.3, 15.7, 62.9, 126, 5.24,
                             31.4, 1.05, 1.05, 2.1, 10.5)[1:N])

    pumpInits <- list(alpha = 1, beta = 1,
                      theta = rep(0.1, pumpConsts$N))

    return(list(pumpData=pumpData,
                pumpCode=pumpCode,
                pumpInit=pumpInits,
                pumpConsts=pumpConsts))

}


test.examples <- list(
    list(name = 'too little data',
         model.input =  makePumpModelInput(4),
         expected.accepted = FALSE),
    list(name = 'rejected',
         model.input =  makePumpModelInput(10),
         expected.accepted = FALSE),
    list(name = 'accepted',
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
}

