RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

context("Testing that SMC functionality throws errors pointing to nimbleSMC.")

test_that("PF MCMCs fail with message", {
    timeModelCode <- nimbleCode({
        x[1] ~ dnorm(mu_0, 1)
        y[1] ~ dnorm(x[1], 1)
        for(i in 2:t){
            x[i] ~ dnorm(x[i-1] * a + b, 1)
            y[i] ~ dnorm(x[i] * c, 1)
        }
        
        a ~ dunif(0, 1)
        b ~ dnorm(0, 1)
        c ~ dnorm(1,1)
        mu_0 ~ dnorm(0, 1)
    })
    
                                        # simulate some data
    t <- 25; mu_0 <- 1
    x <- rnorm(1 ,mu_0, 1)
    y <- rnorm(1, x, 1)
    a <- 0.5; b <- 1; c <- 1
    for(i in 2:t){
        x[i] <- rnorm(1, x[i-1] * a + b, 1)
        y[i] <- rnorm(1, x[i] * c, 1)
    }
    
                                        # build the model
    rTimeModel <- nimbleModel(timeModelCode, constants = list(t = t), 
                              data <- list(y = y), check = FALSE )
    
                                        # Set parameter values and compile the model
    rTimeModel$a <- 0.5
    rTimeModel$b <- 1
    rTimeModel$c <- 1
    rTimeModel$mu_0 <- 1
    
    timeConf <- configureMCMC(rTimeModel)   # default MCMC configuration
    
                                        # Add random walk PMCMC sampler with particle number optimization.
    expect_error(
        timeConf$addSampler(target = c("a", "b", "c", "mu_0"), type = "RW_PF_block",
                            control = list(propCov= diag(4), adaptScaleOnly = FALSE,
                                           latents = "x", pfOptimizeNparticles = TRUE))
    )
    
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
