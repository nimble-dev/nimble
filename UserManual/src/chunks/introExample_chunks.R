## ---- inputPump

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

## ---- explorePump

pump <- nimbleModel(code = pumpCode, name = "pump", constants = pumpConsts,
                    data = pumpData, inits = pumpInits)

pump$getNodeNames()
pump$x
pump$logProb_x
pump$alpha
pump$theta
pump$lambda

## ---- plotPump

pump$plotGraph()

## ---- manipPump

# Show all dependencies of alpha and beta terminating in stochastic nodes
pump$getDependencies(c("alpha", "beta"))
# Now show only the deterministic dependencies
pump$getDependencies(c("alpha", "beta"), determOnly = TRUE)
# Check that the lifted node was initialized. 
pump[["lifted_d1_over_beta"]] # It was.
# Now let's simulate new theta values
set.seed(1) # This makes the simulations here reproducible
pump$simulate("theta")
pump$theta   # the new theta values
# lambda and logProb_x haven't been re-calculated yet
pump$lambda # these are the same values as above
pump$logProb_x
pump$getLogProb("x") # The sum of logProb_x
pump$calculate(pump$getDependencies(c("theta")))
pump$lambda  # Now they have.
pump$logProb_x

## ---- compilePump

Cpump <- compileNimble(pump)
Cpump$theta

## ---- nimbleMCMCpump

mcmc.out <- nimbleMCMC(code = pumpCode, constants = pumpConsts,
                       data = pumpData, inits = pumpInits,
                       nchains = 2, niter = 10000,
                       summary = TRUE, WAIC = TRUE,
                       monitors = c('alpha','beta','theta'))

names(mcmc.out)
mcmc.out$summary
mcmc.out$WAIC  

## ---- mcmcPump
pumpConf <- configureMCMC(pump, print = TRUE)
pumpConf$addMonitors(c("alpha", "beta", "theta"))

pumpMCMC <- buildMCMC(pumpConf)
CpumpMCMC <- compileNimble(pumpMCMC, project = pump)

niter <- 1000
set.seed(1)
samples <- runMCMC(CpumpMCMC, niter = niter)

par(mfrow = c(1, 4), mai = c(.6, .4, .1, .2))
plot(samples[ , "alpha"], type = "l", xlab = "iteration",
     ylab = expression(alpha))
plot(samples[ , "beta"], type = "l", xlab = "iteration",
     ylab = expression(beta))
plot(samples[ , "alpha"], samples[ , "beta"], xlab = expression(alpha),
     ylab = expression(beta))
plot(samples[ , "theta[1]"], type = "l", xlab = "iteration",
     ylab = expression(theta[1]))

acf(samples[, "alpha"]) # plot autocorrelation of alpha sample
acf(samples[, "beta"])  # plot autocorrelation of beta  sample

## ---- mcmcPump2
pumpConf$addSampler(target = c("alpha", "beta"), type = "RW_block",
                    control = list(adaptInterval = 100))
                                     
pumpMCMC2 <- buildMCMC(pumpConf)

# need to reset the nimbleFunctions in order to add the new MCMC
CpumpNewMCMC <- compileNimble(pumpMCMC2, project  = pump,
                              resetFunctions = TRUE)

set.seed(1)
CpumpNewMCMC$run(niter)
samplesNew <- as.matrix(CpumpNewMCMC$mvSamples)

par(mfrow = c(1, 4), mai = c(.6, .4, .1, .2))
plot(samplesNew[ , "alpha"], type = "l", xlab = "iteration",
     ylab = expression(alpha))
plot(samplesNew[ , "beta"], type = "l", xlab = "iteration",
     ylab = expression(beta))
plot(samplesNew[ , "alpha"], samplesNew[ , "beta"], xlab = expression(alpha),
     ylab = expression(beta))
plot(samplesNew[ , "theta[1]"], type = "l", xlab = "iteration",
     ylab = expression(theta[1]))

acf(samplesNew[, "alpha"]) # plot autocorrelation of alpha sample
acf(samplesNew[, "beta"])  # plot autocorrelation of beta  sample

## ---- mcemPump
pump2 <- pump$newModel()

box = list( list(c("alpha","beta"), c(0, Inf)))

pumpMCEM <- buildMCEM(model = pump2, latentNodes = "theta[1:10]",
                      boxConstraints = box)
pumpMLE <- pumpMCEM$run()

pumpMLE

## ---- dont-run-mcemPump

pumpMLE <- c(0.8221657, 1.2589865)

pumpMLE

## ---- nfPump

simNodesMany <- nimbleFunction(
    setup = function(model, nodes) {
        mv <- modelValues(model)
        deps <- model$getDependencies(nodes)
        allNodes <- model$getNodeNames()
    },
    run = function(n = integer()) {
        resize(mv, n)
        for(i in 1:n) {
            model$simulate(nodes)
            model$calculate(deps)
            copy(from = model, nodes = allNodes,
                 to = mv, rowTo = i, logProb = TRUE)
        }
    })

simNodesTheta1to5 <- simNodesMany(pump, "theta[1:5]")
simNodesTheta6to10 <- simNodesMany(pump, "theta[6:10]")

## ---- runPumpSimsR
set.seed(1)  # make the calculation repeatable
pump$alpha <- pumpMLE[1]
pump$beta <- pumpMLE[2]
# make sure to update deterministic dependencies of the altered nodes
pump$calculate(pump$getDependencies(c("alpha","beta"), determOnly = TRUE))
saveTheta <- pump$theta
simNodesTheta1to5$run(10)
simNodesTheta1to5$mv[["theta"]][1:2]
simNodesTheta1to5$mv[["logProb_x"]][1:2]

## ---- runPumpSimsC
CsimNodesTheta1to5 <- compileNimble(simNodesTheta1to5,
                                    project  = pump, resetFunctions = TRUE)
Cpump$alpha <- pumpMLE[1]
Cpump$beta <- pumpMLE[2]
Cpump$calculate(Cpump$getDependencies(c("alpha","beta"), determOnly = TRUE))
Cpump$theta <- saveTheta

set.seed(1)
CsimNodesTheta1to5$run(10)
CsimNodesTheta1to5$mv[["theta"]][1:2]
CsimNodesTheta1to5$mv[["logProb_x"]][1:2]
