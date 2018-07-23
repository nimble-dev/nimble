## @knitr usingModelsExample

mc <- nimbleCode({
    a ~ dnorm(0, 0.001)
    for(i in 1:5) {
        y[i] ~ dnorm(a, sd = 0.1)
        for(j in 1:3)
            z[i,j] ~ dnorm(y[i], tau)
    }
    tau ~ dunif(0, 20)
    y.squared[1:5] <- y[1:5]^2
})

model <- nimbleModel(mc, data = list(z = matrix(rnorm(15), nrow = 5)))

## @knitr understandingNimbleConcepts

{
    a ~ dnorm(0, 0.001)
    for(i in 1:5) {
        y[i] ~ dnorm(a, 0.1)
        for(j in 1:3)
            z[i,j] ~ dnorm(y[i], sd = 0.1)
    }
    y.squared[1:5] <- y[1:5]^2
}


## @knitr usingModelVars

model$a <- 5
model$a
model[["a"]]
model$y[2:4] <- rnorm(3)
model$y
model[["y"]][c(1, 5)] <- rnorm(2)
model$y
model$z[1,]

## @knitr getVarAndNodeNames
model$getVarNames()

model$getNodeNames()

## @knitr usingModelLogProbs

model$logProb_y
model$calculate("y")
model$logProb_y

## @knitr usingNodeNames

model[["y[2]"]]
model[["y[2]"]] <- -5
model$y
model[["z[2, 3]"]]
model[["z[2:4, 1:2]"]][1, 2]
model$z[2, 2]

## @knitr multivariateExpandNodeNames

multiVarCode <- nimbleCode({
    X[1, 1:5] ~ dmnorm(mu[], cov[,])
    X[6:10, 3] ~ dmnorm(mu[], cov[,])
})

multiVarModel <- nimbleModel(multiVarCode, dimensions =
                   list(mu = 5, cov = c(5,5)), calculate = FALSE)

multiVarModel$expandNodeNames("X[1,1:5]")

## @knitr calcSimGLPdemos

mc <- nimbleCode({
    a ~ dnorm(0, 0.001)
    for(i in 1:5) {
        y[i] ~ dnorm(a, 0.1)
        for(j in 1:3)
            z[i,j] ~ dnorm(y[i], sd = 0.1)
    }
    y.squared[1:5] <- y[1:5]^2
})

model <- nimbleModel(mc, data = list(z = matrix(rnorm(15), nrow = 5)))

model$a <- 1
model$y
model$simulate("y[1:3]")
# simulate(model, "y[1:3]")
model$y
model$simulate("y")
model$y
model$z
model$simulate(c("y[1:3]", "z[1:5, 1:3]"))
model$y
model$z
model$simulate(c("z[1:5, 1:3]"), includeData = TRUE)
model$z

## @knitr calcSimGLPdirect

# y2lp <- model$nodes[["y[2]"]]$calculate()
# y2lp
# model$nodes[["y[2]"]]$getLogProb()

## @knitr reinitPumpModel

pumpCode <- nimbleCode({ 
  for (i in 1:N){
      theta[i] ~ dgamma(alpha,beta);
      lambda[i] <- theta[i]*t[i];
      x[i] ~ dpois(lambda[i])
  }
  alpha ~ dexp(1.0);
  beta ~ dgamma(0.1,1.0);
})

pumpConsts <- list(N = 10,
               t = c(94.3, 15.7, 62.9, 126, 5.24,
                 31.4, 1.05, 1.05, 2.1, 10.5))

pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))

pumpInits <- list(alpha = 1, beta = 1,
              theta = rep(0.1, pumpConsts$N))

pump <- nimbleModel(code = pumpCode, name = "pump", constants = pumpConsts,
                    data = pumpData, inits = pumpInits)

## @knitr getVarAndNodeNamesPump

pump$getVarNames()

pump$getNodeNames()
pump$getNodeNames(determOnly = TRUE)
pump$getNodeNames(stochOnly = TRUE)
pump$getNodeNames(dataOnly = TRUE)


## @knitr expandNodeNames

multiVarCode2 <- nimbleCode({
    X[1, 1:5] ~ dmnorm(mu[], cov[,])
    X[6:10, 3] ~ dmnorm(mu[], cov[,])
    for(i in 1:4) 
        Y[i] ~ dnorm(mn, 1)
})

multiVarModel2 <- nimbleModel(multiVarCode2,
                              dimensions = list(mu = 5, cov = c(5,5)),
                             calculate = FALSE)


multiVarModel2$expandNodeNames("Y")

multiVarModel2$expandNodeNames(c("X", "Y"), returnScalarComponents = TRUE)


## @knitr getDependencies

pump$getDependencies("alpha")
pump$getDependencies(c("alpha", "beta"))
pump$getDependencies("theta[1:3]", self = FALSE)
pump$getDependencies("theta[1:3]", stochOnly = TRUE, self = FALSE)
# get all dependencies, not just the direct descendants
pump$getDependencies("alpha", downstream = TRUE)
pump$getDependencies("alpha", downstream = TRUE, dataOnly = TRUE)
