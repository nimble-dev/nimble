# the schools example compiles slowly as it has 12k nodes
# this testing is meant to assess speed of MCMC

useInits <- TRUE

example <- "schools"; vol <- "vol2"

if(!"package:nimble" %in% search()) {
  dir <- file.path("..", "..", "packages", "nimble", "inst", "classic-bugs")
  dir <- file.path(dir, vol, example)
} else {
  dir <- file.path(system.file("classic-bugs", package = "nimble"), vol, example)
}

Rmodel <- readBUGSmodel(example, dir = dir, useInits = useInits)
Cmodel <- compileNimble(Rmodel)


mcmcspec <- MCMCspec(Rmodel, control = list(adaptInterval=100))
vars <- c("theta", "phi", "beta", "alpha", "gamma", "T")
mcmcspec$addMonitors(vars)

Rmcmc <- buildMCMC(mcmcspec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

Rniter <- 5
set.seed(0);
system.time(Rmcmc(Rniter))
Cniter <- 1000
set.seed(0)
system.time(Cmcmc(Cniter))

require(R2jags)
setwd(system.file("classic-bugs", "vol2", "schools"))
jagsTime = system.time({out1 <- jags(data = 'schools-data.R',
  parameters.to.save = vars, n.chains = 1,
  n.iter = 1000, n.burnin = 200, n.thin = 1, model.file = 'schools.bug', DIC = FALSE,
  jags.seed = 0)})

