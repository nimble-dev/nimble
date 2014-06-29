source("../../R/Rcheck.R")
data <- read.jagsdata("jaw-data.R")
inits <- read.jagsdata("jaw-inits.R")
m <- jags.model("jaw-quadratic.bug", data, inits, n.chains=2)
update(m, 1000)
load.module("dic")
x <- coda.samples(m, c("beta0","beta1","beta0.uncentred","beta1.uncentred",
                       "beta2","Sigma2","mu","RSS","deviance"),
                  n.iter=10000)
source("bench-test3.R")
check.fun()
