source("../../R/Rcheck.R")
data <- read.jagsdata("orange-data.R")
inits <- read.jagsdata("orange-inits.R")
m <- jags.model("otree.bug", data, inits, n.chains=2)
update(m, 1000)
x <- coda.samples(m, c("mu","sigma","sigmaC"), n.iter=10000)
source("bench-test1.R")
check.fun()

