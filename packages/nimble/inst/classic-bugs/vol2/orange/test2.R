source("../../R/Rcheck.R")
data <- read.jagsdata("mvotree-data.R")
inits <- read.jagsdata("mvotree-inits.R")
m <- jags.model("mvotree.bug", data, inits, n.chains=2, n.adapt=5000)
update(m, 5000)
x <- coda.samples(m, c("mu","sigma","sigmaC"), n.iter=100000, thin=10)
source("bench-test2.R")
check.fun()


