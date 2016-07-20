source("../../R/Rcheck.R")
d <- read.jagsdata("litters-data.R")
inits <- read.jagsdata("litters-init.R")
m <- jags.model("litters.bug", d, inits, n.chains=2)
update(m, 5000)
x <- coda.samples(m, c("mu","theta"), n.iter=50000, thin=50)
source("bench-test1.R")
check.fun()

