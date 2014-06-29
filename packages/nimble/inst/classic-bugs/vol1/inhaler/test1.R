source("../../R/Rcheck.R")
d <- read.jagsdata("inhaler-data.R")
inits <- read.jagsdata("inhaler-inits.R")
m <- jags.model("inhaler.bug", d, inits, n.chains=2)
update(m, 1000)
x <- coda.samples(m, c("a","beta","kappa","log.sigma","pi","sigma"),
                  n.iter=1000)
source("bench-test1.R")
check.fun()
