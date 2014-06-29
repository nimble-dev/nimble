source("../../R/Rcheck.R")
data <- read.jagsdata("birats-data.R")
inits <- read.jagsdata("birats-inits.R")
m <- jags.model("birats3.bug", data, inits, n.chains=2)
update(m, 1000)
x <- coda.samples(m, c("mu.beta","sigma.beta","alpha0"), thin=10,
                  n.iter=10000)
source("bench-test3.R")
check.fun()

