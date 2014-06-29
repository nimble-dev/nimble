source("../../R/Rcheck.R")
data <- read.jagsdata("stagnant-data.R")
inits <- read.jagsdata("stagnant-inits4.R")
m <- jags.model("stagnant2.bug", data, inits, n.chains=2)
update(m, 1000)
x <- coda.samples(m, c("alpha", "x.change"), thin=10, n.iter=100000)
source("bench-test2.R")
check.fun()

