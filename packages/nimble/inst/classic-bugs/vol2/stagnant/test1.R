source("../../R/Rcheck.R")
data <- read.jagsdata("stagnant-data.R")
inits <- read.jagsdata("stagnant-inits3.R")
m <- jags.model("stagnant2.bug", data, inits, n.chains=2)
update(m, 1000)
x <- coda.samples(m, c("alpha", "x.change"), thin=10, n.iter=100000)
source("bench-test1.R")
check.fun()

