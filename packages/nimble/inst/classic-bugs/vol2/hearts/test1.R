source("../../R/Rcheck.R")
data <- read.jagsdata("hearts-data.R")
inits <- read.jagsdata("hearts-inits.R")
m <- jags.model("hearts.bug", data, inits, n.chains=2)
update(m, 1000)
x <- coda.samples(m, c("theta", "beta","p"), thin=20, n.iter=10000)
source("bench-test1.R")
check.fun()


