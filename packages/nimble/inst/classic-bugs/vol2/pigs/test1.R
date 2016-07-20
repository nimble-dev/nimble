source("../../R/Rcheck.R")
data <- read.jagsdata("pigs-data.R")
inits <- read.jagsdata("pigs-inits.R")
m <- jags.model("pigs.bug", data, inits, n.chains=2)
update(m, 1000)
x <- coda.samples(m, c("p", "i[3]"), thin=5, n.iter=50000)
source("bench-test1.R")
check.fun()

