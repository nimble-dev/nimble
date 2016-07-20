source("../../R/Rcheck.R")
data <- read.jagsdata("eyes-data.R")
inits <- read.jagsdata("eyes-inits.R")
m <- jags.model("eyes2.bug", data, inits, n.chains=2)
update(m, 1000)
x <- coda.samples(m, c("P", "lambda","sigma"), thin=20, n.iter=40000)
source("bench-test1.R")
check.fun()

