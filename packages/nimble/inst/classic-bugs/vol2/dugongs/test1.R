source("../../R/Rcheck.R")
data <- read.jagsdata("dugongs-data.R")
inits <- read.jagsdata("dugongs-inits.R")
m <- jags.model("dugongs.bug", data, inits, n.chains=2)
update(m, 1000)
x <- coda.samples(m, c("U3","alpha","beta","gamma","sigma"), n.iter=10000)
source("bench-test1.R")
check.fun()


