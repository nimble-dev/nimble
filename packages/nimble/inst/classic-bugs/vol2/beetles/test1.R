source("../../R/Rcheck.R")
data <- read.jagsdata("beetles-data.R")
inits <- read.jagsdata("beetles-inits.R")
m <- jags.model("beetles-logit.bug", data, inits, n.chains=2)
update(m, 1000)
x <- coda.samples(m, c("alpha","beta","r.hat","D"), n.iter=10000)
source("bench-test1.R")
check.fun()

