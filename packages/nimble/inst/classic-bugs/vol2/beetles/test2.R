source("../../R/Rcheck.R")
data <- read.jagsdata("beetles-data.R")
inits <- read.jagsdata("beetles-inits.R")
m <- jags.model("beetles-probit.bug", data, inits, n.chains=2)
update(m, 1000)
x <- coda.samples(m, c("alpha","beta","r.hat","D"), n.iter=10000)
source("bench-test2.R")
check.fun()


