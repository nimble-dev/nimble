source("../../R/Rcheck.R")
data <- read.jagsdata("asia2-data.R")
m <- jags.model("asia2.bug", data, n.chains=2)
update(m, 10000)
x <- coda.samples(m, c("theta.b","smoking[2]","smoking[3]"), n.iter=10000)
source("bench-test2.R")
check.fun()

