source("../../R/Rcheck.R")
data <- read.jagsdata("biops-data.R")
inits <- read.jagsdata("biops-inits.R")
m <- jags.model("biops.bug", data, inits, n.chains=2)
update(m, 10000)
x <- coda.samples(m, c("p","error[2,1]", "error[2,2]", "error[3,1]",
                       "error[3,2]", "error[3,3]", "error[4,1]",
                       "error[4,2]", "error[4,3]", "error[4,4]"), n.iter=10000)
source("bench-test1.R")
check.fun()
