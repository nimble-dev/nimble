source("../../R/Rcheck.R")
data <- read.jagsdata("dyes-data.R")
m <- jags.model("dyes.bug", data, n.chains=2)
update(m, 5000)
x <- coda.samples(m, c("theta","sigma2.within","sigma2.between","f.between"),
                   n.iter=100000, thin=50)
source("bench-test1.R")
check.fun()
