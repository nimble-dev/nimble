source("../../R/Rcheck.R")
data <- read.jagsdata("jaw-data.R")
inits <- read.jagsdata("jaw-inits.R")
m <- jags.model("jaw-constant.bug", data, inits, n.chains=2)
update(m, 1000)
load.module("dic")
x <- coda.samples(m, c("beta0", "mu", "Sigma2", "RSS", "deviance"),
                  n.iter=10000)
source("bench-test1.R")
check.fun()

