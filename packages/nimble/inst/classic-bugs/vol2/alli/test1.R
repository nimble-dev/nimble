source("../../R/Rcheck.R")
d <- read.jagsdata("alli-data.R")
inits <- read.jagsdata("alli-inits.R")
m <- jags.model("alli.bug", d, inits, n.chains=2)
update(m, 1000)
x <- coda.samples(m, c("b[1,1]","b[1,2]","b[1,3]","b[1,4]","b[1,5]"), 
                  thin=10, n.iter=10000)
source("bench-test1.R")
check.fun()
