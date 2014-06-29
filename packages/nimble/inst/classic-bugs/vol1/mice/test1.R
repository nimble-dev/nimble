source("../../R/Rcheck.R")
d <- read.jagsdata("mice-data.R")
inits <- read.jagsdata("mice-init.R")
m <- jags.model("mice.bug", d, inits, n.chains=2)
update(m, 10000)
x <- coda.samples(m, c("veh.control","test.sub","pos.control","r","median"),
                  thin=50, n.iter=500000)
source("bench-test1.R")
check.fun()
