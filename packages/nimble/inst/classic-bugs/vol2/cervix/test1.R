load.module("glm")
m <- jags.model("cervix.bug", data=read.data("cervix-data.R"),
                inits=read.data("cervix-inits.R"), n.chain=2)
update(m, 1000)
x <- coda.samples(m, c("beta0C","beta","phi","gamma1","gamma2"),
                  n.iter=10000, thin=10)
source("bench-test1.R")
check.fun()


