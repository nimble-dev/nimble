source("../../R/Rcheck.R")
load.module("glm")
m <- jags.model("epil3.bug", data=read.data("epil-data.R"), n.chains=2)
update(m, 1000)
x <- coda.samples(m, c("alpha0", "alpha.Base", "alpha.Trt", "alpha.BT",
                       "alpha.Age", "alpha.V4", "sigma.b1", "sigma.b"),
                  n.iter=10000, thin=10)
source("bench-test2.R")
check.fun()
