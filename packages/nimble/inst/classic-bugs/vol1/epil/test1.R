source("../../R/Rcheck.R")
load.module("glm")
m <- jags.model("epil2.bug", data=read.data("epil-data.R"), n.chains=2)
update(m, 1000)
x <- coda.samples(m, c("alpha0", "alpha.Base", "alpha.Trt", "alpha.BT",
                       "alpha.Age", "alpha.V4", "sigma.b1"), n.iter=5000,
                  thin=5)
source("bench-test1.R")
check.fun()

