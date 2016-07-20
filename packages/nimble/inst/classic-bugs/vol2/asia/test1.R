source("../../R/Rcheck.R")
data <- read.jagsdata("asia-data.R")
m <- jags.model("asia.bug", data, n.chains=2)
update(m, 10000)
x <- coda.samples(m, c("smoking","tuberculosis","lung.cancer","bronchitis",
                       "either","xray"), n.iter=10000)
source("bench-test1.R")
check.fun()

