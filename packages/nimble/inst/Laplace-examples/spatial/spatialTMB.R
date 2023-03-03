library(TMB)
compile("spatial.cpp")
dyn.load(dynlib("spatial"))
## Read data
source("spatial_data.R")

## Euclidian distance matrix
dd <- as.matrix(dist(Z))
## Make AD functions
obj <- MakeADFun(data = list(n=100, y=y, X=X, dd=dd),
                 parameters=list(
                     b = c(1,1),
                     log_a = log(2), #1.428571,
                     log_sigma = -1 ,#-0.6931472,
                     u = rep(1,n)),
                 DLL = "spatial",
                 random = "u",
                 silent = TRUE
                 )
tmbtime <- system.time(tmbres <- nlminb(obj$par, obj$fn, obj$gr))
tmbtime2 <- system.time(tmbsumm <- summary(sdreport(obj)))
