library(TMB)
compile("spatial.cpp")
dyn.load(dynlib("spatial"))
## Read data
source("spatial_data.R")

## Euclidian distance matrix
dd <- as.matrix(dist(Z))

obj <- MakeADFun(data = list(n=100, y=y, X=X, dd=dd),
                 parameters=list(
                     b = c(1,1),
                     ## a = 2, 
                     log_a = log(2), #1.428571,
                     log_sigma = -1 ,#-0.6931472,
                     u = rep(1,n)),
                 DLL = "spatial",
                 random = "u",
                 ## random.start = rep(1, 100),## Using this uniform start value makes TMB slower
                 ## random.start = expression(last.par[random]), ## Default
                 ## hessian = TRUE,
                 ## inner.method = "BFGS", ## Use optim (BFGS) for inner optimization
                 silent = TRUE
                 )
## We track the parameter values at which fn and gr are evaluated
tmbtime <- system.time(tmbres <- nlminb(obj$par, obj$fn, obj$gr))
tmbtime2 <- system.time(rep <- summary(sdreport(obj), "fixed"))
