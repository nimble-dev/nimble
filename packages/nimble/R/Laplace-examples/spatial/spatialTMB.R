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
                     a = 2, #1.428571,
                     log_sigma = -1 ,#-0.6931472,
                     u = rep(1,n)),
                 DLL = "spatial",
                 random = "u",
                 hessian = TRUE,
                 inner.method = "BFGS",
                 silent = TRUE
                 )

tmbtime <- system.time(
  tmbres <- nlminb(obj$par, obj$fn, obj$gr,
                   lower=c(-100.0, -100.0, 0.01, -3.0),
                   upper=c( 100.0,  100.0, 3.00,  3.0))
  )


## Modify the source code of the h function, a part of MakeADFun in TMB
## to check the Hessian matrix calculation in TMB
# library(Matrix)
# source("../h.R")
# obj$env$h <- h
# suppressMessages(attach(obj$env))
# obj$fn(opt$par)
#

# rep <- sdreport(obj)
# rep
