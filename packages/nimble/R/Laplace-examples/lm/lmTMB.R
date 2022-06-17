library(TMB)
compile("lm.cpp")
dyn.load(dynlib("lm"))
set.seed(123)
data <- list(Y = rnorm(10) + 1:100, x=1:100)
parameters <- list(a=0, b=0, logSigma=0)
obj <- MakeADFun(data, parameters, DLL="lm")
obj$hessian <- TRUE
tmbtime <- system.time(tmbres <- do.call("optim", obj))