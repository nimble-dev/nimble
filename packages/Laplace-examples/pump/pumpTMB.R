library(TMB)
compile("pumpTMB.cpp")
dyn.load(dynlib("pumpTMB"))

## Read data
N <- 10
t <- c(94.3, 15.7, 62.9, 126, 5.24, 31.4, 1.05, 1.05, 2.1, 10.5)
X <- c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)

data <- list(X = X, t = t)
parameters <- list(u = rep(log(0.1), N), alphat = log(0.1), betat = log(0.1))

## Fit model
obj <- MakeADFun(data, parameters, random="u", DLL="pumpTMB")
tmbtime <- system.time(tmbres <- nlminb(obj$par, obj$fn, obj$gr))
## Back transform model parameters
tmbres$par <- exp(tmbres$par) 
names(tmbres$par) <- c("alpha", "beta")
## Return standard errors
tmbsumm <- summary(sdreport(obj), select = "report")
