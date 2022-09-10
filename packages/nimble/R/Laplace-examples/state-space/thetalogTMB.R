library(TMB)
compile("thetalog.cpp")
dyn.load(dynlib("thetalog"))

## Read data
y <- scan("thetalog.dat", skip=3, quiet=TRUE)
data <- list(y=y)

## Parameter initial guess
parameters <- list(
  u = data$y * 0,
  logr0 = 0,
  logpsi = 0,
  logK = 6,
  logQ = 0,
  logR = 0
)

## Fit model
obj <- MakeADFun(data, parameters, random="u", DLL="thetalog")
tmbtime <- system.time(tmbres <- nlminb(obj$par, obj$fn, obj$gr))
tmbsumm <- summary(sdreport(obj), "fixed")
