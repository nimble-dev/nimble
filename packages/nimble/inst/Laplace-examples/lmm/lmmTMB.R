require(TMB)
compile("lmm.cpp")
dyn.load(dynlib("lmm"))

## Test data
set.seed(123)
y <- rep(1900:2010,each=2)
year <- factor(y)
quarter <- factor(rep(1:4,length.out=length(year)))
period <- factor((y > mean(y))+1)
## Random year+quarter effect, fixed period effect:
B_nonsparse <- model.matrix(~year+quarter-1)
A_nonsparse <- model.matrix(~period-1)
B <- as(B_nonsparse,"dgTMatrix")
A <- as(A_nonsparse,"dgTMatrix")
u <- rnorm(ncol(B)) ## logsdu=0
beta <- rnorm(ncol(A))*100
eps <- rnorm(nrow(B),sd=1) ## logsd0=0
x <- as.numeric( A %*% beta + B %*% u + eps )

## Fit model
obj <- MakeADFun(data=list(x=x, B=B, A=A),
                 parameters=list(u=u*0, logsdu=1, logsd0=1, beta=beta*0),
                 random="u",
                 DLL="lmm",
                 ## inner.method = "BFGS", ## Skip using newton
                 silent=TRUE
                 )
tmbtime <- system.time(tmbres <- nlminb(obj$par, obj$fn, obj$gr))
tmbtime2 <- system.time(tmbsumm <- summary(sdreport(obj)))
