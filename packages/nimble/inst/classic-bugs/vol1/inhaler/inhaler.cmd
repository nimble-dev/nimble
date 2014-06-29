model in inhaler.bug
data in inhaler-data.R
compile, nchains(2)
parameters in inhaler-inits.R
initialize
update 1000
monitor a, thin(10)
monitor beta, thin(10)
monitor kappa, thin(10)
monitor log.sigma, thin(10)
monitor pi, thin(10)
monitor sigma, thin(10)
update 10000
coda *
