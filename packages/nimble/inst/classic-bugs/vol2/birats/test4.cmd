model in "birats4.bug"
data in birats-data.R
compile, nchains(2)
inits in birats-inits.R
initialize
update 1000
monitor mu.beta, thin(10) 
monitor sigma.beta, thin(10)
monitor r, thin(10)
monitor alpha0, thin(10)
update 10000 
coda *

