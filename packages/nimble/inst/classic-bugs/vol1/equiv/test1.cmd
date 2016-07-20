model in equiv.bug
data in equiv-data.R
compile, nchains(2) 
inits in equiv-init.R
initialize
update 1000
monitor theta
monitor equivalence
monitor sigma
update 10000 
coda *
