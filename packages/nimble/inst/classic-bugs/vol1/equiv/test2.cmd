model in "equiv.bug"
data in "equivmiss-data.R"
compile, nchains(2) 
inits in "equiv-init.R"
initialize
update 1000
monitor theta
monitor equivalence
monitor sigma
monitor Y[1,1]
monitor Y[3,2]
monitor Y[6,2]
update 10000 
coda *
