/* blocker example with t-distribution for treatment effect */
model in blockert.bug
data in blocker-data.R
compile, nchains(2)
parameters in blockert-init.R
initialize
update 1000 
monitor d, thin(10)
monitor delta.new, thin(10) 
monitor sigma, thin(10) 
update 10000 
coda *
