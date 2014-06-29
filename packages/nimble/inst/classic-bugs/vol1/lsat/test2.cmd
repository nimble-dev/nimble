model in lsat2.bug
data in lsat-data.R
load glm
compile, nchains(2)
inits in lsat2-init.R
initialize
update 1000 
monitor delta, thin(10)
monitor eta, thin(10)
update 10000 
coda *
