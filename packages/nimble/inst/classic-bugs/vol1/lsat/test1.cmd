model in lsat.bug
data in lsat-data.R
load glm
compile, nchains(2)
inits in lsat-init.R
initialize
update 1000 
monitor alpha 
monitor beta 
update 2000 
coda *
