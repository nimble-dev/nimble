model in blocker.bug
data in blocker-data.R
load glm
compile, nchains(2)
parameters in blocker-init.R
initialize
update 3000
monitor d, thin(10)
monitor delta.new, thin(10)
monitor sigma, thin(10)
update 30000 
coda *
