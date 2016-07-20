model in "jaw-constant.bug"
data in "jaw-data.R" 
compile, nchains(2)
parameters in "jaw-inits.R"
initialize
load dic
update 1000
monitor beta0
monitor mu
monitor Sigma2
monitor RSS
monitor deviance
update 10000
coda *

