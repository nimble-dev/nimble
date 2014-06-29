model in "jaw-quadratic.bug"
data in "jaw-data.R" 
compile, nchains(2)
parameters in "jaw-inits.R"
initialize
load dic
update 1000
monitor beta0
monitor beta1
monitor beta0.uncentred
monitor beta1.uncentred
monitor beta2
monitor Sigma2
monitor mu
monitor RSS
monitor deviance
update 10000
coda *

