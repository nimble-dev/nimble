model in "beetles-probit.bug" 
data in "beetles-data.R"
compile, nchains(2)
parameters in "beetles-inits.R"
initialize
update 1000
monitor alpha
monitor beta
monitor r.hat
monitor D
update 10000
coda *

