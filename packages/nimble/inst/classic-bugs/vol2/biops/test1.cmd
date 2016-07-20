model in "biops.bug"
data in "biops-data.R"
compile, nchains(2)
parameters in "biops-inits.R"
initialize
update 1000 
monitor p
monitor error[2,1]
monitor error[2,2]
monitor error[3,1]
monitor error[3,2]
monitor error[3,3]
monitor error[4,1]
monitor error[4,2]
monitor error[4,3]
monitor error[4,4]
update 10000
coda *

