model in "pump.bug"
data in "pump-data.R"
compile, nchains(2)
inits in "pump-init.R"
initialize
update 1000 
monitor theta
monitor alpha
monitor beta
update 5000
coda *
