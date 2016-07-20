model in "rats.bug"
data in "rats-data.R"
compile, nchains(2)
inits in "rats-init.R"
initialize
update 1000 
monitor alpha0 
monitor beta.c 
update 10000
coda *
