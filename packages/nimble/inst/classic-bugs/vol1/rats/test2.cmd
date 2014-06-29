/* RATS example with missing data */
model in "rats.bug"
data in "ratsmiss-data.R"
compile, nchains(2)
inits in "rats-init.R"
initialize
update 1000
monitor alpha0 
monitor beta.c 
monitor Y[26,2]
monitor Y[26,3]
monitor Y[26,4]
monitor Y[26,5]
update 10000 
coda *
