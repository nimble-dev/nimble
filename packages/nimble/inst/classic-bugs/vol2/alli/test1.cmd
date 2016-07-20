model in "alli.bug"
data in "alli-data.R"
compile, nchains(2)
parameters in "alli-inits.R"
initialize
update 1000
monitor b[1,1], thin(10)
monitor b[1,2], thin(10)
monitor b[1,3], thin(10)
monitor b[1,4], thin(10)
monitor b[1,5], thin(10)
/* monitor G2, thin(10) */
update 10000
coda *

