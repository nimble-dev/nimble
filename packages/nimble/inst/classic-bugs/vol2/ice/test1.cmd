model in "icear.bug"
data in "ice-data.R"
compile, nchains(2)
parameters in "ice-inits.R"
initialize
update 10000
monitor sigma, thin(100)
monitor logRR, thin(100)
update 100000
coda *

