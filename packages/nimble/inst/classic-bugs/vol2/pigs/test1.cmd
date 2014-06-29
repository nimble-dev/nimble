model in "pigs.bug"
data in "pigs-data.R"
compile, nchains(2)
parameters in "pigs-inits.R"
initialize
update 1000
monitor p, thin(5)
monitor i[3], thin(5)
update 50000
coda *

