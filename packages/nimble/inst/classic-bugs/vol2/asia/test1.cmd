model in "asia.bug"
data in asia-data.R
compile, nchains(2) 
initialize
update 10000
monitor smoking
monitor tuberculosis 
monitor lung.cancer 
monitor bronchitis 
monitor either 
monitor xray 
update 10000 
coda *

