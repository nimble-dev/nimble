model in "asia2.bug"
data in asia2-data.R
compile, nchains(2) 
initialize
update 10000
monitor theta.b
monitor smoking[2]
monitor smoking[3]
update 10000 
coda *

