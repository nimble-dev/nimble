model in "dugongs.bug"
data in "dugongs-data.R"
compile, nchains(2)
parameters in "dugongs-inits.R"
initialize
update 1000
monitor U3 
monitor alpha 
monitor beta 
monitor gamma 
monitor sigma 
update 10000 
coda *

