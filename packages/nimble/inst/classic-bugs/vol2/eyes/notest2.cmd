model in "eyes.bug" 
data in "eyes-data.R"
load mix
compile, nchains(2)
parameters in "eyes-inits1.R"
initialize
update 5000 
monitor P 
monitor lambda 
monitor sigma 
monitor l0
update 5000 
coda *

