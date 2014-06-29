model in "salm.bug"
data in "salm-data.R"
load glm
compile, nchains(2)
inits in "salm-init.R"
initialize
update 1000 
monitor alpha
monitor beta
monitor gamma
monitor sigma
load dic
monitor deviance
update 5000
coda *
