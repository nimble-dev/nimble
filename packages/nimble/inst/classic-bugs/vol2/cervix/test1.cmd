load glm
model in "cervix.bug" 
data in "cervix-data.R"
compile, nchains(2) 
parameters in "cervix-inits.R"
initialize
update 1000
monitor beta0C, thin(10) 
monitor beta, thin(10) 
monitor phi, thin(10) 
monitor gamma1, thin(10) 
monitor gamma2, thin(10) 
update 10000
coda *

