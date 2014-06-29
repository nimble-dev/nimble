model in "kidney.bug"
data in "kidney-data.R"
compile, nchains(2)
inits in "kidney-init.R"
initialize
update 10000
monitor alpha, thin(50)
monitor beta.age, thin(50)
monitor beta.sex, thin(50)
monitor beta.disease[2], thin(50) 
monitor beta.disease[3], thin(50) 
monitor beta.disease[4], thin(50) 
monitor r, thin(50) 
monitor sigma, thin(50) 
update 50000 
coda *
