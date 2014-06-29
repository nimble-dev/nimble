model in epil3.bug
data in epil-data.R
load glm
compile, nchains(2)
inits in epil-inits.R
initialize
update 1000 
monitor alpha0, thin(10)
monitor alpha.Base, thin(10)
monitor alpha.Trt, thin(10)
monitor alpha.BT, thin(10)
monitor alpha.Age, thin(10)
monitor alpha.V4, thin(10)
monitor sigma.b1, thin(10)
monitor sigma.b, thin(10)
update 10000
coda *
