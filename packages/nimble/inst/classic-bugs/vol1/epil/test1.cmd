model in "epil2.bug"
data in epil-data.R
load glm
compile, nchains(2)
inits in epil-inits.R
initialize
update 1000
monitor alpha0, thin(5)
monitor alpha.Base, thin(5)
monitor alpha.Trt, thin(5)
monitor alpha.BT, thin(5)
monitor alpha.Age, thin(5)
monitor alpha.V4, thin(5)
monitor sigma.b1, thin(5)
update 5000
coda *
