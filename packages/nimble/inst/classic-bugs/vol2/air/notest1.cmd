load glm
model in "air.bug"
data in "air-data.R"
compile, nchains(2)
inits in "air-inits.R"
initialize
update 10000
monitor theta0, thin(10)
monitor theta[2], thin(10)
monitor X, thin(10)
update 50000
coda *

