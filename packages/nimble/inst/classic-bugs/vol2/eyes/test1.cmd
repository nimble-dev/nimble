/* We put a restriction on the mixture model, so that each group has at
   least one observation.  This is done by *pretending* that we have
   some extra data: a censored indicator of the number of observations
   in each group that tells us that there is at least one member.

   This may or may not be cheating.
*/
model in "eyes2.bug" 
data in "eyes-data.R"
compile, nchains(2)
parameters in "eyes-inits.R"
initialize
update 1000 
monitor P, thin(20)
monitor lambda, thin(20) 
monitor sigma, thin(20) 
update 40000 
coda *

