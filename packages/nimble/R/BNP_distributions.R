
#----------------------------------------------#
# Chinesse restaurant process distributions:
#----------------------------------------------#

# random number generator: rCRP
#   input:
#      n: sample size
#      conc: concentration parameter of the DP (related to sampling a new value)
#   output:
#      x: vector of size n
#
#   a new value is sampled with probability conc/(i-1+conc), an existing
#   value is sampled otherwise.

# density evaluation: dCRP
#   input: 
#      x: a vector
#      conc: concentration parameter of the DP
#      log: density returned in log scale or not
#   output:
#      x: a value

rCRP=nimbleFunction(
  run=function(n=integer(0), 
               conc=double(0, default=1))
  {
    returnType(double(1))
    x=numeric(n)# returned vector:
    x[1]=1
    if(n>1){
      for(i in 2:n){
        if(runif(1)<=conc/(conc+i-1)){# a new value
          x[i]=max(x[1:(i-1)])+1 
        }else{# an old value
          index=rcat(n=1, rep(1, i-1)) 
          x[i]=x[index]
        }
      }
    }
    return(x)
  }
)


dCRP=nimbleFunction(
  run=function(x=double(1), 
               conc=double(0, default=1), 
               log=integer(0, default=0))
  {
    returnType(double(0))
    n=length(x) 
    tmpden=numeric(n) 
    
    tmpden[1]=1
    if(n>1){
      for(i in 2:n){
        counts=0 # replaces sum(x[i]==x[1:(i-1)]) (works in nimble?)
        for(j in 1:(i-1)){
          if(x[i]==x[j]){
            counts=counts+1
          }
        }
        if(counts>0){
          tmpden[i]=1/(i-1+conc)
        }else{
          tmpden[i]=conc/(i-1+conc)
        }
      }
    }
    logProb <- sum(log(tmpden))
    if(log) return(logProb)
    else return(exp(logProb)) 
  }
)

