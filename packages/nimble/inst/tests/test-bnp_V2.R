rm(list=ls())
library(nimble)
nimbleOptions(allowDynamicIndexing = TRUE)



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

#--------------------------------------------------------------------
# sampler for a normal model with mean and variance randomly indexed:
MarginalizedG_xi<-nimbleFunction(
  
  contains=sampler_BASE,
  
  setup=function(model, mvSaved, target, control){
    ## note that even in inefficient case, we need to do individual dataNodes[i] <- model$getDependencies(targetElements[i], stochOnly = TRUE) because we are not guaranteed that xi[i] is the cluster membership for y[i]; it could be xi[i] is associated with y[n-i+1], e.g.
    
    calcNodes <- model$getDependencies(target)
    targetElements <- model$expandNodeNames(target, returnScalarComponents=TRUE)
    n <- length(targetElements)# number of observations
    
    dataNodes <- rep("", n)
    type <- 'indivCalcs'
    nInterm <- length(model$getDependencies(targetElements[1], determOnly = TRUE))
    #if(nInterm >= 1)
    intermNodes <- dataNodes
    #if(nInterm == 2)
    intermNodes2 <- dataNodes
    intermNodes3 <- dataNodes
    if(nInterm > 3) type <- "allCalcs"  ## give up and do the inefficient approach
    for(i in seq_len(n)) {
      stochDeps <- model$getDependencies(targetElements[i], stochOnly = TRUE, self=FALSE) #CHANGE
      detDeps <- model$getDependencies(targetElements[i], determOnly = TRUE)
      if(length(stochDeps) != 1) # CHANGE
        stop("Nimble cannot currently assign a sampler to a dCRP node unless each cluster indicator is associated with a single observation.")  ## reason for this is that we do getLogProb(dataNodes[i]), which assumes a single stochastic dependent
      if(length(detDeps) != nInterm) {
        type <- 'allCalcs'  # give up again; should only occur in strange situations
      } else {
        dataNodes[i] <- stochDeps[1] # CHANGE deps
        if(nInterm >= 1)
          intermNodes[i] <- detDeps[1]; intermNodes2[i] <- detDeps[1];
          intermNodes3[i] <- detDeps[1]
          if(nInterm >= 2)
            intermNodes2[i] <- detDeps[2]
          if(nInterm >= 3)
            intermNodes3[i] <- detDeps[3]
      }
    }
    
    curLogProb <- numeric(n) # stores the los probabilities of sampling existing or not indicators
    
  },
  
  
  run = function() {
    conc <- model$getParam(target, 'conc')
    #  -- udating xi:
    for(i in 1:n){ # updates each xi_i at the time , i=1,...,n
      xi_i <- model[[target]][i]
      xi <- model[[target]]
      cond <- sum(xi_i==xi) # if cond=1, xi_i is a singleton
      for(j in 1:n){ # calculate probability of sampling indexes 1,...,n   
        if(i==j){ # index i denotes a new indicator xi_i
          if(cond>1){ # a new parameter has to be created to calculate the prob
            newind <- 1
            mySum <- sum(xi == newind)
            while(mySum>0 & newind <= n) { # need to make sure don't go beyond length of vector
              newind <- newind+1
              mySum <- sum(xi == newind)
            }
            #while(sum(xi == newind)>0 & newind <= n) { # need to make sure don't go beyond length of vector
            #  newind <- newind+1
            #}
            model[[target]][i] <<- newind
            if(type == 'indivCalcs') {
              model$calculate(intermNodes[i])
              if(nInterm >= 2) model$calculate(intermNodes2[i])
              if(nInterm >= 3) model$calculate(intermNodes3[i])
              model$calculate(dataNodes[i])
            } else model$calculate(calcNodes) 
          }else{ # we keep the old parameter as the "new" one
            newind <- xi_i
            model[[target]][i] <<- newind
            if(type == 'indivCalcs') {
              model$calculate(intermNodes[i])
              if(nInterm >= 2) model$calculate(intermNodes2[i])
              if(nInterm >= 3) model$calculate(intermNodes3[i])
              model$calculate(dataNodes[i])
            } else model$calculate(calcNodes) 
          }
          curLogProb[j] <<- log(conc) + model$getLogProb(dataNodes[i]) #<<-
        }else{
          model[[target]][i] <<- model[[target]][j]
          if(type == 'indivCalcs') {
            model$calculate(intermNodes[i])
            if(nInterm >= 2) model$calculate(intermNodes2[i])
            if(nInterm >= 3) model$calculate(intermNodes3[i])
            model$calculate(dataNodes[i])
          } else model$calculate(calcNodes) 
          curLogProb[j] <<- model$getLogProb(dataNodes[i]) #<<-
        }  
        model[[target]][i] <<- xi_i
      } # 
      
      index <- rcat(n=1, exp(curLogProb-max(curLogProb)))#
      if(index==i){# creates a new component: one that is not used
        model[[target]][i] <<- newind
      }else{
        model[[target]][i] <<- model[[target]][index]
      } 
      #model$calculate(calcNodes)
    }
    model$calculate(calcNodes)
    
    
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {})
)


#--------------------------------------------------------------------
Code=nimbleCode(
  {
    for(i in 1:N){
      thetatilde[i] ~ dnorm(mean=mu0, var=tau0) 
      s2tilde[i] ~ dinvgamma(shape=a0, scale=b0) 
    }
    xi[1:N] ~ dCRP(conc)
    
    for(i in 1:N){
      theta[i] <- thetatilde[xi[i]]
      s2[i] <- s2tilde[xi[i]]
      y[i] ~ dnorm(theta[i], var=s2[i])
    }
    conc<-1; a0<-1; b0<-1; mu0<-0; tau0<-sqrt(30)
  }
)

conc<-1; a0<-1; b0<-1; mu0<-0; tau0<-sqrt(30)
Consts=list(N=100)
set.seed(1)
aux=sample(c(1,2), size=Consts$N, replace=TRUE)#c(rep(1,Consts$N/2), rep(2,Consts$N/2))#
Inits=list(xi=aux, thetatilde=rnorm(Consts$N, mu0, sqrt(tau0)),s2tilde=rinvgamma(Consts$N, shape=a0, scale=b0))#list(xi=aux, thetatilde=c(10,-10,rnorm(98, mu0,tau0)),s2tilde=rep(1,100))#

s0=2; s1=2
mu0=5; mu1=-5
Data=list(y=c(rnorm(Consts$N/2,mu0,s0), rnorm(Consts$N/2,mu1,s1)))

model=nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)#, dimensions=list(theta=Consts$N))
cmodel=compileNimble(model)

modelConf=configureMCMC(model)
modelConf$removeSamplers(c("xi"), print="TRUE")
modelConf$addSampler(c("xi"), type="MarginalizedG_xi")
modelConf$printSamplers(c("xi"))
modelMCMC2=buildMCMC(modelConf)

#-- compiling the new samplers:
CmodelNewMCMC=compileNimble(modelMCMC2, project=model,
                            resetFunctions=TRUE, showCompilerOutput = TRUE)

nsave=100
set.seed(2)
t1=proc.time()
CmodelNewMCMC$run(nsave)
proc.time()-t1

#-- the results:
N=Consts$N


samples=as.matrix(CmodelNewMCMC$mvSamples)
xisamples=samples[, (2*N+1):(3*N)]#samples[, 1+(N+1):(2*N)]
s2samples=samples[, 1:N]
thetasamples=samples[, (N+1):(2*N)]
s2i=t(sapply(1:nsave, function(i)s2samples[i, xisamples[i,]]))
thetai=t(sapply(1:nsave, function(i)thetasamples[i, xisamples[i,]]))

apply(xisamples, 1, table)

# posterior predictive:
sec=seq(-20,20,len=100)
#-- predictive 2:
Trunc=25
thetaNew=matrix(0, ncol=2, nrow=Trunc)
countnew=0
vj=c(rbeta(Trunc-1, 1, N+conc), 1)
wj=c(vj[1], vj[2:(Trunc-1)]*cumprod(1-vj[1:(Trunc-2)]), cumprod(1-vj[1:(Trunc-1)])[Trunc-1])
probs=c(rep(1/(N+conc), N), conc/(N+conc))
fsam=matrix(0,ncol=length(sec), nrow=nsave)
for(i in 1:nsave){
  for(j in 1:Trunc){
    index=sample(1:(N+1), size=1, prob=probs)
    if(index==(N+1)){
      thetaNew[j,1]=rnorm(1,mu0, tau0)#
      thetaNew[j,2]=rinvgamma(1, a0,b0)#rnorm(1,mu0, tau0)#
      countnew=countnew+1
    }else{
      thetaNew[j,]=c(thetasamples[i, xisamples[i,index]],s2samples[i, xisamples[i,index]])#thetasample[i, index]
    }
  }
  fsam[i, ]=sapply(sec, function(x)sum(wj*sum(wj*dnorm(x,thetaNew[,1],sqrt(thetaNew[,2])))))#sapply(sec, function(x)sum(wj*sum(wj*dweibull(x,thetaNew,1))))#sapply(sec, function(x)sum(wj*dpois(x, thetaNew)))#
}
fhat=apply(fsam, 2, mean)


#-- Graphics
hist(Data$y, freq=FALSE, breaks=30)
points(sec,fhat, col="black", lwd=2, type="l")





