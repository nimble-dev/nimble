
# tests for BNP models using CRP prior for labels.

rm(list=ls())
library(nimble)
nimbleOptions(allowDynamicIndexing = TRUE)
source("./packages/nimble/R/BNP_distributions.R")
source("./packages/nimble/R/BNP_samplers.R")


#-- Model 1: 
#-- Model : y_i ~ N(theta_i, s2_i), s0 known
#          theta_i <- thetatilde[xi_i]
#          s2_i <- s2tilde[xi_i]
#          xi ~ CRP(alpha)
#          thetatilde_i ~ N(mu0, tau0)
#          s2tilde_i ~ IG(a0,b0)

#-- BUGS definition of the model
Code=nimbleCode(
  {
    for(i in 1:N){
      thetatilde[i] ~ dnorm(mu0, tau0) 
      s2tilde[i] ~ dinvgamma(a0, b0) 
    }
    xi[1:N] ~ dCRP(conc)
    
    for(i in 1:N){
      theta[i] <- thetatilde[xi[i]]
      s2[i] <- s2tilde[xi[i]]
      y[i] ~ dnorm(theta[i], var=s2[i])
    }
    conc<-1
    mu0<-0; tau0<-sqrt(20)
    a0<-1; b0 <- 0.1
  }
)

conc<-1; mu0<-0; tau0<-sqrt(20); a0<-1; b0<-0.1
Consts<-list(N=100)
set.seed(1)
Inits<-list(xi=sample(1:10, size=Consts$N, replace=TRUE), 
           thetatilde=rnorm(Consts$N, mu0, tau0),
           s2tilde=rinvgamma(Consts$N, a0, b0))
Data<-list(y=c(rnorm(Consts$N/2,-5, 2),rnorm(Consts$N/2,5, 2)))

#-- compiling the model:
model<-nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
cmodel<-compileNimble(model)

#-- MCMC configuration:
modelConf<-configureMCMC(model, print=TRUE)
#modelConf$removeSamplers(c("xi"), print="TRUE")
#modelConf$addSampler(c("xi"), type="MarginalizedG_xi")
modelConf$printSamplers(c("xi"))
modelMCMC=buildMCMC(modelConf)

#-- compiling the sampler
CmodelNewMCMC=compileNimble(modelMCMC, project=model,
                     resetFunctions=TRUE, showCompilerOutput = TRUE)

#-- MCMC samples
set.seed(1)
nsave=100
t1=proc.time()
CmodelNewMCMC$run(nsave)
proc.time()-t1

#-- results:
N=Consts$N
samples=as.matrix(CmodelNewMCMC$mvSamples)
xisamples=samples[, (2*N+1):(3*N)]#samples[, 1+(N+1):(2*N)]
s2samples=samples[, 1:N]
thetasamples=samples[, (N+1):(2*N)]
s2i=t(sapply(1:nsave, function(i)s2samples[i, xisamples[i,]]))
thetai=t(sapply(1:nsave, function(i)thetasamples[i, xisamples[i,]]))


apply(xisamples, 1, table)


#-- predictive:
sec=seq(-20,20,len=100)
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
      thetaNew[j,1]=rnorm(1,mu0, tau0)
      thetaNew[j,2]=rinvgamma(1, a0,b0)
      countnew=countnew+1
    }else{
      thetaNew[j,]=c(thetasamples[i, xisamples[i,index]],s2samples[i, xisamples[i,index]])
    }
  }
  fsam[i, ]=sapply(sec, function(x)sum(wj*sum(wj*dnorm(x,thetaNew[,1],sqrt(thetaNew[,2])))))
}
fhat=apply(fsam, 2, mean)

#-- Graphics
hist(Data$y, freq=FALSE, breaks=30)
points(sec,fhat, col="black", lwd=2, type="l")
