
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#-- TESTS FOR BNP MODELS




rm(list=ls())
library(nimble)
#nimbleOptions(allowDynamicIndexing = TRUE)
#source("./packages/nimble/R/BNP_distributions.R")
#source("./packages/nimble/R/BNP_samplers.R")

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#-- tests for BNP models using CRP prior for labels.

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
    for(i in 1:N3){
      thetatilde[i] ~ dnorm(mean=mu0, var=tau20) 
      s2tilde[i] ~ dinvgamma(shape=a0, scale=b0) 
    }
    xi[1:N2] ~ dCRP(conc)
    
    for(i in 1:N){
      theta[i] <- thetatilde[xi[i]]
      s2[i] <- s2tilde[xi[i]]
      y[i] ~ dnorm(theta[i], var=s2[i])
    }
    conc<-1; a0<-1; b0<-0.5; mu0<-0; tau20<-40
  }
)

conc<-1; a0<-1; b0<-0.5; mu0<-0; tau20<-40
Consts=list(N=50, N2=50, N3=50)
set.seed(1)
aux=sample(1:10, size=Consts$N2, replace=TRUE)
Inits=list(xi=aux, thetatilde=rnorm(Consts$N3, mu0, sqrt(tau20)),s2tilde=rinvgamma(Consts$N3, shape=a0, scale=b0))#list(xi=aux, thetatilde=c(10,-10,rnorm(98, mu0,tau0)),s2tilde=rep(1,100))#

s20=4; s21=4
mu01=5; mu11=-5
Data=list(y=c(rnorm(Consts$N/2,mu01,sqrt(s20)), rnorm(Consts$N/2,mu11,sqrt(s21))))

#-- compiling the model:
model<-nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
cmodel<-compileNimble(model)

#----------------------------------------------------------------------------
#-- standarized output:
modelConf<-configureMCMC(model, print=FALSE, thin=10) # less samples!
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


mvSaved=CmodelNewMCMC$mvSamples

# works fine!
SamplerG <- nimble:::sampler_G2(model, mvSaved)
SamplerG$run()

# ploting the results: should be concentrated on 5 and -5
aux=as.matrix(SamplerG$mv)  
trunc=length(aux[1,])/3
for(i in 1:nrow(aux)){
  plot(aux[i, (trunc+1):(2*trunc)], aux[i, 1:trunc], type="h", main=sum(aux[i, 1:trunc])); readline()
}


# has and ERROR:
rSamplerG <- sampler_G2(model, mvSaved)
cSamplerG <- compileNimble(rSamplerG, project=model)


