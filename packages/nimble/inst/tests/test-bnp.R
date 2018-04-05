
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#-- TESTS FOR BNP MODELS




rm(list=ls())
library(nimble)
#nimbleOptions(allowDynamicIndexing = TRUE)
#source("./packages/nimble/R/BNP_distributions.R")
#source("./packages/nimble/R/BNP_samplers.R")


#--------------------------------------------------
# no deterministic nodes
Code=nimbleCode(
  {
    for(i in 1:N3){
      thetatilde[i] ~ dnorm(mean=mu0, var=tau20) 
      #s2tilde[i] ~ dinvgamma(shape=a0, scale=b0) 
    }
    xi[1:N2] ~ dCRP(conc=conc, size=N2)
    
    for(i in 1:N){
      theta[i] <- thetatilde[xi[i]]
      y[i] ~ dnorm(theta[i], var=2)#s2tilde[xi[i]]
    }
    conc<-1;a0<-1 ; b0<- 0.5; mu0<-0; tau20<-40; 
  }
)


conc<-1; a0<-1; b0<-0.5; mu0<-0; tau20<-40
Consts=list(N=50, N2=50, N3=50)
set.seed(1)
aux=sample(1:10, size=Consts$N2, replace=TRUE)
Inits=list(xi=aux, thetatilde=rnorm(Consts$N3, mu0, sqrt(tau20))) #, s2tilde=rinvgamma(Consts$N3, a0, b0)

s20=4; s21=4
mu01=5; mu11=-5
Data=list(y=c(rnorm(Consts$N/2,mu01,sqrt(s20)), rnorm(Consts$N/2,mu11,sqrt(s21))))

model<-nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
cmodel<-compileNimble(model)

# no deterministic nodes
Code=nimbleCode(
  {
    for(i in 1:N3){
      #thetatilde[i] ~ dnorm(mean=mu0, var=tau20) 
      s2tilde[i] ~ dgamma(shape=a0, scale=b0) 
    }
    xi[1:N2] ~ dCRP(conc=conc, size=N2)
    
    for(i in 1:N){
      #theta[i] <- thetatilde[xi[i]]
      s2[i ]<- s2tilde[xi[i]]
      y[i] ~ dnorm(0 , var=s2[i])#
    }
    conc<-1;a0<-1 ; b0<- 0.5; 
  }
)


conc<-1; a0<-1; b0<-0.5;
Consts=list(N=50, N2=50, N3=50)
set.seed(1)
aux=sample(1:10, size=Consts$N2, replace=TRUE)
Inits=list(xi=aux, s2tilde=rgamma(Consts$N3,1,1)) #, s2tilde=rinvgamma(Consts$N3, a0, b0)

s20=4; s21=4
mu01=5; mu11=-5
Data=list(y=c(rnorm(Consts$N/2,mu01,sqrt(s20)), rnorm(Consts$N/2,mu11,sqrt(s21))))

#-- compiling the model:
model<-nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
cmodel<-compileNimble(model)

modelConf<-configureMCMC(model, print=TRUE)
modelConf$printSamplers(c("xi"))
modelMCMC=buildMCMC(modelConf)

#-- compiling the sampler
CmodelNewMCMC=compileNimble(modelMCMC, project=model,
                            resetFunctions=TRUE, showCompilerOutput = TRUE)


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

#-- BUGS definition of the model: no deterministic nodes= works fine
Code=nimbleCode(
  {
    for(i in 1:N3){
      thetatilde[i] ~ dnorm(mean=mu0, var=tau20) 
      s2tilde[i] ~ dinvgamma(shape=a0, scale=b0) 
    }
    xi[1:N2] ~ dCRP(conc, size=N2)
    
    for(i in 1:N){
      y[i] ~ dnorm(thetatilde[xi[i]], var=s2tilde[xi[i]])#
    }
    conc<-1;mu0<-0; tau20<-40; a0<-1; b0<-0.5; 
  }
)

Consts=list(N=50, N2=50, N3=50)

conc<-1; a0<-1; b0<-0.5; mu0<-0; tau20<-40
set.seed(1)
aux=sample(1:10, size=Consts$N2, replace=TRUE)
Inits=list(xi=aux, thetatilde=rnorm(Consts$N3, mu0, sqrt(tau20)),s2tilde=rinvgamma(Consts$N3, shape=a0, scale=b0))#list(xi=aux, thetatilde=c(10,-10,rnorm(98, mu0,tau0)),s2tilde=rep(1,100))#

s20=4; s21=4
mu01=5; mu11=-5
Data=list(y=c(rnorm(Consts$N/2,mu01,sqrt(s20)), rnorm(Consts$N/2,mu11,sqrt(s21))))

#-- compiling the model:
model<-nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
cmodel<-compileNimble(model)

#-- MCMC configuration:
modelConf<-configureMCMC(model, print=FALSE)
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
xisamples=samples[, (N+1):(2*N)]#samples[, 1+(N+1):(2*N)]
s2samples=samples[, 1:N]
thetasamples=samples[, (N+1):(2*N)]
s2i=t(sapply(1:nsave, function(i)s2samples[i, xisamples[i,]]))
thetai=t(sapply(1:nsave, function(i)thetasamples[i, xisamples[i,]]))


apply(xisamples, 1, table)
sapply(1:10, function(i)table(thetasamples[xisamples[i,]]))

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
      thetaNew[j,1]=rnorm(1,mu0, sqrt(tau20))
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



#----------------------------------------------------------------------------
#-- standarized output:
modelConf<-configureMCMC(model, print=FALSE, thin=100) # less samples!
modelConf$printSamplers(c("xi"))
#modelConf$removeSamplers(c("xi"), print="TRUE")
#modelConf$addSampler(c("xi"), type="sampler_MarginalizedG_general")
#modelConf$printSamplers(c("xi"))

modelMCMC=buildMCMC(modelConf)
#-- compiling the sampler
CmodelNewMCMC=compileNimble(modelMCMC, project=model,
                            resetFunctions=TRUE, showCompilerOutput = TRUE)
#-- MCMC samples
set.seed(1)
nsave=1000
t1=proc.time()
CmodelNewMCMC$run(nsave)
proc.time()-t1


mvSaved=modelMCMC$mvSamples
#mvSaved <- CmodelNewMCMC$mvSamples
#samples=as.matrix(cmvSaved)
#SamplerG <- nimble:::sampler_G(model, mvSaved, target, varNames, rndconc)#nimble:::sampler_G(model, mvSaved, target, varNames)## 
SamplerG <- sampler_G2(model, mvSaved)#nimble:::sampler_G3(model, mvSaved)
cSamplerG <- compileNimble(SamplerG, project = model)
cSamplerG$run()
aux=as.matrix(cSamplerG$mv)  ## the mv object is accessed here

trunc=length(aux[1,])/3
for(i in 1:nrow(aux)){
  plot(aux[i, (2*trunc+1):(3*trunc)], aux[i, 1:trunc], type="h", main=sum(aux[i, 1:trunc]),
       xlim=c(-10,10)); readline()
}
for(i in 1:nrow(aux)){
  plot(aux[i, (trunc+1):(2*trunc)], aux[i, 1:trunc], type="h", main=sum(aux[i, 1:trunc])); readline()
}

# less parameters in teh model
trunc=length(aux[1,])/2
for(i in 1:nrow(aux)){
  plot(aux[i, (trunc+1):(2*trunc)], aux[i, 1:trunc], type="h", main=sum(aux[i, 1:trunc]),
       xlim=c(-10,10)); readline()
}


#--------------------------------------------------

Code=nimbleCode(
  {
    for(i in 1:N3){
      thetatilde[i] ~ dnorm(mean=mu0, var=tau20) 
      #s2tilde[i] ~ dinvgamma(shape=a0, scale=b0) 
    }
    xi[1:N2] ~ dCRP(conc, size=N2)
     
    for(i in 1:N){
      #theta[i] <- thetatilde[xi[i]]
      y[i] ~ dnorm(thetatilde[xi[i]] , var=2)#s2tilde[xi[i]]
    }
    conc<-1;a0<-1 ; b0<- 0.5; mu0<-0; tau20<-40; 
  }
)



conc<-1; a0<-1; b0<-0.5; mu0<-0; tau20<-40
Consts=list(N=50, N2=50, N3=50)
set.seed(1)
aux=sample(1:10, size=Consts$N2, replace=TRUE)
Inits=list(xi=aux, thetatilde=rnorm(Consts$N3, mu0, sqrt(tau20))) #, s2tilde=rinvgamma(Consts$N3, a0, b0)

s20=4; s21=4
mu01=5; mu11=-5
Data=list(y=c(rnorm(Consts$N/2,mu01,sqrt(s20)), rnorm(Consts$N/2,mu11,sqrt(s21))))

#-- compiling the model:
model<-nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
cmodel<-compileNimble(model)


Code=nimbleCode(
  {
    for(i in 1:N3){
      lambdatilde[i] ~ dgamma(1,1) 
      #s2tilde[i] ~ dinvgamma(shape=a0, scale=b0) 
    }
    xi[1:N2] ~ dCRP(conc=conc, size=N2)
    
    for(i in 1:N){
      #theta[i] <- thetatilde[xi[i]]
      y[i] ~ dpois(lambdatilde[xi[i]])
    }
    conc<-1
  }
)



conc<-1
Consts=list(N=50, N2=50, N3=50)
set.seed(1)
aux=sample(1:10, size=Consts$N2, replace=TRUE)
Inits=list(xi=aux, lambdatilde=rgamma(Consts$N3,1,1))#list(xi=aux, thetatilde=c(10,-10,rnorm(98, mu0,tau0)),s2tilde=rep(1,100))#

Data=list(y=c(rpois(Consts$N, 3)))

#-- compiling the model:
model<-nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
cmodel<-compileNimble(model)

mConf <- configureMCMC(model, print=TRUE)
mMCMC <- buildMCMC(mConf)
CmMCMC <- compileNimble(mMCMC, project=model, resetFunctions=TRUE, showCompilerOutput = TRUE)





#--------------------------------------------------
# status: done
Code=nimbleCode(
  {
    for(i in 1:N3){
      p[i,1:3] ~ ddirch(alpha=alpha0[1:3])
    }
    xi[1:N2] ~ dCRP(conc, size=N2)
    
    for(i in 1:N){
      y[i,1:3] ~ dmulti(prob=p[xi[i],1:3], size=3)
    }
    conc<-1
  }
)


conc<-1
Consts=list(N=50, N2=50, N3=50, alpha0=c(1,1,1) )
set.seed(1)
p0 <- matrix(0, ncol=3, nrow=Consts$N)
y0 <- matrix(0, ncol=3, nrow=Consts$N)
for(i in 1:Consts$N){
  p0[i,]=rdirch(1, c(1, 1, 1))
  y0[i,] = rmulti(1, prob=c(0.3,0.3,0.4), size=3)
}
aux=sample(1:10, size=Consts$N2, replace=TRUE)
Inits=list(p=p0, xi=aux)
Data=list(y=y0)
model<-nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
modelConf<-configureMCMC(model, print=TRUE)
# conjugate_ddirch_dmulti sampler: p[1:3]



#--------------------------------------------------
# model with random concentration parameter
#-- BUGS definition of the model
Code=nimbleCode(
  {
    for(i in 1:N3){
      thetatilde[i] ~ dnorm(mean=mu0, var=tau20) 
      s2tilde[i] ~ dinvgamma(shape=a0, scale=b0) 
    }
    xi[1:N2] ~ dCRP(conc, size=N2)
    conc ~ dgamma(1, 1)
    
    for(i in 1:N){
      theta[i] <- thetatilde[xi[i]]
      s2[i] <- s2tilde[xi[i]]
      y[i] ~ dnorm(theta[i], var=s2[i])
    }
    a0<-1; b0<-0.5; mu0<-0; tau20<-40
  }
)

conc<-1; a0<-1; b0<-0.5; mu0<-0; tau20<-40
Consts=list(N=50, N2=50, N3=50)
set.seed(1)
aux=sample(1:10, size=Consts$N2, replace=TRUE)
Inits=list(xi=aux, thetatilde=rnorm(Consts$N3, mu0, sqrt(tau20)),
           s2tilde=rinvgamma(Consts$N3, shape=a0, scale=b0),
           conc=1)

s20=4; s21=4
mu01=5; mu11=-5
Data=list(y=c(rnorm(Consts$N/2,mu01,sqrt(s20)), rnorm(Consts$N/2,mu11,sqrt(s21))))

#-- compiling the model:
model<-nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
cmodel<-compileNimble(model)

modelConf<-configureMCMC(model, print=FALSE, thin=100) # less samples!
modelConf$addMonitors('xi')
modelConf$removeSamplers(c("conc"), print="TRUE")
modelConf$addSampler(c("conc"), type="sampler_Augmented_BetaGamma")
modelConf$printSamplers(c("conc"))

modelMCMC=buildMCMC(modelConf)
#-- compiling the sampler
CmodelNewMCMC=compileNimble(modelMCMC, project=model,
                            resetFunctions=TRUE, showCompilerOutput = TRUE)
#-- MCMC samples
set.seed(1)
nsave=1000
t1=proc.time()
CmodelNewMCMC$run(nsave)
proc.time()-t1


mvSaved=modelMCMC$mvSamples

SamplerG <- sampler_G2(model, mvSaved)#nimble:::sampler_G3(model, mvSaved)
cSamplerG <- compileNimble(SamplerG, project = model)
cSamplerG$run()
aux=as.matrix(cSamplerG$mv)  ## the mv object is accessed here

trunc=length(aux[1,])/3
for(i in 1:nrow(aux)){
  plot(aux[i, (trunc+1):(2*trunc)], aux[i, 1:trunc], type="h", main=sum(aux[i, 1:trunc])); readline()
}


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# big model:
codeBNP <- nimbleCode({
  for(i in 1:nStudies) {
    y[i] ~ dbin(size = nStudies, prob = q[i])#dbin(size = m[i], prob = q[i])
    x[i] ~ dbin(size = nStudies, prob = p[i])#dbin(size = n[i], prob = p[i])
    q[i] <- expit(theta + gamma[i])
    p[i] <- expit(gamma[i])
    gamma[i] ~ dnorm(mu[i], var = tau[i])
    mu[i] <- muTilde[xi[i]]
    tau[i] <- tauTilde[xi[i]]
  }
  for(i in 1:nStudies) {
    muTilde[i] ~ dnorm(mu0, sd = sd0)
    tauTilde[i] ~ dinvgamma(a0, b0)
  }
  xi[1:nStudies] ~ dCRP(conc, size = nStudies)
  conc <- 1 #~ dgamma(1, 1)
  mu0 <- 0 #~ dflat()
  sd0 <- 10 #~ dunif(0, 100)
  a0 <- 1 #~ dunif(0, 100)
  b0 <- 1 #~ dunif(0, 100)
  theta <- 0 #~ dflat()
})

Consts=list(nStudies=10)
set.seed(1)
Inits=list(gamma=rep(1,10),
           muTilde=rep(1,10),
           tauTilde=rep(1,10),
           xi=rep(1,10))

Data=list(y=rbinom(10, 10, 0.5), x=rbinom(10, 10, 0.5))

model<-nimbleModel(codeBNP, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
cmodel<-compileNimble(model)


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# other big model
library(emplik)
data(myeloma)

n <- nrow(myeloma)
time <-  myeloma[ , 1]
vstatus <- myeloma[ , 2]  # 0 = alive (i.e., censored)
alive <- vstatus == 0
cens_time <- rep(NA, n)
cens_time[alive] <- time[alive]
cens_time[!alive] <- Inf
time[alive] <- NA

logBUN <- myeloma[ , 3]
HGB <- myeloma[ , 4]
logBUN <- (logBUN - mean(logBUN)) / sd(logBUN)
HGB <- (HGB - mean(HGB)) / sd(HGB)

## accelerated failure time model per https://www4.stat.ncsu.edu/~ghosal/papers/PMR.pdf for Bayesian semiparametric AFT models
codeAFT <- nimbleCode({
  for(i in 1:n) {
    x[i] ~ dweib(alpha, 1+exp(lambda[i]))   # 'data' nodes
    is_cens[i] ~ dinterval(x[i], c[i])
    lambda[i] <-  inprod(Z[i, 1:p], delta[1:p]) + eta[i] 
    eta[i] <- etaTilde[xi[i]]
  }
  xi[1:n] ~ dCRP(conc, size = n)
  conc ~ dgamma(1, 1)
  for(i in 1:n)
    etaTilde[i] ~ dunif(b0, B0)
  alpha ~ dunif(a0, A0)
  for(j in 1:p)
    delta[j] ~ dflat()
})

constants = list(b0 = -10, B0 = 10, a0 = 0.1, A0 = 10, p = 2, n = n, c
                 = cens_time, Z = cbind(logBUN, HGB))
data = list(is_cens = as.numeric(alive), x = time)
xInit <- rep(NA, n)
xInit[alive] <- cens_time[alive] + 10
inits = list(alpha = 1, delta = c(0, 0), conc = 1, etaTilde = runif(n,
                                                                    constants$b0, constants$B0),
             xi = sample(1:3, n, replace = TRUE), x = xInit)

model <- nimbleModel(codeAFT, constants = constants, data = data, inits = inits)
conf = configureMCMC(model, print = TRUE)
mcmc = buildMCMC(conf)


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
Code=nimbleCode(
  {
    for(i in 1:N){
      s2tilde[i] ~ dinvgamma(shape=a0, scale=b0) 
    }
    xi[1:N] ~ dCRP(conc, size=N)
    conc ~ dgamma(1, 1)
    
    for(i in 1:N){
      s2[i] <- s2tilde[xi[i]]
      y[i] ~ dnorm(mu1 + mu2, s2[i])
    }
    mu1 ~ dnorm(0,1)
    mu2 ~ dnorm(0,1)
    a0<-1; b0<-0.5 
  }
)

Code=nimbleCode(
  {
    for(i in 1:N){
      #thetatilde[i] ~ dnorm(mean=mu0, var=tau20) 
      s2tilde[i] ~ dinvgamma(shape=a0, scale=b0) 
    }
    xi[1:N] ~ dCRP(conc, size=N)
    conc ~ dgamma(1, 1)
    
    for(i in 1:N){
      #theta[i] <- thetatilde[xi[i]]
      #s2[i] <- s2tilde[xi[i]]
      y[i] ~ dnorm(mu1 + mu2, exp(s2tilde[xi[i]]))#y[i] ~ dnorm(theta[i], var=s2[i])
    }
    mu1 ~ dnorm(0,1)
    mu2 ~ dnorm(0,1)
    a0<-1; b0<-0.5 # tau20<-40 #; mu0<-0;
  }
)

conc<-1; a0<-1; b0<-0.5; mu0<-0; tau20<-40
Consts=list(N=50)
set.seed(1)
aux=sample(1:10, size=Consts$N, replace=TRUE)
Inits=list(xi=aux, #thetatilde=rnorm(Consts$N3, mu0, sqrt(tau20)),
           s2tilde=rinvgamma(Consts$N, shape=a0, scale=b0),
           conc=1, mu1=1, mu2=2)

s20=4; s21=4
mu01=5; mu11=-5
Data=list(y=c(rnorm(Consts$N/2,mu01,sqrt(s20)), rnorm(Consts$N/2,mu11,sqrt(s21))))

#-- compiling the model:
model<-nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
cmodel<-compileNimble(model)


Code=nimbleCode(
  {
    for(i in 1:N){
      thetatilde[i] ~ dnorm(mean=mu0, var=tau20) 
      #s2tilde[i] ~ dinvgamma(shape=a0, scale=b0) 
    }
    xi[1:N] ~ dCRP(conc, size=N)
    conc ~ dgamma(1, 1)
    
    for(i in 1:N){
      #theta[i] <- thetatilde[xi[i]]
      #s2[i] <- s2tilde[xi[i]]
      y[i] ~ dnorm( mu1+ thetatilde[xi[i]] , exp(s2))#y[i] ~ dnorm(theta[i], var=s2[i])
    }
    mu1 ~ dnorm(0,1)
    s2 ~ dinvgamma(1,1)
    tau20<-40 ; mu0<-0; # 
  }
)

conc<-1; a0<-1; b0<-0.5; mu0<-0; tau20<-40
Consts=list(N=50)
set.seed(1)
aux=sample(1:10, size=Consts$N, replace=TRUE)
Inits=list(xi=aux, #thetatilde=rnorm(Consts$N3, mu0, sqrt(tau20)),
           thetatilde=rnorm(Consts$N, 0,1),
           conc=1, mu1=1, s2=2)

s20=4; s21=4
mu01=5; mu11=-5
Data=list(y=c(rnorm(Consts$N/2,mu01,sqrt(s20)), rnorm(Consts$N/2,mu11,sqrt(s21))))

#-- compiling the model:
model<-nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
cmodel<-compileNimble(model)


# OK: y[i] ~ dnorm(mu1 + mu2, sigma[xi[i]])
# OK: y[i] ~ dnorm(mu1 + mu2, exp(sigma[xi[i]]))
y[i] ~ dnorm(mu + theta[xi[i]], sigma)


#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#-- stick_breaking representation and BUGS Code

#-- stick breaking nimble function:

#' The Stick Breaking function
#'
#'   Based on the stick breaking construction, weights are computed 
#'   with the argument vector.
#' 
#' @name StickBreakingFunction 
#' 
#' @param z vector argument.
#' @param log logical; if TRUE, weights are returned on the log scale.
#' @author Claudia Wehrhahn
#' @details values in \code{z} have to be between \deqn{[0,1)]}. If one of the components 
#' is equal to 1, then the returned vector has length smaller than \code{z}. If one of the
#' components is smaller than 0 or greater than 1, NaNs are returned.
#' @references Sethuraman, J. (1994) A constructive definition of Dirichlet 
#' priors. \emph{Statistica sinica}, 639-650.
#' @examples
#' z <- rbeta(5, 1,1)
#' stick_breaking(z)
#' 
#' cstick_breaking <- compileNimble(stick_breaking)
#' cstick_breaking(z)
NULL

#' @rdname StickBreakingFunction
#' @export
stick_breaking=nimbleFunction(
  run=function(z=double(1),
               log=integer(0, default=0)){ # z values must be different of 1, otherwise the process is truncated to a smaller value tha N; never use z[N]!
    returnType(double(1))
    
    cond <- sum(z < 0)
    if(cond > 0){
      print('values in vector z have to be in [0,1)')
      cond1 <- 1
    }else{
      cond1 <- 0
    }
    cond <- sum(z > 1)
    if(cond > 0){
      print('values in vector z have to be in [0,1)')
      cond2 <- 1
    }else{
      cond2 <- 0
    }
    
    cond <- sum(z == 1)
    if(cond > 0){ # el vector de probabilidades es mas chico
      print('length of returned vector is less than length(z)')
      N <- 1
      while(z[N] != 1){
        N <- N + 1
      }
    }else{
      N<-length(z)
    }
    
    if(cond1 + cond2 > 0){
      return(c(NaN))
    }else{
      x<-numeric(N) # returned vector
      tmp<-0#1
      
      x[1]<-log(z[1]) #z[1]
      for(i in 2:(N-1)){
        tmp=tmp+log(1-z[i-1]) #tmp*(1-z[i-1])
        x[i]<-log(z[i])+tmp #z[i]*tmp
      }
      x[N]<-tmp+log(1-z[N-1])#tmp*(1-z[N-1])
      if(log) return(x)
      else return(exp(x))
    }
  }
)

z <- rbeta(5, 1,1)
stick_breaking(z)

cstick_breaking <- compileNimble(stick_breaking)
cstick_breaking(z)

#-- BUGS code:
Code=nimbleCode(
  {
    for(i in 1:Trunc){
      thetatilde[i] ~ dnorm(mean=mu0, var=tau20) 
      s2tilde[i] ~ dinvgamma(shape=a0, scale=b0) 
      z[i] ~ dbeta(1, conc) # could be dbeta(a_i,b_i)
    }
    w[1:Trunc] <- stick_breaking(z[1:Trunc])
    
    for(i in 1:N){
      xi[i] ~ dcat(w[1:Trunc])
      theta[i] <- thetatilde[xi[i]]
      s2[i] <- s2tilde[xi[i]]
      y[i] ~ dnorm(theta[i], var=s2[i])
    }
    conc<-1; a0<-1; b0<-0.5; mu0<-0; tau20<-40
  }
)

conc<-1; a0<-1; b0<-0.5; mu0<-0; tau20<-40
Consts=list(N=50, Trunc=25)
set.seed(1)
Inits=list(thetatilde=rnorm(Consts$Trunc, mu0, sqrt(tau20)),
           s2tilde=rinvgamma(Consts$Trunc, shape=a0, scale=b0),
           z=rbeta(Consts$Trunc, 1, conc),
           xi=sample(1:10, size=Consts$N, replace=TRUE))
s20=4; s21=4
mu01=5; mu11=-5
Data=list(y=c(rnorm(Consts$N/2,mu01,sqrt(s20)), rnorm(Consts$N/2,mu11,sqrt(s21))))

#-- compiling the model:
model=nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)#, dimensions=list(theta=Consts$N))
cmodel=compileNimble(model)

#-- MCMC configuration:
modelConf=configureMCMC(model, print=FALSE)
modelConf$printSamplers(c('xi[1]', 'z[1]'))
modelMCMC2=buildMCMC(modelConf)

#-- compiling the sampler
CmodelNewMCMC=compileNimble(modelMCMC2, project=model,
                            resetFunctions=TRUE, showCompilerOutput = TRUE)

#-- MCMC samples
set.seed(1)
nsave=100
t1=proc.time()
CmodelNewMCMC$run(nsave)
proc.time()-t1

#-- results:
samples=as.matrix(CmodelNewMCMC$mvSamples)
s2tildepost=samples[, 1:25]
thetatildepost=samples[, 26:50]
Zpost=samples[, 51:75]
Tr=Consts$Trunc
Wpost=t(apply(Zpost, 1, function(x)c(x[1], x[2:(Tr-1)]*cumprod(1-x[1:(Tr-2)]), cumprod(1-x[1:(Tr-1)])[N=Tr-1])))

# grid:
ngrid=102
grid=seq(-10, 25,len=ngrid)

# posterior preddictives
predSB=matrix(0, ncol=ngrid, nrow=nsave)
for(i in 1:nsave){
  predSB[i, ]=sapply(1:ngrid, function(j)sum(Wpost[i, ]*dnorm(grid[j], thetatildepost[i,],sqrt(s2tildepost[i,]))))
}

hist(Data$y, freq=FALSE, xlim=c(min(grid), max(grid)))
points(grid, apply(predSB, 2, mean), col="blue", type="l", lwd=2)




##########################################################################
# testing sampler_G
##########################################################################



Code=nimbleCode(
  {
    mu0 ~ dnorm(0,1)
    for(i in 1:N3){
      thetatilde[i] ~ dnorm(mean=mu0, var=40) 
      s2tilde[i] ~ dinvgamma(shape=1, scale=0.5) 
    }
    xi[1:N2] ~ dCRP(conc1, size=N2)
    conc1 <- beta + alpha + 3
    alpha ~ dgamma(1,3)
    beta ~ dgamma(4,5)
    
    for(i in 1:N){
      s2[i] <- s2tilde[xi[i]]
      y[i] ~ dnorm(thetatilde[xi[i]], var=s2[i])
    }
  }
)

conc<-1; a0<-1; b0<-0.5; mu0<-0; tau20<-40
Consts=list(N=50, N2=50, N3=50)
set.seed(1)
aux=sample(1:10, size=Consts$N2, replace=TRUE)
Inits=list(xi=aux, thetatilde=rnorm(Consts$N3, mu0, sqrt(tau20)),
           s2tilde=rinvgamma(Consts$N3, a0, b0),
           mu0 =1, alpha=1, beta=1)

s20=4; s21=4
mu01=5; mu11=-5
Data=list(y=c(rnorm(Consts$N/2,mu01,sqrt(s20)), rnorm(Consts$N/2,mu11,sqrt(s21))))

model<-nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
cmodel<-compileNimble(model)

modelConf<-configureMCMC(model, print=TRUE, thin=100) # less samples!
modelConf$printSamplers(c("xi"))
#modelConf$removeSamplers(c("xi"), print="TRUE")
#modelConf$addSampler(c("xi"), type="sampler_MarginalizedG_general")
#modelConf$printSamplers(c("xi"))

modelMCMC=buildMCMC(modelConf)
#-- compiling the sampler
CmodelNewMCMC=compileNimble(modelMCMC, project=model,
                            resetFunctions=TRUE, showCompilerOutput = TRUE)
#-- MCMC samples
set.seed(1)
nsave=1000
t1=proc.time()
CmodelNewMCMC$run(nsave)
proc.time()-t1


mvSaved=modelMCMC$mvSamples
#mvSaved <- CmodelNewMCMC$mvSamples
#samples=as.matrix(cmvSaved)
#SamplerG <- nimble:::sampler_G(model, mvSaved, target, varNames, rndconc)#nimble:::sampler_G(model, mvSaved, target, varNames)## 
SamplerG <- sampler_G2(model, mvSaved)#nimble:::sampler_G3(model, mvSaved)
cSamplerG <- compileNimble(SamplerG, project = model)
cSamplerG$run()
aux=as.matrix(cSamplerG$mv)  ## the mv object is accessed here

trunc=length(aux[1,])/3
for(i in 1:nrow(aux)){
  plot(aux[i, (2*trunc+1):(3*trunc)], aux[i, 1:trunc], type="h", main=sum(aux[i, 1:trunc]),
       xlim=c(-10,10)); readline()
}
for(i in 1:nrow(aux)){
  plot(aux[i, (trunc+1):(2*trunc)], aux[i, 1:trunc], type="h", main=sum(aux[i, 1:trunc])); readline()
}

# less parameters in teh model
trunc=length(aux[1,])/2
for(i in 1:nrow(aux)){
  plot(aux[i, (trunc+1):(2*trunc)], aux[i, 1:trunc], type="h", main=sum(aux[i, 1:trunc]),
       xlim=c(-10,10)); readline()
}




##########################################################################
# conjugate samplers that are/will be added the BNP models:
##########################################################################


#--------------------------------------------------
# status: not identified by nimble
Code=nimbleCode(
  {
    p ~ dbeta(1, 1)
    for(i in 1:N){
      y[i] ~ dbin( p, 5) # dbinom(5, p)
    }
  }
)


Consts=list(N=50)
set.seed(1)
Inits=list(p=rbeta(1, 1, 1))
Data=list(y=rbinom(Consts$N, size=5, prob=0.5))
model<-nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
modelConf<-configureMCMC(model, print=TRUE)
# there are warnings and the assigned sampler is a RW sampler.

#--------------------------------------------------
# status: 
Code=nimbleCode(
  {
    l ~ dgamma(shape=1, rate=1)
    for(i in 1:N){
      y[i] ~ dgamma(shape=1, rate=l)
    }
  }
)


Consts=list(N=50)
set.seed(1)
Inits=list(l=rgamma(1, shape=1, rate=1))
Data=list(y=rgamma(Consts$N, shape=1 ,rate=1))
model<-nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
modelConf<-configureMCMC(model, print=TRUE)
# conjugate_dgamma_dgamma sampler: l



#--------------------------------------------------
# status: 
Code=nimbleCode(
  {
    l ~ dgamma(shape=1, rate=1)
    for(i in 1:N){
      y[i] ~ dexp(rate=l)
    }
  }
)


Consts=list(N=50)
set.seed(1)
Inits=list(l=rgamma(1, 1, 1))
Data=list(y=rexp(Consts$N, rate=1))
model<-nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
modelConf<-configureMCMC(model, print=TRUE)
# conjugate_dgamma_dexp sampler: l


#--------------------------------------------------
# status: done
Code=nimbleCode(
  {
    p[1:3] ~ ddirch(alpha=alpha0[1:3])
    #for(i in 1:N){
      y[1:3] ~ dmulti(prob=p[1:3], size=3)
    #}
    #alpha0 <- c(1,1,1)  
  }
)


Consts=list(N=1, alpha0=c(1,1,1)  )
set.seed(1)
Inits=list(p=rdirch(1, c(1, 1, 1)))
Data=list(y=rmulti(Consts$N, prob=c(0.3,0.3,0.4), size=3))
model<-nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
modelConf<-configureMCMC(model, print=TRUE)
# conjugate_ddirch_dmulti sampler: p[1:3]


#--------------------------------------------------
# status: done
Code=nimbleCode(
  {
    p ~ dbeta(1, 1)
    for(i in 1:N){
      y[i] ~ dbern(p)
    }
  }
)


Consts=list(N=50)
set.seed(1)
Inits=list(p=rbeta(1, 1, 1))
Data=list(y=rbinom(Consts$N, prob=0.5))
model<-nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
modelConf<-configureMCMC(model, print=TRUE)
# conjugate_dbeta_dbern sampler: p


#--------------------------------------------------
# status: done
Code=nimbleCode(
  {
    mu ~ dnorm(0, 10)
    for(i in 1:N){
      y[i] ~ dnorm(mu, 1)
    }
  }
)
Consts=list(N=50)
set.seed(1)
Inits=list(mu=rnorm(1, 0, 10))
Data=list(y=rnorm(Consts$N, 0, 1))
model<-nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
modelConf<-configureMCMC(model, print=TRUE)

#--------------------------------------------------
# status: done 
Code=nimbleCode(
  {
    lambda ~ dgamma(1,1)
    for(i in 1:N){
      y[i] ~ dpois(lambda)
    }
  }
)


Consts=list(N=50)
set.seed(1)
Inits=list(lambda=dgamma(1, 1, 1))
Data=list(y=rpois(Consts$N, 1))
model<-nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
modelConf<-configureMCMC(model, print=TRUE)
# conjugate_dgamma_dpois sampler: lambda






