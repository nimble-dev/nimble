
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

#-- MCMC configuration:
modelConf<-configureMCMC(model, print=FALSE)
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
#target='xi'
#varNames=c('thetatilde', 's2tilde') #  in this order!!!!
#rndconc='FALSE'


#SamplerG <- nimble:::sampler_G(model, mvSaved, target, varNames, rndconc)#nimble:::sampler_G(model, mvSaved, target, varNames)## 
SamplerG <- sampler_G2(model, mvSaved)#nimble:::sampler_G3(model, mvSaved)
cSamplerG <- compileNimble(SamplerG, project = model)
cSamplerG$run()
aux=as.matrix(cSamplerG$mv)  ## the mv object is accessed here

trunc=length(aux[1,])/3
for(i in 1:nrow(aux)){
  plot(aux[i, (trunc+1):(2*trunc)], aux[i, 1:trunc], type="h", main=sum(aux[i, 1:trunc])); readline()
}



#--------------------------------------------------
# model with random concentration parameter
#-- BUGS definition of the model
Code=nimbleCode(
  {
    for(i in 1:N3){
      thetatilde[i] ~ dnorm(mean=mu0, var=tau20) 
      s2tilde[i] ~ dinvgamma(shape=a0, scale=b0) 
    }
    xi[1:N2] ~ dCRP(conc)
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

#-- standarized output:
modelConf<-configureMCMC(model, print=FALSE, thin=10) # less samples!

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
nsave=100
t1=proc.time()
CmodelNewMCMC$run(nsave)
proc.time()-t1


mvSaved=CmodelNewMCMC$mvSamples

SamplerG <- nimble:::sampler_G2(model, mvSaved)
SamplerG$run()
aux=as.matrix(SamplerG$mv)  ## the mv object is accessed here

trunc=length(aux[1,])/3
for(i in 1:nrow(aux)){
  plot(aux[i, (trunc+1):(2*trunc)], aux[i, 1:trunc], type="h", main=sum(aux[i, 1:trunc])); readline()
}


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

