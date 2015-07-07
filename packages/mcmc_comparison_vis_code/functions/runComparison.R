runComparison<-function(allModels,mcmcs,niter=10000,thin=1,MCMCdefs,standir,BUGSdir,stanNameMapsList,stanModelNamesList,plot_on=FALSE,...){
  library(nimble)
  library(rjags)
  library(rstan)
  library(xtable)
  library(ggplot2)
  library(gridExtra)
  library(igraph)
  
  x=list()
  for (i in 1:length(allModels)){
      thisBUGSdir <- if(missing(BUGSdir)) getBUGSexampleDir(allModels[i]) else BUGSdir 
    x[[i]]<-readBUGSmodel(model=allModels[i],
                          dir=thisBUGSdir,
                          returnModelComponentsOnly=TRUE)
  }
  
  mods=list()

  noConjDef <- list(noConj = quote({ configureMCMC(Rmodel, useConjugacy=FALSE) }))
  if(missing(MCMCdefs)) MCMCdefs <- noConjDef
  else MCMCdefs <- c(MCMCdefs, noConjDef)

  if(missing(stanNameMapsList)) stanNameMapsList <- NULL
  if(missing(stanModelNamesList)) stanModelNamesList <- NULL
  
  for (i in 1:length(allModels)){
      print(allModels[i])
      
      stanModelName <- stanModelNamesList[[ allModels[i] ]]
      if(is.null(stanModelName)) stanModelName <- allModels[i]
      stanNameMaps <- stanNameMapsList[[ stanModelName ]]

      if(is.null(stanNameMaps)) stanNameMaps <- list()
    suite_output <- MCMCsuite(x[[i]]$model, constants = x[[i]]$data, inits = x[[i]]$inits, 
                               MCMCs = mcmcs,makePlot=plot_on,savePlot=plot_on,niter=niter,thin=thin
                              ,summaryStats=c('mean','median','sd','CI95_low','CI95_upp','effectiveSize')
                              #change
                              ,MCMCdefs = MCMCdefs ##list(noConj = quote({ configureMCMC(Rmodel, useConjugacy=FALSE) }))
                             ,stan_model=if('stan' %in% mcmcs) paste(standir,stanModelName,'/',stanModelName,'.stan',sep="") else ""
                              , stanNameMaps = stanNameMaps, ...
    )
    mods[allModels[i]][[1]]=list(suite_output$summary,create_time_df(suite_output$timing,length(mcmcs)))
  }
  return(mods)
}
  


