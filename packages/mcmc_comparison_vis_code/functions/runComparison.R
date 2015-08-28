runComparison<-function(allModels,mcmcs,niter=10000,thin=1,MCMCdefs,standir, BUGSdir, stanNameMapsList,stanModelNamesList, stanDataFileList, stanCodeFileList, plot_on=FALSE,...){
  library(nimble)
  library(rjags)
  library(rstan)
  library(xtable)
  library(ggplot2)
  library(gridExtra)
  library(igraph)

  if(is.character(allModels)) {
      x=list()
      for (i in 1:length(allModels)){
          thisBUGSdir <- if(missing(BUGSdir)) getBUGSexampleDir(allModels[i]) else BUGSdir 
          x[[i]]<-readBUGSmodel(model=allModels[i],
                                dir=thisBUGSdir,
                                returnModelComponentsOnly=TRUE)
      }
  } else {
      if(!is.list(allModels)) stop('allModels must be a list if it is not a vector of BUGS example names')
      if(!is.list(allModels[[1]])) x <- list(allModels)
      else x <- allModels
      modelCounter <- 1
      allModels <- unlist(lapply(x, function(z) {ans <- z[['name']]; if(is.null(ans)) {ans <- paste0('model',modelCount); modelCounter <<- modelCounter + 1}; ans} ) )
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

      if(missing(stanCodeFileList)) stanCodeFile <- stanModelName
      else stanCodeFile <- stanCodeFileList[[ allModels[i] ]]
      
      if(missing(stanDataFileList)) stanDataFile <- NULL
      else stanDataFile <- file.path(standir, stanModelName, stanDataFileList[[ allModels[i] ]]) ## it's ok if this is NULL
      
      if(is.null(stanNameMaps)) stanNameMaps <- list()
    suite_output <- MCMCsuite(x[[i]]$model, constants = x[[i]]$data, inits = x[[i]]$inits, 
                               MCMCs = mcmcs,makePlot=plot_on,savePlot=plot_on,niter=niter,thin=thin
                              ,summaryStats=c('mean','median','sd','CI95_low','CI95_upp','effectiveSize')
                              #change
                              ,MCMCdefs = MCMCdefs ##list(noConj = quote({ configureMCMC(Rmodel, useConjugacy=FALSE) }))
                             ,stan_model=if('stan' %in% mcmcs) file.path(standir,stanModelName,paste0(stanCodeFile,'.stan')) else ""
                            , stanNameMaps = stanNameMaps
                            , stan_data = stanDataFile, ...
    )
    mods[allModels[i]][[1]]=list(suite_output$summary,create_time_df(suite_output$timing,length(mcmcs)))
  }
  return(mods)
}
  


