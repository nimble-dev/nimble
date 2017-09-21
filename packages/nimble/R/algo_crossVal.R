library(nimble)
library(coda)
library(parallel)

crossValCalculate <- function(MCMCOut, realData){
  discrep <- sum(colMeans((MCMCOut - saveData)^2))
  return(discrep)
}

getSD <- function(MCMCout, sampNum, NbootReps, saveData){
  l <- ceiling(min(1000, (sampNum)/20)) ##length of each
  q <- sampNum - l + 1
  h <- ceiling(sampNum/l)
  blockcvValues <- numeric(NbootReps) ## block estimates of cv
  for(r in 1:NbootReps){
    randNum <- runif(h, 0, 1)
    randIndexStart <- ceiling(randNum*q)
    indexes <- cbind(randIndexStart, randIndexStart + l - 1)
    chunks <- c(apply(indexes, 1, function(x){
      seq(x[1], x[2], 1)}
    ))
    blockcvValues[r] <- crossValCalculate(MCMCout[chunks,], saveData)
  }
  return(sd(blockcvValues))
}

runCrossValidate <- function(model,
                             dataNames,
                             MCMCIter,
                             burnInProp, thin,
                             leaveOutIndex,
                             MCMCdefs=NULL,
                             returnSamples=FALSE,
                             NbootReps=200){
  # simpSetMCMCDefs <- function(Rmodel, MCMCdefs, MCMCname) {
  #   eval(MCMCdefs[[MCMCname]])
  # }
  # if(!is.null(MCMCdefs)){
    ## eval(MCMCdefs.opt2[['nimbleOpt2']])
    ## customSpec <- simpSetMCMCDefs(occ.R.model, MCMCdefs.opt2, 'nimbleOpt2')
  # }
  if(!inherits(model, "RmodelBaseClass")){
    stop("model is not an Rmodel")
  }
  if(thin == 0 | thin > MCMCIter){
    stop("thin needs to be an integer > 0 and < MCMCIter")
  }
  if(burnInProp >= 1 | burnInProp < 0){
    stop("burnInProp needs to be between 0 and 1")
  }
  testDataNames <- try(model[[dataNames]], silent=TRUE)
  if(inherits(testDataNames, "try-error")){
    stop(paste("dataNames", dataNames,
               "is not the name of the data in model"))
  } else{
    test2DataNames <- all(model$expandNodeNames(dataNames) %in%
                            model$getNodeNames(dataOnly=TRUE))
    if(test2DataNames == FALSE){
      stop(paste("dataNames", dataNames,
                 "is not the name of the data in model"))
    }
  }
  allDataNodeNames <- unique(sapply(dataNames, function(x){return(model$expandNodeNames(nodes = x))}))
  approxNumPerFold <- ceiling(length(allDataNodeNames)/k)
  leaveOutNames <- list()
  currentDataNames <- list()
  reducedDataNodeNames <- allDataNodeNames
  for(i in 1:(k-1)){
    leaveOutNames[[i]] <- sample(reducedDataNodeNames, approxNumPerFold)
    reducedDataNodeNames <- reducedDataNodeNames[-which(reducedDataNodeNames == leaveOutNames[[i]])]
    currentDataNames[[i]] <- allDataNodeNames[-which(reducedDataNodeNames == leaveOutNames[[i]])]
  }
  leaveOutNames[[k]] <- reducedDataNodeNames
  currentDataNames[[k]] <- allDataNodeNames[-which(reducedDataNodeNames == leaveOutNames[[k]])]
  
  prepModel <- model$newModel()
  calcCrossVal <- function(i,
                           leaveOutNames,
                           currentDataNames,
                           mcmcIter,
                           burnInProp,
                           returnSamples){
    saveData <- values(prepModel, leaveOutNames[[i]])
    newModel <- prepModel$newModel()
    newModel$resetData()
    values(newModel, leaveOutNames[[i]]) <- NA
    newModel$setData(currentDataNames[[i]])
    compileNimble(newModel)
    modelMCMCConf <- configureMCMC(newModel,
                                   monitors = leaveOutNames[[i]],
                                   thin = thin)
    modelMCMC <- buildMCMC(modelMCMCConf)
    C.modelMCMC <- compileNimble(modelMCMC,
                                 project = newModel,
                                 resetFunctions = (i != 1))
    C.modelMCMC$run(MCMCIter)
    MCMCout <- as.matrix(C.modelMCMC$mvSamples)[,leaveOutNames[[i]]]
    sampNum <- dim(MCMCout)[1]
    startIndex <- max(c(ceiling(burnInProp*MCMCIter), 1))
    message(paste("dropping data block", i))
    crossValAveSD <- getSD(MCMCout=MCMCout, sampNum=sampNum,
                           NbootReps=NbootReps,
                           dataDimensions= dataDimensions,
                           saveData=saveData)
    crossVal <- crossValCalculate(MCMCout[startIndex:sampNum,], saveData)
    nimble:::clearCompiled(C.modelMCMC)
    return(list(crossValAverage = crossVal,
                crossValAveSD = crossValAveSD,
                samples= if(returnSamples) MCMCout else NA))
  }
  crossValOut <- mclapply(1:numBlocks, calcCrossVal,
                          leaveOutNames,
                          currentDataNames,
                          mcmcIter,
                          burnInProp,
                          returnSamples)
  crossVal <- sum(sapply(crossValOut, function(x) x$crossValAverage),
                  na.rm=TRUE)
  crossValSD <- sqrt(sum(sapply(crossValOut,
                                function(x) x$crossValAveSD^2),
                         na.rm=TRUE))
  out <- list(crossVal=crossVal,
              crossValSD=crossValSD)
  if(!returnSamples){
    out$samples <- NULL
  }
  else{
    out$samples <- lapply(crossValOut, function(x) x$samples)
  }
  
  return(out)
}
