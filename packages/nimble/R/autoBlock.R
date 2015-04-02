



#' Automated parameter blocking for efficient MCMC sampling
#'
#' Automated blocking description
#'
#' @details		
#' Automated blocking details
#' 
#' @param Rmodel
#'
#' @return info
#'  
#' @examples
#' ab <- autoBlock(Rmodel)
#' 
#' @author Daniel Turek
#' @export
autoBlock <- function(Rmodel,
                      autoIt = 20000,
                      run = list('all', 'default'),
                      verbose = FALSE,
                      round = TRUE ) {
    if(autoIt < 10000) stop('Minimum auto-blocking iterations is 10,000')
    control <- list(niter=autoIt, verbose=verbose)
    ab <- autoBlockClass(Rmodel, control)
    if(!'auto' %in% run) run <- c(run, 'auto')  ## always use 'autoBlock' routine
    ab$run(run)
    lastAutoInd <- max(grep('^auto', ab$naming))   ## index of final 'auto' iteration
    lastAutoGrouping <- ab$grouping[[lastAutoInd]]  ## grouping of final 'auto' iteration
    nonTrivialGroups <- lastAutoGrouping[unlist(lapply(lastAutoGrouping, function(x) length(x)>1))]
    abList <- list(ab)
    names(abList)[1] <- 'model'
    df <- createDFfromABlist(abList, autoIt)
    dfmin <- reduceDF(df, round = round)
    cat('\nAuto-Blocking summary:\n')
    print(dfmin)
    if(length(nonTrivialGroups) > 0) {
        cat('\nAuto-Blocking converged on the node groupings:\n')
        for(i in seq_along(nonTrivialGroups)) {
            group <- nonTrivialGroups[[i]]
            cat(paste0('[', i, '] '))
            cat(paste0(group, collapse = ', '))
            cat('\n')
        }
    } else cat('\nAuto-Blocking converged on all scalar (univariate) sampling\n')
    cat('\n')
    ## create a new MCMC spec with the autoBlock groupings:
    spec <- configureMCMC(Rmodel, nodes = NULL)
    for(nodeGroup in lastAutoGrouping) addSamplerToSpec(Rmodel, spec, nodeGroup)
    retList <- list(summary=dfmin, autoGroups=nonTrivialGroups, spec=spec)
    return(invisible(retList))
}



autoBlockModel <- setRefClass(
    Class = 'autoBlockModel',
    fields = list(
        Rmodel = 'ANY',
        md = 'ANY',
        scalarNodeVector = 'character',
        nodeGroupScalars = 'list',
        nodeGroupAllBlocked = 'list',
        monitorsVector = 'character',
        initialMCMCspec = 'ANY'
        ),
    methods = list(
        initialize = function(Rmodel_orig) {
            require(nimble)
            md <<- Rmodel_orig$modelDef
            Rmodel <<- md$newModel()
            nimCopy(from = Rmodel_orig, to = Rmodel, logProb = TRUE)
            for(var in ls(Rmodel_orig$isDataEnv)) Rmodel$isDataEnv[[var]] <<- Rmodel_orig$isDataEnv[[var]]  ## copies data flags to the new model
            scalarNodeVector <<- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE, returnScalarComponents=TRUE)
            nodeGroupScalars <<- lapply(scalarNodeVector, function(x) x)
            nodeGroupAllBlocked <<- list(scalarNodeVector)
            stochNodeVector <- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE, returnScalarComponents=FALSE)
            monitorsVector <<- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE)
        },
        ## here is where the initial MCMC spec is created, for re-use -- for new version
        createInitialMCMCspec = function(runList) {
            initialMCMCspec <<- configureMCMC(Rmodel)
            nInitialSamplers <- length(initialMCMCspec$samplerSpecs)
            initialMCMCspec$addSampler('RW',       control = list(targetNode  = scalarNodeVector[1]), print=FALSE)  ## add one RW sampler
            initialMCMCspec$addSampler('RW_block', control = list(targetNodes = scalarNodeVector[1]), print=FALSE)  ## add one RW_block sampler
            addCustomizedSamplersToInitialMCMCspec(runList)
            initialMCMCspec$addMonitors(monitorsVector, print=FALSE)
            RinitialMCMC <- buildMCMC(initialMCMCspec)
            Cmodel <- compileNimble(Rmodel)
            CinitialMCMC <- compileNimble(RinitialMCMC, project = Rmodel)   ## (new version) yes, we need this compileNimble call -- this is the whole point!
            initialMCMCspec$setSamplers(1:nInitialSamplers, print=FALSE)  ## important for new version: removes all news samplers added to initial MCMC spec
        },
        addCustomizedSamplersToInitialMCMCspec = function(runListCode) {
            if(is.list(runListCode)) { lapply(runListCode, function(el) addCustomizedSamplersToInitialMCMCspec(el)); return() }
            if(is.call(runListCode)) {
                if(is.call(runListCode[[1]]) && length(runListCode[[1]])==3 && runListCode[[1]][[3]]=='addSampler') {
                    runListCode[[1]][[2]] <- as.name('initialMCMCspec')
                    eval(substitute(RUNLISTCODE, list(RUNLISTCODE=runListCode)))
                    return()
                }
                lapply(runListCode, function(el) addCustomizedSamplersToInitialMCMCspec(el))
                return()
            }
        },
        createGroups = function(listOfBlocks = list()) {
            listOfBlocks <- lapply(listOfBlocks, function(blk) Rmodel$expandNodeNames(blk, returnScalarComponents=TRUE))
            nodes <- scalarNodeVector
            nodes <- setdiff(nodes, unlist(listOfBlocks))
            nodeList <- lapply(nodes, function(x) x)
            for(ng in listOfBlocks) nodeList[[length(nodeList)+1]] <- ng
            return(nodeList)
        }
    )
)



autoBlockParamDefaults <- function() {
    list(
        makePlots = FALSE,
        niter = 20000,
        verbose = FALSE
        )
}


autoBlockClass <- setRefClass(

    Class = 'autoBlockClass',

    fields = list(
        
        ## special
        abModel = 'ANY',
        it = 'numeric',

        ## overall control
        makePlots = 'logical',
        niter = 'numeric',
        verbose = 'logical',

        ## persistant lists of historical data
        naming = 'list',
        candidateGroups = 'list',
        grouping = 'list',
        groupSizes = 'list',
        groupIDs = 'list',
        samplers = 'list',
        timing = 'list',
        ess = 'list',
        essPT = 'list',
        empCov = 'list',
        empCor = 'list',
        distMatrix = 'list',
        hTree = 'list'
        ),

    methods = list(

        initialize = function(Rmodel, control=list()) {
            require(lattice)
            require(coda)
            require(nimble)
            abModel <<- autoBlockModel(Rmodel)
            defaultsList <- autoBlockParamDefaults()
            for(i in seq_along(defaultsList)) if(is.null(control[[names(defaultsList)[i]]])) control[[names(defaultsList)[i]]] <- defaultsList[[i]]
            for(i in seq_along(control)) eval(substitute(verbose <<- VALUE, list(verbose=as.name(names(control)[i]), VALUE=control[[i]])))
            it <<- 0
        },

        run = function(runList) {
            if(!is.list(runList)) stop('runList argument should be a list')
            if(is.null(names(runList))) names(runList) <- rep('', length(runList))

            abModel$createInitialMCMCspec(runList)  ## here is where the initial MCMC spec is created, for re-use -- for new version
            
            for(i in seq_along(runList)) {
                runListElement <- runList[[i]]
                runListName <- names(runList)[i]
                if(is.character(runListElement)) {
                    type <- runListElement
                } else if(is.list(runListElement)) {
                    type <- 'blocks'
                } else if(class(runListElement) == '{') {
                    type <- 'spec'
                } else stop('don\'t understand element in run list')
                switch(type,
                       
                       none =    { specList <- list(createSpecFromGroups(abModel$nodeGroupScalars))
                                   runSpecListAndSaveBest(specList, 'none') },

                       all =     { specList <- list(createSpecFromGroups(abModel$nodeGroupAllBlocked))
                                   runSpecListAndSaveBest(specList, 'all') },

                       default = { specList <- list(configureMCMC(oldSpec = abModel$initialMCMCspec))
                                   runSpecListAndSaveBest(specList, 'default') },

                       blocks =  { specList <- list(createSpecFromGroups(abModel$createGroups(runListElement)))
                                   name <- if(runListName == '') 'customBlocks' else runListName
                                   runSpecListAndSaveBest(specList, name) },

                       spec =    { Rmodel <- abModel$Rmodel  ## just hoping that the customSpec will find this
                                   specList <- list(eval(runListElement, envir=environment()))
                                   name <- if(runListName == '') 'customSpec' else runListName
                                   runSpecListAndSaveBest(specList, name) },

                       auto =    { autoIt <- 0
                                   while((autoIt < 2) || ((!groupingsEquiv(grouping[[it]], grouping[[it-1]])) && (min(essPT[[it]]) > min(essPT[[it-1]])))) {
                                       candidateGroupsList <- if(autoIt==0) list(abModel$nodeGroupScalars)  else determineCandidateGroupsFromCurrentSample()
                                       specList <- lapply(candidateGroupsList, function(groups) createSpecFromGroups(groups))
                                       runSpecListAndSaveBest(specList, paste0('auto',autoIt), auto=TRUE)
                                       autoIt <- autoIt + 1
                                   }
                               },
                       stop('don\'t understand element in run list'))
            }

            names(candidateGroups) <<- naming
            names(grouping) <<- naming
            names(groupSizes) <<- naming
            names(groupIDs) <<- naming
            names(samplers) <<- naming
            names(timing) <<- naming
            names(ess) <<- naming
            names(essPT) <<- naming
        },

        determineCandidateGroupsFromCurrentSample = function() {
            cutree_heights <- seq(0, 1, by=0.1)
            cutreeList <- lapply(cutree_heights, function(height) cutree(hTree[[it]], h = height))
            names(cutreeList) <- paste0('cut', cutree_heights)
            uniqueCutreeList <- unique(cutreeList)
            for(i in seq_along(uniqueCutreeList)) { for(j in seq_along(cutreeList)) { if(all(uniqueCutreeList[[i]]==cutreeList[[j]])) { names(uniqueCutreeList)[i] <- names(cutreeList)[j]; break } } }
            candidateGroupsList <- lapply(uniqueCutreeList, function(ct) determineGroupsFromCutree(ct))
            return(candidateGroupsList)
        },
        
        determineGroupsFromCutree = function(ct) {
            return(lapply(unique(ct), function(x) names(ct)[ct==x]))
        },
        
        runSpecListAndSaveBest = function(specList, name, auto=FALSE) {
            lapply(specList, function(spec) checkOverMCMCspec(spec))
            RmcmcList <- lapply(specList, function(spec) buildMCMC(spec))
            CmcmcList <- compileNimble(RmcmcList, project = abModel$Rmodel)
            if(!is.list(CmcmcList)) CmcmcList <- list(CmcmcList)  ## make sure compileNimble() returns a list...
            timingList <- essList <- essPTList <- essPTminList <- list()
            for(i in seq_along(CmcmcList)) {
                timingList[[i]] <- as.numeric(system.time(CmcmcList[[i]]$run(niter))[3])
                burnedSamples <- extractAndBurnSamples(CmcmcList[[i]])
                essList[[i]]   <- apply(burnedSamples, 2, effectiveSize)
                essPTList[[i]] <- essList[[i]] / timingList[[i]]
                essPTminList[[i]] <- sort(essPTList[[i]])[1]
            }
            bestInd <- as.numeric(which(unlist(essPTminList) == max(unlist(essPTminList))))
            if(!is.null(names(specList))) name <- paste0(name, '-', names(specList)[bestInd])
            
            it <<- it + 1
            naming[[it]] <<- name
            candidateGroups[[it]] <<- lapply(specList, function(spec) determineGroupsFromSpec(spec))
            grouping[[it]] <<- candidateGroups[[it]][[bestInd]]
            groupSizes[[it]] <<- determineNodeGroupSizesFromGroups(grouping[[it]])
            groupIDs[[it]] <<- determineNodeGroupIDsFromGroups(grouping[[it]])
            samplers[[it]] <<- determineSamplersFromGroupsAndSpec(grouping[[it]], specList[[bestInd]])
            timing[[it]] <<- timingList[[bestInd]]
            ess[[it]] <<- essList[[bestInd]]
            essPT[[it]] <<- sort(essPTList[[bestInd]])
            
            if(auto) {
                burnedSamples <- extractAndBurnSamples(CmcmcList[[bestInd]])
                empCov[[it]] <<- cov(burnedSamples)
                empCor[[it]] <<- cov2cor(empCov[[it]])
                distMatrix[[it]] <<- as.dist(1 - abs(empCor[[it]]))
                hTree[[it]] <<- hclust(distMatrix[[it]])
            }
            
            if(verbose) printCurrent(name, specList[[bestInd]])
            if(makePlots && auto) makeCurrentPlots(name)
        },

        extractAndBurnSamples = function(Cmcmc) {
            samples <- as.matrix(Cmcmc$mvSamples)
            ## make sure we don't keep samples from deterministic nodes
            namesToKeep <- setdiff(dimnames(samples)[[2]], abModel$Rmodel$getNodeNames(determOnly=TRUE, returnScalarComponents=TRUE))
            burnedSamples <- samples[(floor(niter/2)+1):niter, namesToKeep]
            burnedSamples
        },

        determineGroupsFromSpec = function(spec) {
            groups <- list()
            for(ss in spec$samplerSpecs) {
                if(ss$type == 'RW_block') {
                    nodes <- ss$control$targetNodes
                } else if(ss$type == 'crossLevel') {
                    topNodes <- ss$control$topNodes
                    lowNodes <- spec$model$getDependencies(topNodes, self=FALSE, stochOnly=TRUE, includeData=FALSE)
                    nodes <- c(topNodes, lowNodes)
                } else if(!is.null(ss$control$targetNode)) {
                    nodes <- ss$control$targetNode
                } else stop('don\'t understand sampler type')
                groups[[length(groups)+1]] <- spec$model$expandNodeNames(nodes, returnScalarComponents=TRUE)
            }
            return(groups)
        },

        determineNodeGroupSizesFromGroups = function(groups) {
            groupSizeVector <- numeric(0)
            for(gp in groups) for(node in gp) groupSizeVector[[node]] <- length(gp)
            return(groupSizeVector)
        },

        determineNodeGroupIDsFromGroups = function(groups) {
            groupIDvector <- numeric(0)
            for(i in seq_along(groups)) for(node in groups[[i]]) groupIDvector[[node]] <- i
            return(groupIDvector)
        },

        determineSamplersFromGroupsAndSpec = function(groups, spec) {
            samplerSpecs <- spec$samplerSpecs
            if(length(groups) != length(samplerSpecs)) stop('something wrong')
            samplerVector <- character(0)
            for(i in seq_along(groups)) for(node in groups[[i]]) samplerVector[[node]] <- samplerSpecs[[i]]$type
            return(samplerVector)
        },
        
        createSpecFromGroups = function(groups) {
            ##spec <- configureMCMC(Rmodel, nodes=NULL, monitors=character(0)) ## original version
            spec <- configureMCMC(oldSpec = abModel$initialMCMCspec)  ## new version
            spec$setSamplers()  ## new version -- removes all the samplers from initalMCMCspec
            for(nodeGroup in groups) addSamplerToSpec(abModel$Rmodel, spec, nodeGroup)
            return(spec)
        },
        
        checkOverMCMCspec = function(spec) {
            warn <- FALSE
            for(ss in spec$samplerSpecs) {
                ## if(ss$type == 'end') {
                ##     msg <- 'using \'end\' sampler may lead to results we don\'t want'
                ##     cat(paste0('\nWARNING: ', msg, '\n\n')); warning(msg)
                ## }
                if(grepl('^conjugate_', ss$type) && nimble:::nimbleOptions$verifyConjugatePosterior) {
                    ##msg <- 'conjugate sampler running slow due to checking the posterior'
                    ##cat(paste0('\nWARNING: ', msg, '\n\n')); warning(msg)
                    warn <- TRUE
                }
            }
            if(warn) {
                msg <- 'Conjugate sampler functions in \'default\' spec are running slow due to verifying the posterior;\nThis behaviour can be changed using a NIMBLE package option.'
                warning(msg, call. = FALSE)
            }
        },

        printCurrent = function(name, spec) {
            cat(paste0('\n################################\nbegin iteration ', it, ': ', name, '\n################################\n'))
            if(length(candidateGroups[[it]]) > 1) { cat('\ncandidate groups:\n'); cg<-candidateGroups[[it]]; for(i in seq_along(cg)) { cat(paste0('\n',names(cg)[i],':\n')); printGrouping(cg[[i]]) } }
            cat('\ngroups:\n'); printGrouping(grouping[[it]])
            cat('\nsamplers:\n'); spec$getSamplers()
            cat(paste0('\nMCMC runtime: ', round(timing[[it]], 1), ' seconds\n'))
            cat('\nESS:\n'); print(round(ess[[it]], 0))
            cat('\nESS/time:\n'); print(round(essPT[[it]], 1))
            cat(paste0('\n################################\nend iteration ', it, ': ', name, '\n################################\n\n'))
        },

        makeCurrentPlots = function(name) {
            dev.new()
            if(inherits(try(plot(as.dendrogram(hTree[[it]]), ylim=c(0,1), main=name), silent=TRUE), 'try-error')) dev.off()
        },

        printGrouping = function(g) {
            for(i in seq_along(g)) cat(paste0('[', i, '] ', paste0(g[[i]], collapse=', '), '\n'))
        },

        groupingsEquiv = function(grouping1, grouping2) {
            grouping1 <- lapply(grouping1, sort)
            grouping2 <- lapply(grouping2, sort)
            while(length(grouping1) > 0) {
                grp1 <- grouping1[[1]]
                found <- FALSE
                for(i in seq_along(grouping2)) {
                    grp2 <- grouping2[[i]]
                    if(identical(grp1, grp2)) {
                        found <- TRUE
                        grouping1[1] <- grouping2[i] <- NULL
                        break
                    }
                }
                if(!found) return(FALSE)
            }
            if(length(grouping2) == 0) return(TRUE) else return(FALSE)
        }
        )
)

addSamplerToSpec <- function(Rmodel, spec, nodeGroup) {
    if(length(nodeGroup) > 1) {
        spec$addSampler(type = 'RW_block', control = list(targetNodes=nodeGroup), print = FALSE); return()
    }
    if(!(nodeGroup %in% Rmodel$getNodeNames())) {
        spec$addSampler(type = 'RW', control = list(targetNode=nodeGroup), print = FALSE); return()
    }
    if(nodeGroup %in% Rmodel$getMaps('nodeNamesEnd')) {
        ##cat(paste0('warning: using \'end\' sampler for node ', nodeGroup, ' may lead to results we don\'t want\n\n'))
        spec$addSampler(type = 'end', control = list(targetNode=nodeGroup), print = FALSE); return()
    }
    ## conjugacyResult <- Rmodel$checkConjugacy(nodeGroup)
    ## if((!is.null(conjugacyResult)) && conjOveride) {
    ##     spec$addSampler(type = conjugacyResult$samplerType, control = conjugacyResult$control, print = FALSE); return()
    ## }
    if(Rmodel$getNodeInfo()[[nodeGroup]]$isDiscrete()) {
        spec$addSampler(type = 'slice', control = list(targetNode=nodeGroup), print = FALSE); return()
    }
    if(length(Rmodel$expandNodeNames(nodeGroup, returnScalarComponents = TRUE)) > 1) {
        spec$addSampler(type = 'RW_block', control = list(targetNodes=nodeGroup), print = FALSE); return()
    }
    spec$addSampler(type = 'RW', control = list(targetNode=nodeGroup), print = FALSE); return()
}

createDFfromABlist <- function(lst, niter) {
    df <- data.frame(model=character(), blocking=character(), timing=numeric(), node=character(), groupSize = numeric(), groupID = numeric(), sampler = character(), ess=numeric(), essPT=numeric(), stringsAsFactors=FALSE)
    for(iAB in seq_along(lst)) {
        ab <- lst[[iAB]]
        abName <- names(lst)[iAB]
        for(iBlock in seq_along(ab$naming)) {
            blocking <- ab$naming[[iBlock]]
            timing <- ab$timing[[iBlock]]
            ess <- ab$ess[[iBlock]]
            nodes <- names(ess)
            essPT <- ab$essPT[[iBlock]][nodes]            ## sort
            groupSizes <- ab$groupSizes[[iBlock]][nodes]  ##
            groupIDs <- ab$groupIDs[[iBlock]][nodes]      ##
            samplers <- ab$samplers[[iBlock]][nodes]      ##
            newIndDF <- (1:length(nodes)) + dim(df)[1]
            df[newIndDF,] <- NA
            df[newIndDF,]$model <- abName
            df[newIndDF,]$blocking <- blocking
            df[newIndDF,]$timing <- timing
            df[newIndDF,]$node <- nodes
            df[newIndDF,]$groupSize <- groupSizes
            df[newIndDF,]$groupID <- groupIDs
            df[newIndDF,]$sampler <- samplers
            df[newIndDF,]$ess <- ess
            df[newIndDF,]$essPT <- essPT
        }
    }
    df$timePer10k <- df$timing * 10000/niter
    df$essPer10k  <- df$ess    * 10000/niter * 2
    df$Efficiency <- df$essPer10k / df$timePer10k
    df$mcmc <- gsub('-.+', '', df$blocking)
    return(df)
}



plotABS <- function(df, xlimToMin=FALSE, together) {
    models <- unique(df$model)
    nModels <- length(models)
    if(missing(together)) together <- if(nModels <= 5) TRUE else FALSE
    nVertPlots <- if(together) nModels*2 else nModels
    xVarNames <- c('ess', 'essPT')
    parCmd <- quote(par(mfrow=c(nVertPlots,1),mar=c(1,0,1,0),tcl=-.1,mgp=c(3,0,0),cex.axis=.7))
    if(together) { quartz(); eval(parCmd) }
    for(xVarName in xVarNames) {
        if(!together) { quartz(); eval(parCmd) }
        maxMinXVar<-0; for(mod in models) {dfMod<-df[df$model==mod,]; blks<-unique(dfMod$blocking); for(blk in blks) {maxMinXVar<-max(maxMinXVar,min(dfMod[dfMod$blocking==blk,xVarName]))}}
        maxXVar <- if(xlimToMin) maxMinXVar else max(df[, xVarName])
        xlim <- c(maxXVar*-0.05, maxXVar)
        maxTiming <- max(df[, 'timing'])
        for(mod in models) {
            dfMod <- df[df$model==mod,]
            blockings <- unique(dfMod$blocking)
            nBlockings <- length(blockings)
            bestBlk<-''; bestEssPT<-0; for(blk in blockings) { if(min(dfMod[dfMod$blocking==blk,'essPT'])>bestEssPT) {bestEssPT<-min(dfMod[dfMod$blocking==blk,'essPT']); bestBlk<-blk} }
            plot(-100,-100,xlim=xlim,ylim=c(0,nBlockings+1),xlab='',ylab='',main=paste0(xVarName, ' for model ', mod))
            for(iBlocking in 1:nBlockings) {
                blocking <- blockings[iBlocking]
                dfModBlock <- dfMod[dfMod$blocking==blocking,]
                xVarValues <- dfModBlock[,xVarName]
                groupSizes <- dfModBlock[,'groupSize']
                timing <- dfModBlock[,'timing'][1]   # only first element
                timingOnXaxis <- timing/maxTiming * xlim[2]
                yCoord <- nBlockings+1-iBlocking
                lines(x=c(0,timingOnXaxis), y=rep(yCoord,2), lty=1, lwd=2, col='lightgrey')
                col <- if(blocking == bestBlk) 'green' else 'black'
                text(x=xVarValues, y=yCoord, labels=groupSizes, cex=0.7, col=col)
                col <- if(blocking == bestBlk) 'green' else 'blue'
                text(x=xlim[1], y=yCoord, labels=blocking, col=col)
                if(timing==maxTiming) text(xlim[2], yCoord+1, paste0('t = ',round(timing,1)))
            }
        }
    }
}


printMinTimeABS <- function(df, round=TRUE, addAutoMax=TRUE, sortOutput=FALSE) {
    namesToRemove <- intersect(c('groupID', 'sampler'), names(df))
    for(name in namesToRemove) { ind <- which(names(df)==name); df <- df[, -ind] }
    models <- unique(df$model)
    cat('\n')
    dfReturn <- data.frame()
    for(mod in models) {
        dfMod <- df[df$model == mod, ]
        blockings <- unique(dfMod$blocking)
        dfOut <- dfMod[numeric(0), ]
        for(blk in blockings) {
            dfModBlk <- dfMod[dfMod$blocking == blk, ]
            ind <- which(dfModBlk$essPT == min(dfModBlk$essPT))[1]
            dfOut[dim(dfOut)[1] + 1, ] <- dfModBlk[ind, ]
        }
        if(sortOutput) dfOut <- dfOut[sort(dfOut$essPT,index.return=TRUE)$ix, ]
        dimnames(dfOut)[[1]] <- 1:(dim(dfOut)[1])
        if(round) {
            dfOut$timing     <- round(dfOut$timing, 2)
            dfOut$timePer10k <- round(dfOut$timePer10k, 2)
            dfOut$ess        <- round(dfOut$ess, 1)
            dfOut$essPer10k  <- round(dfOut$essPer10k, 1)
            dfOut$essPT      <- round(dfOut$essPT, 1)
            dfOut$Efficiency <- round(dfOut$Efficiency, 1)
        }
        if(addAutoMax && ('auto0' %in% blockings)) {
            autoBlockings <- blockings[grepl('^auto', blockings)]
            dfAuto <- dfOut[dfOut$blocking %in% autoBlockings,]
            maxEffInd <- which(dfAuto$Efficiency == max(dfAuto$Efficiency))
            nextInd <- dim(dfOut)[1] + 1
            dfOut[nextInd,] <- dfAuto[maxEffInd,]
            dfOut[nextInd, 'blocking'] <- dfOut[nextInd, 'mcmc'] <- 'autoMax'
        }
        print(dfOut)
        cat('\n')
        dfReturn <- rbind(dfReturn, dfOut)
    }
    return(invisible(dfReturn))
}


reduceDF <- function(df, addAutoMax=TRUE, sortOutput=TRUE, round=TRUE) {
    df = data.frame(mcmc=df$mcmc, node=df$node, S=df$essPer10k, C=df$timePer10k, Efficiency=df$Efficiency, stringsAsFactors=FALSE)
    dfOut <- df[numeric(), ]
    mcmcs <- unique(df$mcmc)
    for(mcmc in mcmcs) {
        dfBlk <- df[df$mcmc == mcmc, ]
        ind <- which(dfBlk$Efficiency == min(dfBlk$Efficiency))[1]
        dfOut[dim(dfOut)[1]+1, ] <- dfBlk[ind, ]
    }
    dfOut[dfOut$mcmc=='auto0', 'mcmc'] <- 'All Scalar'
    dfOut[dfOut$mcmc=='all', 'mcmc'] <- 'All Blocked'
    dfOut[dfOut$mcmc=='default', 'mcmc'] <- 'Default'
    if(addAutoMax) {
        autoBlockings <- dfOut$mcmc[grepl('^auto', dfOut$mcmc)]
        autoLast <- autoBlockings[length(autoBlockings)]
        ## replace autoLast with 'autoMax'
        dfOut[dfOut$mcmc==autoLast, 'mcmc'] <- 'Auto-Blocking'
        ## remove any remaining 'auto#' entries
        dfOut <- dfOut[!dfOut$mcmc %in% autoBlockings,]
    }
    if(sortOutput) dfOut <- dfOut[sort(dfOut$Efficiency,index.return=TRUE)$ix, ]
    dimnames(dfOut)[[1]] <- 1:(dim(dfOut)[1])
    if(round) {
        dfOut$S          <- round(dfOut$S, 2)
        dfOut$C          <- round(dfOut$C, 2)
        dfOut$Efficiency <- round(dfOut$Efficiency, 2)
    }
    return(dfOut)
}






