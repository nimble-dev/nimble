#' Automated parameter blocking procedure for efficient MCMC sampling
#' 
#' Runs NIMBLE's automated blocking procedure for a given model object, to dynamically determine a blocking scheme of the continuous-valued model nodes.  This blocking scheme is designed to produce efficient MCMC sampling (defined as number of effective samples generated per second of algorithm runtime).  See Turek, et al (2015) for details of this algorithm.  This also (optionally) compares this blocked MCMC against several static MCMC algorithms, including all univariate sampling, blocking of all continuous-valued nodes, NIMBLE's default MCMC configuration, and custom-specified blockings of parameters.
#' 
#' This method allows for fine-tuned usage of the automated blocking procedure.  However, the main entry point to the automatic blocking procedure is intended to be through either buildMCMC(..., autoBlock = TRUE), or configureMCMC(..., autoBlock = TRUE).
#' 
#' @author Daniel Turek
#'
#' @seealso configureMCMC buildMCMC
#'
#' @param Rmodel A NIMBLE model object, created from \code{\link{nimbleModel}}.
#'
#' @param autoIt The number of MCMC iterations to run intermediate MCMC algorithms, through the course of the procedure.  Default 20,000.
#'
#' @param run List of additional MCMC algorithms to compare against the automated blocking MCMC.  These may be specified as: the character string 'all' to denote blocking all continuous-valued nodes; the character string 'default' to denote NIMBLE's default MCMC configuration; a named list element consisting of a quoted code block, which when executed returns an MCMC configuration object for comparison; a custom-specificed blocking scheme, specified as a named list element which itself is a list of character vectors, where each character vector specifies the nodes in a particular block.  Default is c('all', 'default').
#'
#' @param verbose Logical specifying whether to output considerable details of the automated block procedure, through the course of execution.  Default FALSE.
#' 
#' @param setSeed Logical specificying whether to call set.seed(0) prior to beginning the blocking procedure.  Default TRUE.
#'
#' @param makePlots Logical specifying whether to plot the hierarchical clustering dendrograms, through the course of execution.  Default FALSE.
#'
#' @param round Logical specifying whether to round the final output results to two decimal places.  Default TRUE.
#' 
#' @return Returns a named list containing elements:
#' \itemize{
#' \item \code{summary}: A data frame containing a numerical summary of the performance of all MCMC algorithms (including that from automated blocking)
#' \item \code{autoGroups}: A list specifying the parameter blockings converged on by the automated blocking procedure
#' \item \code{conf}: A NIMBLE MCMC configuration object corresponding to the results of the automated blocking procedure
#' }
#' 
#' @references
#'
#' Turek, D., de Valpine, P., Paciorek, C., and Anderson-Bergman, C. (2015). Automated Parameter Blocking for Efficient Markov-Chain Monte Carlo Sampling. <arXiv:1503.05621>. 
#'
#' @export
autoBlock <- function(Rmodel,
                      autoIt = 20000,
                      run = list('all', 'default'),
                      setSeed = TRUE,
                      verbose = FALSE,
                      makePlots = FALSE,
                      round = TRUE ) {
    if(autoIt < 10000) stop('Minimum auto-blocking iterations is 10,000')
    control <- list(niter=autoIt, setSeed=setSeed, verbose=verbose, makePlots=makePlots)
    ab <- autoBlockClass(Rmodel, control)
    if(!'auto' %in% run) run <- c(run, 'auto')  ## always use 'autoBlock' routine
    ab$run(run)
    abList <- list(ab)
    names(abList)[1] <- 'model'
    df <- createDFfromABlist(abList, autoIt)
    dfmin <- reduceDF(df, round = round)
    cat('\nAuto-Blocking summary:\n')
    print(dfmin)
    lastAutoInd <- max(grep('^auto', ab$naming))   ## index of final 'auto' iteration
    lastAutoGrouping <- ab$grouping[[lastAutoInd]]  ## grouping of final 'auto' iteration
    nonTrivialGroups <- lastAutoGrouping[unlist(lapply(lastAutoGrouping, function(x) length(x)>1))]
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
## create a new MCMC conf with the autoBlock groupings:
    conf <- configureMCMC(Rmodel, nodes = NULL)
    for(nodeGroup in lastAutoGrouping) addSamplerToConf(Rmodel, conf, nodeGroup)
    retList <- list(summary=dfmin, autoGroups=nonTrivialGroups, conf=conf)
    return(invisible(retList))
}



autoBlockModel <- setRefClass(
    Class = 'autoBlockModel',
    fields = list(
        Rmodel_orig = 'ANY',
        Rmodel = 'ANY',
        Cmodel = 'ANY',
        md = 'ANY',
        scalarNodeVector = 'character',
        scalarNodeVectorCont = 'character',
        scalarNodeVectorDisc = 'character',
        nodeGroupScalars = 'list',
        nodeGroupAllBlocked = 'list',
        monitorsVector = 'character',
        initialMCMCconf = 'ANY'
    ),
    methods = list(
        initialize = function(Rmodel_orig) {
            Rmodel_orig <<- Rmodel_orig
            md <<- Rmodel_orig$modelDef
            Rmodel <<- Rmodel_orig$newModel(replicate = TRUE, check = FALSE)
            ##nimCopy(from = Rmodel_orig, to = Rmodel, logProb = TRUE)
            ##for(var in ls(Rmodel_orig$isDataEnv)) Rmodel$isDataEnv[[var]] <<- Rmodel_orig$isDataEnv[[var]]  ## copies data flags to the new model
            scalarNodeVector <<- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE, returnScalarComponents=TRUE)
            discreteInd <- sapply(scalarNodeVector, function(n) Rmodel$isDiscrete(n), USE.NAMES=FALSE)
            scalarNodeVectorCont <<- scalarNodeVector[!discreteInd]   ## making work with discrete nodes
            scalarNodeVectorDisc <<- scalarNodeVector[ discreteInd]   ## making work with discrete nodes
            if(length(scalarNodeVectorCont) == 0) stop('autoBlocking only works with one or more continuous-valued model nodes')   ## making work with discrete nodes
            nodeGroupScalars <<- c(unique(lapply(scalarNodeVectorDisc, Rmodel$expandNodeNames)), scalarNodeVectorCont)   ## making work with discrete nodes, and also with dmulti distributions
            ##nodeGroupAllBlocked <<- list(scalarNodeVector)   ## making work with discrete nodes
            ##nodeGroupAllBlocked <<- c(lapply(scalarNodeVectorDisc, function(x) x), list(scalarNodeVectorCont))   ## making work with discrete nodes
            nodeGroupAllBlocked <<- c(unique(lapply(scalarNodeVectorDisc, Rmodel$expandNodeNames)), list(scalarNodeVectorCont))   ## making work with discrete nodes, and also with dmulti distributions
            monitorsVector <<- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE)
        },
        ## here is where the initial MCMC conf is created, for re-use -- for new version
        createInitialMCMCconf = function(runList) {
            initialMCMCconf <<- configureMCMC(Rmodel)
            nInitialSamplers <- length(initialMCMCconf$samplerConfs)
            initialMCMCconf$addSampler(target = scalarNodeVectorCont[1], type = 'slice',    print=FALSE)  ## add one slice sampler
            initialMCMCconf$addSampler(target = scalarNodeVectorCont[1], type = 'RW',       print=FALSE)  ## add one RW sampler
            initialMCMCconf$addSampler(target = scalarNodeVectorCont[1], type = 'RW_block', print=FALSE, silent = TRUE)  ## add one RW_block sampler
            addCustomizedSamplersToInitialMCMCconf(runList)
            initialMCMCconf$addMonitors(monitorsVector, print=FALSE)
            RinitialMCMC <- buildMCMC(initialMCMCconf)
            Cmodel <<- compileNimble(Rmodel)
            CinitialMCMC <- compileNimble(RinitialMCMC, project = Rmodel)   ## (new version) yes, we need this compileNimble call -- this is the whole point!
            initialMCMCconf$setSamplers(1:nInitialSamplers, print=FALSE)  ## important for new version: removes all news samplers added to initial MCMC conf
        },
        addCustomizedSamplersToInitialMCMCconf = function(runListCode) {
            if(is.list(runListCode)) { lapply(runListCode, function(el) addCustomizedSamplersToInitialMCMCconf(el)); return() }
            if(is.call(runListCode)) {
                if(is.call(runListCode[[1]]) && length(runListCode[[1]])==3 && runListCode[[1]][[3]]=='addSampler') {
                    runListCode[[1]][[2]] <- as.name('initialMCMCconf')
                    eval(substitute(RUNLISTCODE, list(RUNLISTCODE=runListCode)))
                    return()
                }
                lapply(runListCode, function(el) addCustomizedSamplersToInitialMCMCconf(el))
                return()
            }
        },
        createGroups = function(listOfBlocks = list()) {
            listOfBlocks <- lapply(listOfBlocks, function(blk) Rmodel$expandNodeNames(blk, returnScalarComponents=TRUE))
            if(any(unlist(listOfBlocks) %in% scalarNodeVectorDisc)) stop('cannot put block sampler on discrete-valued model nodes')
            nodes <- scalarNodeVector
            nodes <- setdiff(nodes, unlist(listOfBlocks))
            nodeList <- lapply(nodes, function(x) x)
            for(ng in listOfBlocks) nodeList[[length(nodeList)+1]] <- ng
            return(nodeList)
        },
        resetCmodelInitialValues = function() {
            nimCopy(from = Rmodel_orig, to = Cmodel, logProb = TRUE)
            calculate(Cmodel)
        }
    )
)



autoBlockParamDefaults <- function() {
	list(
            makePlots = FALSE,
            niter = 20000,
            setSeed = TRUE,
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
        setSeed = 'logical',
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
            abModel <<- autoBlockModel(Rmodel)
            defaultsList <- autoBlockParamDefaults()
            for(i in seq_along(defaultsList)) if(is.null(control[[names(defaultsList)[i]]])) control[[names(defaultsList)[i]]] <- defaultsList[[i]]
            for(i in seq_along(control)) eval(substitute(verbose <<- VALUE, list(verbose=as.name(names(control)[i]), VALUE=control[[i]])))
            it <<- 0
        },

        run = function(runList) {
            if(!is.list(runList)) stop('runList argument should be a list')
            if(is.null(names(runList))) names(runList) <- rep('', length(runList))

            abModel$createInitialMCMCconf(runList)  ## here is where the initial MCMC conf is created, for re-use -- for new version
            
            for(i in seq_along(runList)) {
                runListElement <- runList[[i]]
                runListName <- names(runList)[i]
                if(is.character(runListElement)) {
                    type <- runListElement
                } else if(is.list(runListElement)) {
                    type <- 'blocks'
                } else if(class(runListElement) == '{') {
                    type <- 'conf'
                } else stop('don\'t understand element in run list')
                switch(type,
                       
                       none =    { confList <- list(createConfFromGroups(abModel$nodeGroupScalars))
                                   runConfListAndSaveBest(confList, 'none') },

                       all =     { confList <- list(createConfFromGroups(abModel$nodeGroupAllBlocked))
                                   runConfListAndSaveBest(confList, 'all') },

                       default = { ##confList <- list(configureMCMC(oldConf = abModel$initialMCMCconf))
                                   ## forcing this processing through createConfFromGroups()
                                   ## in order to always standardize the ordering of samplers;
                                   ## even though this might result in a different sampler ordering
                                   ## than the true NIMBLE 'default' MCMC conf
                                   ##groups <- determineGroupsFromConf(abModel$initialMCMCconf)
                                   groups <- lapply(determineGroupsFromConf(abModel$initialMCMCconf), function(nodes) unique(abModel$Rmodel$expandNodeNames(nodes)))  ## making work with dmulti distribution
                                   confList <- list(createConfFromGroups(groups))
                                   runConfListAndSaveBest(confList, 'default') },

                       blocks =  { confList <- list(createConfFromGroups(abModel$createGroups(runListElement)))
                                   name <- if(runListName == '') 'customBlocks' else runListName
                                   runConfListAndSaveBest(confList, name) },

                       conf =    { Rmodel <- abModel$Rmodel  ## just hoping that the customConf will find this
                                   confList <- list(eval(runListElement, envir=environment()))
                                   name <- if(runListName == '') 'customConf' else runListName
                                   runConfListAndSaveBest(confList, name) },

                       auto =    { autoIt <- 0
                                   while((autoIt < 2) || ((!groupingsEquiv(grouping[[it]], grouping[[it-1]])) && (min(essPT[[it]]) > min(essPT[[it-1]])))) {
                                       candidateGroupsList <- if(autoIt==0) list(abModel$nodeGroupScalars)  else determineCandidateGroupsFromCurrentSample()
                                       confList <- lapply(candidateGroupsList, function(groups) createConfFromGroups(groups))
                                       runConfListAndSaveBest(confList, paste0('auto',autoIt), auto=TRUE)
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
            groupsContOnly <- lapply(unique(ct), function(x) names(ct)[ct==x])   ## making work with discrete nodes
            ##groupsAllNodes <- c(lapply(abModel$scalarNodeVectorDisc, function(x) x), groupsContOnly)   ## making work with discrete nodes
            groupsAllNodes <- c(unique(lapply(abModel$scalarNodeVectorDisc, abModel$Rmodel$expandNodeNames)), groupsContOnly)   ## making work with discrete nodes and dmulti distribution
            return(groupsAllNodes)   ## making work with discrete nodes
        },
        
        runConfListAndSaveBest = function(confList, name, auto=FALSE) {
            lapply(confList, function(conf) checkOverMCMCconf(conf))
            RmcmcList <- lapply(confList, function(conf) buildMCMC(conf))
            CmcmcList <- compileNimble(RmcmcList, project = abModel$Rmodel)
            if(!is.list(CmcmcList)) CmcmcList <- list(CmcmcList)  ## make sure compileNimble() returns a list...
            timingList <- essList <- essPTList <- essPTminList <- list()
            for(i in seq_along(CmcmcList)) {
                if(setSeed) set.seed(0)
                abModel$resetCmodelInitialValues()
                timingList[[i]] <- as.numeric(system.time(CmcmcList[[i]]$run(niter, progressBar = FALSE))[3])
                burnedSamples <- extractAndBurnSamples(CmcmcList[[i]])
                essList[[i]] <- apply(burnedSamples, 2, effectiveSize)
                essList[[i]] <- essList[[i]][essList[[i]] > 0]  ## exclude nodes with ESS=0 -- for discrete nodes which are fixed to a certain value; making work with discrete nodes
                essPTList[[i]] <- essList[[i]] / timingList[[i]]
                essPTminList[[i]] <- sort(essPTList[[i]])[1]
            }
            bestInd <- as.numeric(which(unlist(essPTminList) == max(unlist(essPTminList))))
            if(length(bestInd) > 1) stop('there should never be an exact tie for the best...')
            if(!is.null(names(confList))) name <- paste0(name, '-', names(confList)[bestInd])
            
            it <<- it + 1
            naming[[it]] <<- name
            candidateGroups[[it]] <<- lapply(confList, function(conf) determineGroupsFromConf(conf))
            grouping[[it]] <<- candidateGroups[[it]][[bestInd]]
            groupSizes[[it]] <<- determineNodeGroupSizesFromGroups(grouping[[it]])
            groupIDs[[it]] <<- determineNodeGroupIDsFromGroups(grouping[[it]])
            samplers[[it]] <<- determineSamplersFromGroupsAndConf(grouping[[it]], confList[[bestInd]])
            timing[[it]] <<- timingList[[bestInd]]
            ess[[it]] <<- essList[[bestInd]]
            essPT[[it]] <<- sort(essPTList[[bestInd]])
            
            if(auto) {
                burnedSamples <- extractAndBurnSamples(CmcmcList[[bestInd]])
                burnedSamples <- burnedSamples[, abModel$scalarNodeVectorCont]   ## making work with discrete nodes

                ##empCov[[it]] <<- cov(burnedSamples)
                e <- try(cov(burnedSamples))
                if(inherits(e, 'try-error')) {
                    message('try-error, going into browser'); browser(); 1; 2
                } else empCov[[it]] <<- e

                ##empCor[[it]] <<- cov2cor(empCov[[it]])
                e <- try(cov2cor(empCov[[it]]))
                if(inherits(e, 'try-error')) {
                    message('try-error, going into browser'); browser(); 3; 4
                } else empCor[[it]] <<- e

                distMatrix[[it]] <<- as.dist(1 - abs(empCor[[it]]))

                ##hTree[[it]] <<- hclust(distMatrix[[it]])
                e <- try(hclust(distMatrix[[it]]))
                if(inherits(e, 'try-error')) {
                    message('try-error, going into browser'); browser(); 5; 6
                } else hTree[[it]] <<- e
            }
            
            if(verbose) printCurrent(name, confList[[bestInd]])
            if(makePlots && auto) makeCurrentPlots(name)
        },

        extractAndBurnSamples = function(Cmcmc) {
            samples <- as.matrix(Cmcmc$mvSamples)
            ## make sure we don't keep samples from deterministic nodes
            namesToKeep <- setdiff(dimnames(samples)[[2]], abModel$Rmodel$getNodeNames(determOnly=TRUE, returnScalarComponents=TRUE))
            burnedSamples <- samples[(floor(niter/2)+1):niter, namesToKeep]
            burnedSamples
        },

        determineGroupsFromConf = function(conf) {
            groups <- list()
            for(ss in conf$samplerConfs) {
                if(ss$name == 'crossLevel') {
                    topNodes <- ss$target
                    lowNodes <- conf$model$getDependencies(topNodes, self=FALSE, stochOnly=TRUE, includeData=FALSE)
                    nodes <- c(topNodes, lowNodes)
                } else {
                    nodes <- ss$target
                }
                groups[[length(groups)+1]] <- conf$model$expandNodeNames(nodes, returnScalarComponents=TRUE)
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

        determineSamplersFromGroupsAndConf = function(groups, conf) {
            samplerConfs <- conf$samplerConfs
            if(length(groups) != length(samplerConfs)) stop('something wrong')
            samplerVector <- character(0)
            for(i in seq_along(groups)) for(node in groups[[i]]) samplerVector[[node]] <- samplerConfs[[i]]$name
            return(samplerVector)
        },
        
        createConfFromGroups = function(groups) {
            groups <- sortGroups(groups)
            ##conf <- configureMCMC(Rmodel, nodes=NULL, monitors=character(0)) ## original version
            conf <- configureMCMC(oldConf = abModel$initialMCMCconf)  ## new version
            conf$setSamplers()  ## new version -- removes all the samplers from initalMCMCconf
            for(nodeGroup in groups) addSamplerToConf(abModel$Rmodel, conf, nodeGroup)
            return(conf)
        },
        
        sortGroups = function(groups) {
            eachGroupSorted <- lapply(groups, sort)
            groupsAsStrings <- lapply(eachGroupSorted, function(grp) paste0(grp, collapse = '_'))
            sortedInd <- sort(unlist(groupsAsStrings), index.return = TRUE)$ix
            sortedGroups <- eachGroupSorted[sortedInd]
            return(sortedGroups)
        },
        
        checkOverMCMCconf = function(conf) {
            warn <- FALSE
            for(ss in conf$samplerConfs) {
                ## if(ss$name == 'posterior_predictive') {
                ##     msg <- 'using \'posterior_predictive\' sampler may lead to results we don\'t want'
                ##     cat(paste0('\nWARNING: ', msg, '\n\n')); warning(msg)
                ## }
                if(grepl('^conjugate_', ss$name) && getNimbleOption('verifyConjugatePosteriors')) {
                    ##msg <- 'conjugate sampler running slow due to checking the posterior'
                    ##cat(paste0('\nWARNING: ', msg, '\n\n')); warning(msg)
                    warn <- TRUE
                }
            }
            if(warn) {
                msg <- 'Conjugate sampler functions in \'default\' conf are running slow due to verifying the posterior;\nThis behaviour can be changed using a NIMBLE package option.'
                warning(msg, call. = FALSE)
            }
        },

        printCurrent = function(name, conf) {
            cat(paste0('\n################################\nbegin iteration ', it, ': ', name, '\n################################\n'))
            if(length(candidateGroups[[it]]) > 1) { cat('\ncandidate groups:\n'); cg<-candidateGroups[[it]]; for(i in seq_along(cg)) { cat(paste0('\n',names(cg)[i],':\n')); printGrouping(cg[[i]]) } }
            cat('\ngroups:\n'); printGrouping(grouping[[it]])
            cat('\nsamplers:\n'); conf$getSamplers()
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

addSamplerToConf <- function(Rmodel, conf, nodeGroup) {
    if(length(nodeGroup) > 1) {
        conf$addSampler(target = nodeGroup, type = 'RW_block', print = FALSE, silent = TRUE); return()
    }
    if(!(nodeGroup %in% Rmodel$getNodeNames()) && !Rmodel$isDiscrete(nodeGroup)) {
        conf$addSampler(target = nodeGroup, type = 'RW', print = FALSE); return()
    }
    if(nodeGroup %in% Rmodel$getMaps('nodeNamesEnd')) {
        ##cat(paste0('warning: using \'posterior_predictive\' sampler for node ', nodeGroup, ' may lead to results we don\'t want\n\n'))
        conf$addSampler(target = nodeGroup, type = 'posterior_predictive', print = FALSE); return()
    }
    ## conjugacyResult <- Rmodel$checkConjugacy(nodeGroup)
    ## if((!is.null(conjugacyResult)) && conjOveride) {
    ##     conf$addSampler(target = ??????, type = conjugacyResult$samplerType, control = conjugacyResult$control, print = FALSE); return()
    ## }
    if(Rmodel$isBinary(nodeGroup)) {
        conf$addSampler(target = nodeGroup, type = 'binary', print = FALSE); return()
    }
    if(Rmodel$isDiscrete(nodeGroup)) {
        if(Rmodel$getDistribution(nodeGroup) == 'dmulti') {
            conf$addSampler(target = nodeGroup, type = 'RW_multinomial', print = FALSE); return()
        }
        conf$addSampler(target = nodeGroup, type = 'slice', print = FALSE); return()
    }
    if(length(Rmodel$expandNodeNames(nodeGroup, returnScalarComponents = TRUE)) > 1) {
        conf$addSampler(target = nodeGroup, type = 'RW_block', print = FALSE, silent = TRUE); return()
    }
    conf$addSampler(target = nodeGroup, type = 'RW', print = FALSE); return()
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
    if(together) { eval(parCmd) }
    for(xVarName in xVarNames) {
        if(!together) { eval(parCmd) }
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



