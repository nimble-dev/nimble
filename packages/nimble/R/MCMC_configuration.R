
samplerConf <- setRefClass(
    Class = 'samplerConf',
    fields = list(
        name            = 'ANY',
        samplerFunction = 'ANY',
        baseClassName   = 'ANY',
        target          = 'ANY',
        control         = 'ANY',
        targetAsScalar  = 'ANY'
    ),
    methods = list(
        initialize = function(name, samplerFunction, target, control, model) {
            baseClassName <<- environment(environment(samplerFunction)$contains)$className
            if(is.null(baseClassName) || (baseClassName != 'sampler_BASE')) warning('MCMC sampler nimbleFunctions should inherit from (using "contains" argument) base class sampler_BASE.')
            setName(name)
            setSamplerFunction(samplerFunction)
            setTarget(target, model)
            setControl(control)
            if(name == 'crossLevel')   control <<- c(control, list(dependent_nodes = model$getDependencies(target, self = FALSE, stochOnly = TRUE)))  ## special case for printing dependents of crossLevel sampler (only)
        },
        setName = function(name) name <<- name,
        setSamplerFunction = function(fun) samplerFunction <<- fun,
        setTarget = function(target, model) {
            target <<- target
            targetAsScalar <<- model$expandNodeNames(target, returnScalarComponents = TRUE, sort = TRUE)
        },
        setControl = function(control) control <<- control,
        buildSampler = function(model, mvSaved) {
            samplerFunction(model=model, mvSaved=mvSaved, target=target, control=control)
        },
        toStr = function(displayControlDefaults=FALSE, displayNonScalars=FALSE, displayConjugateDependencies=FALSE) {
            tempList <- list()
            tempList[[paste0(name, ' sampler')]] <- paste0(target, collapse = ', ')
            infoList <- c(tempList, control)
            mcmc_listContentsToStr(infoList, displayControlDefaults, displayNonScalars, displayConjugateDependencies)
        },
        show = function() {
            cat(toStr())
        }
    )
)


## NOTE: methods are documented as a "docstring" with each method - see 'removeSamplers' below. roxygen will automatically grab info from these docstrings and inject into the Rd in the Methods Section
## NOTE: including the name of the class in @aliases is important because by default we only get help("MCMCconf-class") and not help(MCMCconf)
## NOTE: the empty lines are important in the final formatting, so please don't remove any of them in your own help info

#' Class \code{MCMCconf}
#' @aliases MCMCconf addSampler removeSamplers setSamplers printSamplers getSamplers setSamplerExecutionOrder getSamplerExecutionOrder addMonitors addMonitors2 setMonitors setMonitors2 resetMonitors getMonitors getMonitors2 printMonitors setThin setThin2
#' @export
#' @description
#' Objects of this class configure an MCMC algorithm, specific to a particular model.  Objects are normally created by calling \code{\link{configureMCMC}}.
#' Given an MCMCconf object, the actual MCMC function can be built by calling \code{\link{buildMCMC}(conf)}.
#' See documentation below for method initialize() for details of creating an MCMCconf object.
#' @author Daniel Turek
#' @seealso \code{\link{configureMCMC}}
#' @examples
#' code <- nimbleCode({
#'  mu ~ dnorm(0, 1)
#'  x ~ dnorm(mu, 1)
#' })
#' Rmodel <- nimbleModel(code)
#' conf <- configureMCMC(Rmodel)
#' conf$setSamplers(1)
#' conf$addSampler(target = 'x', type = 'slice', control = list(adaptInterval = 100))
#' conf$addMonitors('mu')
#' conf$addMonitors2('x')
#' conf$setThin(5)
#' conf$setThin2(10)
#' conf$printMonitors()
#' conf$printSamplers()
MCMCconf <- setRefClass(
    
    Class = 'MCMCconf',                           
    
    fields = list(
        model               = 'ANY',
        monitors            = 'ANY',
        monitors2           = 'ANY',
        thin                = 'ANY',
        thin2               = 'ANY',
        enableWAIC          = 'ANY',
        controlWAIC         = 'ANY',
        samplerConfs        = 'ANY',
        samplerExecutionOrder = 'ANY',
        controlDefaults     = 'ANY',
        unsampledNodes      = 'ANY',
        postPredSamplerDownstreamNodes = 'ANY',
        ##namedSamplerLabelMaker = 'ANY',  ## usage long since deprecated (Dec 2020)
        mvSamples1Conf      = 'ANY',
        mvSamples2Conf      = 'ANY'
    ),
    
    methods = list(
        
        initialize = function(model, nodes, control = list(), ##rules,
            monitors,                thin  = 1,
            monitors2 = character(), thin2 = 1,
            useConjugacy = getNimbleOption('MCMCuseConjugacy'),
            onlyRW = FALSE,
            onlySlice = FALSE,
            multivariateNodesAsScalars = getNimbleOption('MCMCmultivariateNodesAsScalars'),
            enableWAIC = getNimbleOption('MCMCenableWAIC'), controlWAIC = list(),
            print = TRUE, ...) {
            '
Creates a MCMC configuration for a given model.  The resulting object is suitable as an argument to buildMCMC.

Arguments:

model: A NIMBLE model object, created from nimbleModel(...)

nodes: An optional character vector, specifying the nodes for which samplers should be created.
Nodes may be specified in their indexed form, \'y[1, 3]\', or nodes specified without indexing will be expanded fully, e.g., \'x\' will be expanded to \'x[1]\', \'x[2]\', etc.
If missing, the default value is all non-data stochastic nodes.
If NULL, then no samplers are added.

control: An optional list of control arguments to sampler functions.  If a control list is provided, the elements will be provided to all sampler functions which utilize the named elements given.
For example, the standard Metropolis-Hastings random walk sampler (sampler_RW) utilizes control list elements \'adaptive\', \'adaptInterval\', \'scale\'.
The default values for control list arguments for samplers (if not otherwise provided as an argument to configureMCMC() or addSampler()) are contained in the setup code of each sampling algorithm.

monitors: A character vector of node names or variable names, to record during MCMC sampling.
This set of monitors will be recorded with thinning interval \'thin\', and the samples will be stored into the \'mvSamples\' object.
The default value is all top-level stochastic nodes of the model -- those having no stochastic parent nodes.

monitors2: A character vector of node names or variable names, to record during MCMC sampling.
This set of monitors will be recorded with thinning interval \'thin2\', and the samples will be stored into the \'mvSamples2\' object.
The default value is an empty character vector, i.e. no values will be recorded.

thin: The thinning interval for \'monitors\'.  Default value is one.

thin2: The thinning interval for \'monitors2\'.  Default value is one.

useConjugacy: A logical argument, with default value TRUE.  If specified as FALSE, then no conjugate samplers will be used, even when a node is determined to be in a conjugate relationship.

onlyRW: A logical argument, with default value FALSE.  If specified as TRUE, then Metropolis-Hastings random walk samplers will be assigned for all non-terminal continuous-valued nodes nodes. Discrete-valued nodes are assigned a slice sampler, and terminal nodes are assigned a posterior_predictive sampler.

onlySlice: A logical argument, with default value FALSE.  If specified as TRUE, then a slice sampler is assigned for all non-terminal nodes. Terminal nodes are still assigned a posterior_predictive sampler.

multivariateNodesAsScalars: A logical argument, with default value FALSE.  If specified as TRUE, then non-terminal multivariate stochastic nodes will have scalar samplers assigned to each of the scalar components of the multivariate node.  The default value of FALSE results in a single block sampler assigned to the entire multivariate node.  Note, multivariate nodes appearing in conjugate relationships will be assigned the corresponding conjugate sampler (provided useConjugacy == TRUE), regardless of the value of this argument.

enableWAIC: A logical argument, specifying whether to enable WAIC calculations for the resulting MCMC algorithm.  Defaults to the value of nimbleOptions(\'MCMCenableWAIC\'), which in turn defaults to FALSE.  Setting nimbleOptions(\'MCMCenableWAIC\' = TRUE) will ensure that WAIC is enabled for all calls to \`configureMCMC\` and \`buildMCMC\`.

controlWAIC A named list of inputs that control the behavior of the WAIC calculation, passed as the \'control\' input to \'buildWAIC\'. See \'help(waic)\`.

print: A logical argument specifying whether to print the montiors and samplers.  Default is TRUE.

...: Additional named control list elements for default samplers, or additional arguments to be passed to the autoBlock function when autoBlock = TRUE.
'
            if(is(model, 'RmodelBaseClass')) {
                model <<- model
            } else if(is(model, 'CmodelBaseClass')) {
                model <<- model$Rmodel
            } else stop('\'model\' must be a compiled or un-compiled NIMBLE model object')
            monitors  <<- character()
            monitors2 <<- character()
            addMonitors( monitors,  print = FALSE)
            addMonitors2(monitors2, print = FALSE)
            setThin( thin,  print = FALSE)
            setThin2(thin2, print = FALSE)
            enableWAIC <<- enableWAIC
            controlWAIC <<- controlWAIC
            samplerConfs <<- list()
            samplerExecutionOrder <<- numeric()
            controlDefaults <<- list(...)
            ##namedSamplerLabelMaker <<- labelFunctionCreator('namedSampler')  ## usage long since deprecated (Dec 2020)
            for(i in seq_along(control))     controlDefaults[[names(control)[i]]] <<- control[[i]]
            if(missing(nodes)) {
                nodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
                # Check of all(model$isStoch(nodes)) is not needed in this case
            } else {
                if(is.null(nodes) || length(nodes)==0)     nodes <- character(0)
            }
            
            addDefaultSampler(nodes = nodes,
                              useConjugacy = useConjugacy,
                              onlyRW = onlyRW,
                              onlySlice = onlySlice,
                              multivariateNodesAsScalars = multivariateNodesAsScalars,
                              print = FALSE)

            if(print)   show()    ##printSamplers()
        },

        addDefaultSampler = function(nodes = character(),
                                     control = list(),
                                     useConjugacy = getNimbleOption('MCMCuseConjugacy'),
                                     onlyRW = FALSE,
                                     onlySlice = FALSE,
                                     multivariateNodesAsScalars = getNimbleOption('MCMCmultivariateNodesAsScalars'),
                                     print = TRUE,
                                     ...) {
            '
For internal use.  Adds default MCMC samplers to the specified nodes.
'
            useNewConfigureMCMC <- isTRUE(nimbleOptions("useNewConfigureMCMC"))
            
            controlDefaultsArg <- list(...)
            for(i in seq_along(control))     controlDefaultsArg[[names(control)[i]]] <- control[[i]]

            if(!is.character(nodes))   stop('nodes argument must be a character vector of model nodes or variables')
            nl_checkVarNamesInModel(model, removeIndexing(nodes))
            nodes <- model$expandNodeNames(nodes)
            if(useNewConfigureMCMC) {
                if(!(all(model$isStoch(nodes)))) {
                    stop('assigning samplers to non-stochastic nodes: ',
                         paste0(nodes[!model$isStoch(nodes)],
                                collapse=', ')) }    ## ensure all target node(s) are stochastic
            }
            nodes <- model$topologicallySortNodes(nodes)   ## topological sort
            if(getNimbleOption('MCMCorderPosteriorPredictiveSamplersLast')) {
                postPredBool <- nodes %in% model$getNodeNames(predictiveOnly = TRUE)
                nodes <- c(nodes[!postPredBool], nodes[postPredBool])  ## put posterior predictive nodes at the end
            }

            if(!useNewConfigureMCMC) {
                if(!(all(model$isStoch(nodes)))) { stop('assigning samplers to non-stochastic nodes: ', paste0(nodes[!model$isStoch(nodes)], collapse=', ')) }
                if(useConjugacy) conjugacyResultsAll <- model$checkConjugacy(nodes)
            } else {
                ## convert to node IDs:
                nodeIDs <- model$expandNodeNames(nodes, returnType = 'ids')
                nodeIDsOrig <- nodeIDs
                
                ## determine which posterior predictive nodes should be sampled with posterior_predictive sampler.
                ## this requires some care, because it's only those nodes for which all
                ## downstream dependents are also slated for sampling.
                ## upon finding one, the root node and all downstream dependents are then
                ## removed from the set of node (ids) for later sampler assignment.
                ## controlled by system option: MCMCusePosteriorPredictiveSampler
                nodesForPosteriorPredictiveSampler <- character()
                if(getNimbleOption('MCMCusePosteriorPredictiveSampler')) {
                    predictiveNodeIDsToSample <- numeric()
                    predictiveNodeIDs <- model$getPredictiveNodeIDs()
                    predictiveNodeIDsToCheck <- intersect(predictiveNodeIDs, nodeIDsOrig)    ## only check predictive nodes which are slated for sampling
                    isEndNodeBool <- model$isEndNode(predictiveNodeIDsToCheck)
                    predictiveNodeIDsToCheckEnd <- predictiveNodeIDsToCheck[isEndNodeBool]
                    predictiveNodeIDsToCheckNonEnd <- predictiveNodeIDsToCheck[!isEndNodeBool]
                    ## the following while-loop checks downstream dependencies of *non-end* predictive nodes
                    nPredictiveNodesToCheck <- length(predictiveNodeIDsToCheckNonEnd)   ## only loop over and check *non-end* predictive nodes
                    nextPredNodeInd <- 1
                    while(nextPredNodeInd <= nPredictiveNodesToCheck) {
                        nid <- as.numeric(predictiveNodeIDsToCheckNonEnd[nextPredNodeInd])
                        downstreamNoSelfIDs <- model$getDependencies(nid, self = FALSE, stochOnly = TRUE, downstream = TRUE, returnType = 'ids')
                        ## quick reality check:
                        if(!all(downstreamNoSelfIDs %in% predictiveNodeIDs))   stop('predictive node IDs in model appear to be set wrong')
                        ## skip nodes if the entire downstream network wasn't slated for sampling:
                        if(!all(c(nid, downstreamNoSelfIDs) %in% nodeIDsOrig))   { nextPredNodeInd <- nextPredNodeInd + 1;   next }
                        ## found a posterior predictive (non-end) node to sample:
                        predictiveNodeIDsToSample <- c(predictiveNodeIDsToSample, nid)
                        ## remove all downstream nodes from the nodeIDs for future sampler assignment:
                        nodeIDs <- setdiff(nodeIDs, downstreamNoSelfIDs)
                        ## remove up-to-and-including this node from the set of predictive (non-end) nodes to check:
                        predictiveNodeIDsToCheckNonEnd <- predictiveNodeIDsToCheckNonEnd[-(1:nextPredNodeInd)]
                        ## remove all downstream dependencies of this node from *both* sets (end, and non-end) of predictive nodes:
                        predictiveNodeIDsToCheckEnd    <- setdiff(predictiveNodeIDsToCheckEnd,    downstreamNoSelfIDs)
                        predictiveNodeIDsToCheckNonEnd <- setdiff(predictiveNodeIDsToCheckNonEnd, downstreamNoSelfIDs)
                        nPredictiveNodesToCheck <- length(predictiveNodeIDsToCheckNonEnd)
                        nextPredNodeInd <- 1
                    }
                    ## now, any remaining *end* predictive nodes should necessarily receive posterior_predictive samplers:
                    predictiveNodeIDsToSample <- c(predictiveNodeIDsToSample, predictiveNodeIDsToCheckEnd)
                    ## convert back to node names:
                    nodes <- model$modelDef$maps$graphID_2_nodeName[nodeIDs]
                    nodesForPosteriorPredictiveSampler <- model$modelDef$maps$graphID_2_nodeName[predictiveNodeIDsToSample]
                }
                
                if(useConjugacy) conjugacyResultsAll <- nimble:::conjugacyRelationshipsObject$checkConjugacy(model, nodeIDs) ## Later, this can go through model$checkConjugacy if we make it check whether nodes are already nodeIDs.  To isolate changes, I am doing it directly here.
                nodeDeclIDs <- model$modelDef$maps$graphID_2_declID[nodeIDs] ## Below, nodeDeclIDs[i] gives the nodeDeclID.  We could add an interface to get this.
                nodeDeclID_2_nodes <- split(nodes, nodeDeclIDs)
                
                uniqueNodeDeclIDs <- unique(nodeDeclIDs)
                nodeTraits <- lapply(uniqueNodeDeclIDs,
                                     function(x) {
                                         declInfo <- model$modelDef$declInfo[[x]]
                                         dist <- declInfo$distributionName
                                         distInfo <- getDistributionInfo(dist)
                                         discrete <- distInfo$discrete
                                         ## Following can be replaced by an efficiency version model$isBinary 
                                         binary <- dist == 'dbern'
                                         ## If dist == 'dbin', then binary-ness will be checked for each node, below
                                         ## This could be improved to see if they all have a literal "1", for example.
                                         ## The checking below is inefficient!
                                         ## For nodeScalarComponents, we will check a single node
                                         ## We could use returnType = 'ids', but we have a warning generated in that case,
                                         ## for future investigation.
                                         
                                         ## Determining nodeLength is a bit tricky.
                                         ## The only purpose is to determine scalar vs. non-scalar.
                                         ## In the future, we may want to make this available in distributionInfo.
                                         ## A difficulty is whether it is possible to declare a scalar
                                         ##     case of a multivariate node.
                                         ## If so, then the status (scalar vs non-scalar) needs to be checked
                                         ## node by node, not just once for the declaration.
                                         ## Here we use a system that marks scalar as scalar.  The alternative is really
                                         ## "maybe non-scalar", in which case it is checked node-by-node below.
                                         ## Unfortunately, this processing works from nodeNames.
                                         allNodeScalarComponents <- model$expandNodeNames(
                                             nodeDeclID_2_nodes[[as.character(x)]],
                                             returnScalarComponents = TRUE)
                                         nodeLength <- if(length(allNodeScalarComponents) ==
                                                          length(nodeDeclID_2_nodes[[as.character(x)]]))
                                                           1 else 2 ## 2 indicates "maybe non-scalar"
                                         list(dist = dist,
                                              discrete = discrete,
                                              binary = binary,
                                              nodeLength = nodeLength)
                                     }
                                     )
                names(nodeTraits) <- as.character(uniqueNodeDeclIDs)
                
                allDists <- unlist(lapply(model$modelDef$declInfo, `[[`, 'distributionName'))
                allDists <- allDists[!is.na(allDists)]
                check_dCRP <- any(allDists == "dCRP")
                
                clusterNodeInfo <- NULL; dcrpNode <- NULL; numCRPnodes <- 0; clusterNodeParams <- NULL

                for(i in seq_along(nodes)) {
                    node <- nodes[i]
                    if(!useNewConfigureMCMC) {
                        discrete <- model$isDiscrete(node)
                        binary <- model$isBinary(node)
                        nodeDist <- model$getDistribution(node)
                        nodeScalarComponents <- model$expandNodeNames(node, returnScalarComponents = TRUE)
                        nodeLength <- length(nodeScalarComponents)
                    } else {
                        nodeDeclID <- nodeDeclIDs[i]
                        nodeTrait <- nodeTraits[[as.character(nodeDeclID)]] ## from split, the names are nodeDeclIds
                        nodeDist <- nodeTrait$dist
                        discrete <- nodeTrait$discrete
                        if(nodeDist != "dbin") {
                            binary <- nodeTrait$binary
                        } else {
                            ## Check dbin case one by one, since the paramExpr may have come from
                            ## a constants replacement, with a 1 for some nodes of a declaration but not others.
                            binary <- model$getParamExpr(node, 'size') == 1
                        }
                        nodeLength <- nodeTrait$nodeLength
                        if(nodeLength == 2) { ## code for "maybe non-scalar", so we check each one
                            nodeScalarComponents <- model$expandNodeNames(node, returnScalarComponents = TRUE)
                            nodeLength <- length(nodeScalarComponents)
                        }
                    }

                    ## if node is the root of a posterior predictive (entirely non-data) network of nodes, assign 'posterior_predictive' sampler
                    if(node %in% nodesForPosteriorPredictiveSampler) { addSampler(target = node, type = 'posterior_predictive', control = controlDefaultsArg);     next }
                    
                    ## for multivariate nodes, either add a conjugate sampler, RW_multinomial, or RW_block sampler
                    if(nodeLength > 1) {
                        if(useConjugacy) {
                            conjugacyResult <- conjugacyResultsAll[[node]]
                            if(!is.null(conjugacyResult)) {
                                addConjugateSampler(conjugacyResult = conjugacyResult,
                                                    dynamicallyIndexed = model$modelDef$varInfo[[model$getVarNames(nodes=node)]]$anyDynamicallyIndexed);     next }
                        }
                        if(nodeDist == 'dmulti')              { addSampler(target = node, type = 'RW_multinomial', control = controlDefaultsArg);     next }
                        if(nodeDist == 'ddirch')              { addSampler(target = node, type = 'RW_dirichlet',   control = controlDefaultsArg);     next }
                        if(nodeDist == 'dwish')               { addSampler(target = node, type = 'RW_wishart',     control = controlDefaultsArg);     next }
                        if(nodeDist == 'dinvwish')            { addSampler(target = node, type = 'RW_wishart',     control = controlDefaultsArg);     next }
                        if(nodeDist == 'dlkj_corr_cholesky')  {
                            if(nodeLength >= 9) {
                                addSampler(target = node, type = 'RW_block_lkj_corr_cholesky', control = controlDefaultsArg)
                            } else {
                                if(nodeLength == 4) {
                                    addSampler(target = node, type = 'RW_lkj_corr_cholesky', control = controlDefaultsArg)  ## only a scalar free param in 2x2 case
                                } else warning("Not assigning sampler to dlkj_corr_cholesky node for 1x1 case.")
                            }
                            next
                        }
                        if(nodeDist == 'dcar_normal')         { addSampler(target = node, type = 'CAR_normal', control = controlDefaultsArg);         next }
                        if(nodeDist == 'dcar_proper')         { addSampler(target = node, type = 'CAR_proper', control = controlDefaultsArg);         next }
                        if(nodeDist == 'dCRP')                {
                            numCRPnodes <- numCRPnodes + 1
                            clusterNodeInfo[[numCRPnodes]] <- findClusterNodes(model, node)
                            controlCRP <- controlDefaultsArg
                            controlCRP$checkConjugacy <- useConjugacy
                            controlCRP$clusterVarInfo <- clusterNodeInfo[[numCRPnodes]]
                            addSampler(target = node, type = 'CRP', control = controlCRP)
                            dcrpNode[numCRPnodes] <- node
                            next
                        }
                        if(multivariateNodesAsScalars) {
                            for(scalarNode in nodeScalarComponents) {
                                if(onlySlice) addSampler(target = scalarNode, type = 'slice', control = controlDefaultsArg)
                                else          addSampler(target = scalarNode, type = 'RW',    control = controlDefaultsArg)    };     next }
                        addSampler(target = node, type = 'RW_block', silent = TRUE, control = controlDefaultsArg);     next }
                    
                    if(onlyRW && !discrete)   { addSampler(target = node, type = 'RW',    control = controlDefaultsArg);     next }
                    if(onlySlice)             { addSampler(target = node, type = 'slice', control = controlDefaultsArg);     next }
                    
                    ## if node passes checkConjugacy(), assign 'conjugate_dxxx' sampler
                    if(useConjugacy) {
                        conjugacyResult <- conjugacyResultsAll[[node]]
                        if(!is.null(conjugacyResult)) {
                            addConjugateSampler(conjugacyResult = conjugacyResult,
                                                dynamicallyIndexed = model$modelDef$varInfo[[model$getVarNames(nodes=node)]]$anyDynamicallyIndexed);     next }
                    }
                    
                    ## if node is discrete 0/1 (binary), assign 'binary' sampler
                    if(binary) { addSampler(target = node, type = 'binary', control = controlDefaultsArg);     next }
                    
                    ## for categorical nodes, assign a 'categorical' sampler
                    if(nodeDist == 'dcat') { addSampler(target = node, type = 'categorical', control = controlDefaultsArg);     next }
                    
                    ## if node distribution is discrete, assign 'slice' sampler
                    if(discrete) { addSampler(target = node, type = 'slice', control = controlDefaultsArg);     next }
                    
                    ## if node distribution is dgamma and its dependency is dCRP, assign 'CRP_concentration' sampler
                    if(check_dCRP) {
                        if(nodeDist == 'dgamma'){
                            depNode <- model$getDependencies(node, self=FALSE)
                            if(length(depNode) == 1) {
                                depNodeDist <- model$getDistribution(depNode)
                                if(!is.na(depNodeDist[1]) & depNodeDist[1] == 'dCRP'){  ## depNodeDist should be length 1
                                    addSampler(target = node, type = 'CRP_concentration', control = controlDefaultsArg)
                                    next
                                }
                            }
                        }
                    }
                    
                    ## default: 'RW' sampler
                    addSampler(target = node, type = 'RW', control = controlDefaultsArg);     next
                }

                ## For CRP-based models, wrap samplers for cluster parameters so not sampled if cluster is unoccupied,
                ## and check for non-fixed (hyper)parameters of cluster parameters and assign special slice sampler that knows how
                ## to determine dependencies dynamically.
                ## If anything contraindicates wrapping, we avoid it. E.g., a dangerous case is if a single hyperparameter is involved
                ## in multiple sets of cluster parameters, but another hyperparameter is involved in one of those sets.
                ## In that case, the processing of second hyperparameter could turn on the wrapping for some cluster parameters
                ## which would mess up sampling of the first hyperparameter. 
                wrap <- TRUE
                if(!is.null(clusterNodeInfo)) {
                    allClusterNodes <- lapply(clusterNodeInfo, function(x) x$clusterNodes)
                    allClusterNodesVec <- unlist(allClusterNodes)
                    for(k in seq_along(clusterNodeInfo)) {
                        for(idx in seq_along(clusterNodeInfo[[k]]$clusterNodes)) {

                            clusterNodes <- clusterNodeInfo[[k]]$clusterNodes[[idx]]
                            clusterNodeParams <- model$getParents(clusterNodes, stochOnly = TRUE, includeData = FALSE)
                            clusterNodeParams <- clusterNodeParams[!clusterNodeParams %in% allClusterNodesVec]

                            ## Avoid cases where there is other stochastic indexing (that might or might not use dCRP) but also
                            ## indexes the cluster nodes, e.g., mu[xi[i]] and mu[eta[i]] as hard to determine if cluster is occupied.
                            clusterNodeDeps <- model$getDependencies(clusterNodes, stochOnly = TRUE, self = FALSE)
                            clusterNodeDeps <- clusterNodeDeps[!clusterNodeDeps %in% allClusterNodesVec]
                            if(!all(clusterNodeDeps %in%
                                    model$getDependencies(dcrpNode[[k]], stochOnly = TRUE, self = FALSE)))
                                wrap <- FALSE
                                
                            ## For now avoid wrapper if any overlap of clusterNodes, as hard to determine if cluster is occupied.
                            ## We'll need to come back to this to handle the mu[xi[i],eta[j]] case if we want to
                            ## avoid sampling empty clusters in that case.
                            if(length(allClusterNodes) > 1 && any(clusterNodes %in% unlist(allClusterNodes[-k])))
                                wrap <- FALSE

                            ## Avoid cases where a parameter is involved in multiple sets of cluster nodes.
                            for(cnt in seq_along(clusterNodeParams)) {
                                depNodes <- model$getDependencies(clusterNodeParams[cnt],
                                                                  stochOnly = TRUE, self = FALSE, includeData = TRUE)
                                if(!identical(sort(depNodes), sort(clusterNodes)))
                                    wrap <- FALSE
                            }
                        }
                    }
                    if(wrap) {
                        for(k in seq_along(clusterNodeInfo)) {
                            for(idx in seq_along(clusterNodeInfo[[k]]$clusterNodes)) {

                                clusterNodes <- clusterNodeInfo[[k]]$clusterNodes[[idx]]
                                clusterNodeParams <- model$getParents(clusterNodes, stochOnly = TRUE, includeData = FALSE)
                                clusterNodeParams <- clusterNodeParams[!clusterNodeParams %in% allClusterNodesVec]

                                for(cnt in seq_along(clusterNodeParams)) {
                                    removeSamplers(clusterNodeParams[cnt])
                                    controlCRP <- controlDefaultsArg
                                    controlCRP$dcrpNode <- dcrpNode[[k]]
                                    controlCRP$clusterNodes <- clusterNodes
                                    controlCRP$clusterIDs <- clusterNodeInfo[[k]]$clusterIDs[[idx]]
                                    addSampler(clusterNodeParams[cnt], 'slice_CRP_base_param', control = controlCRP)
                                }

                                samplers <- getSamplers(clusterNodes)
                                removeSamplers(clusterNodes)
                                for(i in seq_along(samplers)) {
                                    node <- samplers[[i]]$target
                                    ## getSamplers() returns samplers in order of configuration not in order of input.
                                    clusterID <- which(clusterNodes == node)
                                    if(length(clusterID) != 1)
                                        stop("Cannot determine wrapped sampler for cluster parameter ", node, ".")
                                    controlCRP <- controlDefaultsArg
                                    controlCRP$wrapped_type <- samplers[[i]]$name
                                    controlCRP$wrapped_conf <- samplers[[i]]
                                    controlCRP$dcrpNode <- dcrpNode[[k]]
                                    controlCRP$clusterID <- clusterNodeInfo[[k]]$clusterIDs[[idx]][clusterID]
                                    addSampler(target = node, type = 'CRP_cluster_wrapper', control = controlCRP)
                                }
                            }
                        }
                    }
                }
            }
            
            setUnsampledNodes()
            if(print)   printSamplers(byType = TRUE)   ##show()    ##printSamplers()
        },

        addConjugateSampler = function(conjugacyResult, dynamicallyIndexed = FALSE, print = FALSE) {
            ## update May 2016: old (non-dynamic) system is no longer supported -DT
            ##if(!getNimbleOption('useDynamicConjugacy')) {
            ##    addSampler(target = conjugacyResult$target, type = conjugacyResult$type, control = conjugacyResult$control)
            ##    return(NULL)
            ##}
            prior <- conjugacyResult$prior
            dependentCounts <- sapply(conjugacyResult$control, length)
            names(dependentCounts) <- gsub('^dep_', '', names(dependentCounts))
            ## we now have separate sampler functions for the same conjugacy
            ## when there are both dynamically and non-dynamically indexed
            ## nodes with that conjugacy, as the dynamically-indexed
            ## sampler has a check for zero contributions for each dependent
            conjSamplerName <- createDynamicConjugateSamplerName(prior = prior, dependentCounts = dependentCounts, dynamicallyIndexed = dynamicallyIndexed)
            if(!dynamicConjugateSamplerExists(conjSamplerName)) {
                conjSamplerDef <- conjugacyRelationshipsObject$generateDynamicConjugateSamplerDefinition(prior = prior, dependentCounts = dependentCounts, doDependentScreen = dynamicallyIndexed)
                dynamicConjugateSamplerAdd(conjSamplerName, conjSamplerDef)
            }
            conjSamplerFunction <- dynamicConjugateSamplerGet(conjSamplerName)
            nameToPrint <- gsub('^sampler_', '', conjSamplerName)
            addSampler(target = conjugacyResult$target, type = conjSamplerFunction, control = conjugacyResult$control, print = print, name = nameToPrint)
        },
        
        addSampler = function(target = character(),    ## target argument is *not* expanded (unless expandTarget = TRUE)
                              type = 'RW',
                              control = list(),
                              print = NULL,        ## default value print is TRUE when default=TRUE, and FALSE otherwise
                              name,
                              scalarComponents = FALSE,
                              expandTarget = FALSE,
                              silent = FALSE,
                              default = FALSE,
                              useConjugacy = getNimbleOption('MCMCuseConjugacy'),
                              onlyRW = FALSE,
                              onlySlice = FALSE,
                              multivariateNodesAsScalars = getNimbleOption('MCMCmultivariateNodesAsScalars'),
                              ...) {
            '
Adds a sampler to the list of samplers contained in the MCMCconf object.

Arguments:

target: The target node or nodes to be sampled.  This may be specified as a character vector of model node and/or variable names.  For univariate samplers, only a single target node should be provided (unless \'expandTarget\' is TRUE).  For multivariate samplers, one instance of the multivariate sampler will be assigned to all nodes specified.  Nodes are specified in combination with the \'expandTarget\' and \'scalarComponents\' arguments.  See details.

type: When \'default\' is FALSE, specifies the type of sampler to add, specified as either a character string or a nimbleFunction object.  If the character argument type=\'newSamplerType\', then either newSamplerType or sampler_newSamplertype must correspond to a nimbleFunction (i.e. a function returned by nimbleFunction, not a specialized nimbleFunction).  Alternatively, the type argument may be provided as a nimbleFunction itself rather than its name.  In that case, the \'name\' argument may also be supplied to provide a meaningful name for this sampler.  The default value is \'RW\' which specifies scalar adaptive Metropolis-Hastings sampling with a normal proposal distribution. This default will result in an error if \'target\' specifies more than one target node (unless \'expandTarget\' is TRUE).  This argument is not used when the \'default\' argument is TRUE.

control: An optional list of control arguments to sampler functions.  These will override those specified in the control list argument to configureMCMC.  If a control list is provided, the elements will be provided to all sampler functions which utilize the named elements given. For example, the standard Metropolis-Hastings random walk sampler (sampler_RW) utilizes control list elements \'adaptive\', \'adaptInterval\', \'scale\'. The default values for control list arguments for samplers (if not otherwise provided as an argument to configureMCMC or addSampler) are contained in the setup code of each sampling algorithm.

print: Logical argument, specifying whether to print the details of newly added sampler(s).

name: Optional character string name for the sampler, which is used by the printSamplers method.  If \'name\' is not provided, the \'type\' argument is used to generate the sampler name.

scalarComponents: Logical argument, with default value FALSE.  When \'expandTarget\' is TRUE and therefore \'target\' undergoes expansion via \'expandNodeNames\', this argument is passed as the \'returnScalarComponents\' argument to \'expandNodeNames\'.  This has the effect of returning the constituent scalar components which make up any multivariate nodes in \'target\', in which case a separate instance of the specified sampler type will be added to each scalar component.  This argument is only used when \'expandTarget\' is TRUE.  See details.

expandTarget: Logical argument, indicating whether to expand the \'target\' argument using \'expandNodeNames\' into the underlying constitutent model nodes.  These may include univariate and/or multivariate nodes, although additional control is provided using the \'scalarComponents\' argument.  When \'expandTarget\' is TRUE and this expansion to individual nodes occurs, samplers will be assigned independently to each (univariate or multivariate) node component of \'target\'.  That is, distinct instances of the specified sampler type will be assigned to each of the nodes which comprise \'target\'.  If \'target\' is comprised of multiple univariate nodes, then a univariate sampler should be specified, and similarly when \'target\' expands to one or more multivariate nodes, then a multivariate sampler should be specified.

silent: Logical argument, specifying whether to print warning messages when assigning samplers.

default: Logical argument, with default value FALSE.  When FALSE, the \'type\' argument dictates what sampling algorithm is assigned to the specified nodes.  When TRUE, default samplers will be assigned to the specified nodes following the same logic as the configureMCMC method, and also using the \'useConjugacy\', \'onlyRW\', \'onlySlice\' and \'multivariateNodesAsScalars\' arguments.

useConjugacy: Logical argument, with default value TRUE.  If specified as FALSE, then no conjugate samplers will be used, even when a node is determined to be in a conjugate relationship.  This argument is only used when the \'default\' argument is TRUE.

onlyRW: Logical argument, with default value FALSE.  If specified as TRUE, then Metropolis-Hastings random walk samplers will be assigned for all non-terminal continuous-valued nodes nodes. Discrete-valued nodes are assigned a slice sampler, and terminal nodes are assigned a posterior_predictive sampler.  This argument is only used when the \'default\' argument is TRUE.

onlySlice: Logical argument, with default value FALSE.  If specified as TRUE, then a slice sampler is assigned for all non-terminal nodes. Terminal nodes are still assigned a posterior_predictive sampler.  This argument is only used when the \'default\' argument is TRUE.

multivariateNodesAsScalars: Logical argument, with default value FALSE.  If specified as TRUE, then non-terminal multivariate stochastic nodes will have scalar samplers assigned to each of the scalar components of the multivariate node.  The default value of FALSE results in a single block sampler assigned to the entire multivariate node.  Note, multivariate nodes appearing in conjugate relationships will be assigned the corresponding conjugate sampler (provided \'useConjugacy\' is TRUE), regardless of the value of this argument.  This argument is only used when the \'default\' argument is TRUE.

...: Additional named arguments passed through ... will be used as additional control list elements.

Details:

Samplers will be assigned to nodes specified by the \'target\' argument.  The \'target\' argument does not undergo expansion, unless \'expandTarget\' is TRUE, in which case \'target\' is expanded to the underlying (univariate or multivariate) nodes, and a distinct sampler is added to each.  When \'expandTarget\' is TRUE and \'scalarComponents\' is TRUE, then any multivariate nodes are further decomponsed into their underlying scalar components, and a separate sampler is added to each.

Samplers are added added to the end of the list of samplers for this MCMCconf object, and do not replace any exisiting samplers.  Samplers are removed using the removeSamplers method.

Invisibly returns a list of the current sampler configurations, which are samplerConf reference class objects.
'
            if(length(target) == 0)   return(invisible(samplerConfs))

            if(default) {
                addDefaultSampler(nodes = target,
                                  control = control,
                                  useConjugacy = useConjugacy,
                                  onlyRW = onlyRW,
                                  onlySlice = onlySlice,
                                  multivariateNodesAsScalars = multivariateNodesAsScalars,
                                  print = if(is.null(print)) TRUE else print,   ## default of print is TRUE when adding default sampler
                                  ...)
                return(invisible(samplerConfs))
            }

            print <- if(is.null(print)) FALSE else print                        ## default of print is FALSE otherwise
            
            nameProvided <- !missing(name)
            if(is.character(type)) {
                if(type == 'conjugate') {
                    conjugacyResult <- model$checkConjugacy(target)[[target]]
                    if(!is.null(conjugacyResult)) {
                        varName <- model$getVarNames(nodes = target)
                        if(length(varName) > 1) stop("MCMCconf: conjugate sampler for more than one variable: ", target, ".")
                        return(addConjugateSampler(conjugacyResult = conjugacyResult,
                                                   dynamicallyIndexed = model$modelDef$varInfo[[varName]]$anyDynamicallyIndexed,
                                                   print = print))
                    } else stop(paste0('Cannot assign conjugate sampler to non-conjugate node: \'', target, '\''))
                }
                thisSamplerName <- if(nameProvided) name else gsub('^sampler_', '', type)   ## removes 'sampler_' from beginning of name, if present
                if(thisSamplerName == 'RW_block' && !silent) {
                    messageIfVerbose('  [Note] Assigning an RW_block sampler to nodes with very different scales can result in low MCMC efficiency.  If all nodes assigned to RW_block are not on a similar scale, we recommend providing an informed value for the \"propCov\" control list argument, or using the AFSS sampler instead.')
                }
                if(thisSamplerName %in% c("RW_PF", "RW_PF_block")) {
                    if (!("nimbleSMC" %in% (installed.packages()[,"Package"]))) {
                        stop(paste0("Particle filters have been moved to the `nimbleSMC` package. ",
                                    "Install and load `nimbleSMC` to use them."))
                    } else if (!("nimbleSMC" %in% .packages())) {
                        stop("`nimbleSMC` must be loaded to use particle filtering samplers.")
                    }
                }
                if(exists(type, inherits = TRUE) && is.nfGenerator(eval(as.name(type)))) {   ## try to find sampler function 'type'
                    samplerFunction <- eval(as.name(type))
                } else {
                    sampler_type <- paste0('sampler_', type)   ## next, try to find sampler function 'sampler_type'
                    if(exists(sampler_type) && is.nfGenerator(eval(as.name(sampler_type)))) {   ## try to find sampler function 'sampler_type'
                        samplerFunction <- eval(as.name(sampler_type))
                    } else stop(paste0('cannot find sampler type \'', type, '\''))
                }
            } else if(is.function(type)) {
                if(nameProvided) {
                    thisSamplerName <- name
                } else {
                    typeArg <- substitute(type)
                    if(is.name(typeArg)) {
                        thisSamplerName <- gsub('^sampler_', '', deparse(typeArg))
                    } else {
                        thisSamplerName <- 'custom_function'
                    }
                }
                samplerFunction <- type
            } else stop('sampler type must be character name or function')
            if(!is.character(thisSamplerName)) stop('sampler name should be a character string')
            if(!is.function(samplerFunction)) stop('sampler type does not specify a function')

            if(!(all(model$isStoch(target)))) { warning(paste0('No sampler assigned to non-stochastic node: ', paste0(target,collapse=', '))); return(invisible(samplerConfs)) }   ## ensure all target node(s) are stochastic

            ##libraryTag <- if(nameProvided) namedSamplerLabelMaker() else thisSamplerName   ## unique tag for each 'named' sampler, internal use only  ## usage long since deprecated (Dec 2020)
            ##if(is.null(controlNamesLibrary[[libraryTag]]))   controlNamesLibrary[[libraryTag]] <<- mcmc_findControlListNamesInCode(samplerFunction)   ## populate control names library
            ##requiredControlNames <- controlNamesLibrary[[libraryTag]]
            controlArgs <- c(control, list(...))
            thisControlList <- mcmc_generateControlListArgument(control=controlArgs, controlDefaults=controlDefaults)  ## should name arguments
            
            if(!expandTarget) {
                ## no node expansion takes place when expandTarget is FALSE
                addSamplerOne(thisSamplerName, samplerFunction, target, thisControlList, print)
            } else {
                ## assign sampler type to each component of target after expanding node names,
                ## also using scalarComponents argument here, to control the node expansion
                targetExpanded <- model$expandNodeNames(target, returnScalarComponents = scalarComponents, sort = TRUE)
                for(i in seq_along(targetExpanded)) {
                    addSamplerOne(thisSamplerName, samplerFunction, targetExpanded[i], thisControlList, print)
                }
            }
            
            return(invisible(samplerConfs))
        },

        addSamplerOne = function(thisSamplerName, samplerFunction, targetOne, thisControlList, print) {
            '
For internal use only
'
            newSamplerInd <- length(samplerConfs) + 1
            samplerConfs[[newSamplerInd]] <<- samplerConf(name=thisSamplerName, samplerFunction=samplerFunction, target=targetOne, control=thisControlList, model=model)
            samplerExecutionOrder <<- c(samplerExecutionOrder, newSamplerInd)
            if(print) printSamplers(newSamplerInd)
        },
        
        removeSamplers = function(..., ind, print = FALSE) {
            '
Removes one or more samplers from an MCMCconf object.

This function also has the side effect of resetting the sampler execution ordering so as to iterate over the remaining set of samplers, sequentially, executing each sampler once.

Arguments:

...: Character node names or numeric indices.  Character node names specify the node names for samplers to remove, or numeric indices can provide the indices of samplers to remove.

ind: A numeric vector or character vector specifying the samplers to remove.  A numeric vector may specify the indices of the samplers to be removed.  Alternatively, a character vector may be used to specify a set of model nodes and/or variables, and all samplers whose \'target\' is among these nodes will be removed.  If omitted, then all samplers are removed.

print: A logical argument specifying whether to print the current list of samplers once the removal has been done (default FALSE).
'
            if(missing(ind)) {
                ind <- list(...)
                ind <- unname(unlist(ind))
                if(is.null(ind))   ind <- seq_along(samplerConfs)
            }
            if(is.character(ind))   ind <- findSamplersOnNodes(ind)
            if(length(ind) > 0 && max(ind) > length(samplerConfs)) stop('MCMC configuration doesn\'t have that many samplers')
            samplerConfs[ind] <<- NULL
            samplerExecutionOrder <<- seq_along(samplerConfs)
            if(print) printSamplers()
            return(invisible(NULL))
        },

        removeSampler = function(...) {
            '
Alias for removeSamplers method
'
            removeSamplers(...)
        },
        
        ## Note: function prototype is identical to addSampler
        replaceSamplers = function(...) {
            '
Replaces one or more samplers from an MCMCconf object with newly specified sampler(s).  Operation and arguments are identical to the \'addSampler\' method, with the additional side effect of first removing any existing samplers which operate on the specified node(s).

This function also has the side effect of resetting the sampler execution ordering so as to iterate over the remaining set of samplers, sequentially, executing each sampler once.

See \'addSamplers\' for a description of the arguments.

This function also has the side effect of resetting the sampler execution ordering so as to iterate over the newly specified set of samplers, sequentially, executing each sampler once.
'
            samplerConfs_save <- samplerConfs
            samplerExecutionOrder_save <- samplerExecutionOrder
            m <- match.call(definition = addSampler)
            removeSamplers(eval(m$target), eval(m$nodes))   ## both unnamed arguments accepted as ... argument
            e <- try(addSampler(...), silent = TRUE)   ## pass all arguments along to addSampler
            ## if an error occurred in addSampler, then restore original samplerConfs
            if(inherits(e, 'try-error')) {
                samplerConfs <<- samplerConfs_save
                samplerExecutionOrder <<- samplerExecutionOrder_save
                errorMessage <- sub('^Error.+?: ', '', e[1])
                stop(errorMessage)
            }
        },

        replaceSampler = function(...) {
            '
Alias for replaceSamplers method
'
            replaceSamplers(...)
        },
        
        setSamplers = function(..., ind, print = FALSE) {
            '
Sets the ordering of the list of MCMC samplers.

This function also has the side effect of resetting the sampler execution ordering so as to iterate over the specified set of samplers, sequentially, executing each sampler once.

Arguments:

...: Chracter strings or numeric indices.  Character names may be used to specify the node names for samplers to retain.  Numeric indices may be used to specify the indicies for the new list of MCMC samplers, in terms of the current ordered list of samplers.

ind: A numeric vector or character vector.  A numeric vector may be used to specify the indicies for the new list of MCMC samplers, in terms of the current ordered list of samplers.
For example, if the MCMCconf object currently has 3 samplers, then the ordering may be reversed by calling MCMCconf$setSamplers(3:1), or all samplers may be removed by calling MCMCconf$setSamplers(numeric(0)).

Alternatively, a character vector may be used to specify a set of model nodes and/or variables, and the sampler list will modified to only those samplers acting on these target nodes.

As another alternative, a list of samplerConf objects may be used as the argument, in which case this ordered list of samplerConf objects will define the samplers in this MCMC configuration object, completely over-writing the current list of samplers.  No checking is done to ensure the validity of the contents of these samplerConf objects; only that all elements of the list argument are, in fact, samplerConf objects.

print: A logical argument specifying whether to print the new list of samplers (default FALSE).
'   
            if(missing(ind)) {
                ind <- list(...)
                ind <- unname(unlist(ind))
                if(length(ind) == 0)   ind <- numeric(0)
            }
            if(is.list(ind)) {
                if(!all(sapply(ind, class) == 'samplerConf')) stop('item in list argument to setSamplers is not a samplerConf object')
                samplerConfs <<- ind
                return(invisible(NULL))
            }
            if(is.character(ind))   ind <- findSamplersOnNodes(ind)
            if(length(ind) > 0 && max(ind) > length(samplerConfs)) stop('MCMC configuration doesn\'t have that many samplers')
            samplerConfs <<- samplerConfs[ind]
            samplerExecutionOrder <<- seq_along(samplerConfs)
            if(print) printSamplers()
            return(invisible(NULL))
        },

        setSampler = function(...) {
            '
Alias for setSamplers method
'
            setSamplers(...)
        },
        
        printSamplers = function(..., ind, type, displayControlDefaults = FALSE, displayNonScalars = FALSE, displayConjugateDependencies = FALSE, executionOrder = FALSE, byType = FALSE) {
            '
Prints details of the MCMC samplers.

Arguments:

...: Character node or variable names, or numeric indices.  Numeric indices may be used to specify the indices of the samplers to print, or character strings may be used to indicate a set of target nodes and/or variables, for which all samplers acting on these nodes will be printed. For example, printSamplers(\'x\') will print all samplers whose target is model node \'x\', or whose targets are contained (entirely or in part) in the model variable \'x\'.  If omitted, then all samplers are printed.

ind: A numeric vector or character vector.  A numeric vector may be used to specify the indices of the samplers to print, or a character vector may be used to indicate a set of target nodes and/or variables, for which all samplers acting on these nodes will be printed. For example, printSamplers(\'x\') will print all samplers whose target is model node \'x\', or whose targets are contained (entirely or in part) in the model variable \'x\'.  If omitted, then all samplers are printed.

type: a character vector containing sampler type names.  Only samplers with one of these specified types, as printed by this printSamplers method, will be displayed.  Standard regular expression mathing using is also applied.

displayConjugateDependencies: A logical argument, specifying whether to display the dependency lists of conjugate samplers (default FALSE).

displayNonScalars: A logical argument, specifying whether to display the values of non-scalar control list elements (default FALSE).

executionOrder: A logical argument, specifying whether to print the sampler functions in the (possibly modified) order of execution (default FALSE).

byType: A logical argument, specifying whether the nodes being sampled should be printed, sorted and organized according to the type of sampler (the sampling algorithm) which is acting on the nodes (default FALSE).
'
            if(missing(ind)) {
                ind <- list(...)
                ind <- unname(unlist(ind))
                if(length(ind) == 0)   ind <- seq_along(samplerConfs)
            }
            if(is.character(ind))   ind <- findSamplersOnNodes(ind)
            if(length(ind) > 0 && max(ind) > length(samplerConfs)) stop('MCMC configuration doesn\'t have that many samplers')
            if(!missing(type)) {
                if(!is.character(type)) stop('type argument must have type character')
                ## find sampler indices with 'name' matching anything in 'type' argument:
                typeInd <- unique(unname(unlist(lapply(type, grep, x = lapply(samplerConfs, `[[`, 'name')))))
                ind <- intersect(ind, typeInd)
            }
            if(byType) {
                printSamplersByType(ind)
                return(invisible(NULL))
            }
            makeSpaces <- if(length(ind) > 0) newSpacesFunction(max(ind)) else NULL
            if(executionOrder)      ind <- samplerExecutionOrder[samplerExecutionOrder %in% ind]
            for(i in ind) {
                info <- paste0('[', i, '] ', makeSpaces(i), samplerConfs[[i]]$toStr(displayControlDefaults, displayNonScalars, displayConjugateDependencies))
                if(samplerConfs[[i]]$name %in% c('CRP', 'CRP_moreGeneral')) {
                    if(exists('checkConjugacy', samplerConfs[[i]]$control) &&
                       samplerConfs[[i]]$control$checkConjugacy) {
                        conjInfo <- checkCRPconjugacy(model, samplerConfs[[i]]$target)
                        if(is.null(conjInfo)) conjInfo <- "non-conjugate"
                    } else conjInfo <- "non-conjugate"
                    info <- paste0(info, ",  ", conjInfo)
                }
                cat(paste0(info, "\n"))
            }
            if(!executionOrder && !identical(as.numeric(samplerExecutionOrder), as.numeric(seq_along(samplerConfs)))) {
                messageIfVerbose('\n  [Note] Samplers have a modified order of execution.')
                messageIfVerbose('  [Note] To print samplers in the modified order of execution, use printSamplers(executionOrder = TRUE).\n')
            }
            return(invisible(NULL))
        },

        printSamplersByType = function(ind) {
            if(length(ind) == 0) return(invisible(NULL))
            indent <- '  - '
            samplerTypes <- unlist(lapply(ind, function(i) samplerConfs[[i]]$name))
            samplerTypes <- gsub('^conjugate_.+', 'conjugate', samplerTypes)
            uniqueSamplerTypes <- sort(unique(samplerTypes), decreasing = TRUE)
            nodesSortedBySamplerType <- lapply(uniqueSamplerTypes, function(type) sapply(samplerConfs[which(samplerTypes == type)], `[[`, 'target', simplify = FALSE))
            names(nodesSortedBySamplerType) <- uniqueSamplerTypes
            for(i in seq_along(nodesSortedBySamplerType)) {
                theseSampledNodes <- nodesSortedBySamplerType[[i]]
                cat(paste0(names(nodesSortedBySamplerType)[i], ' sampler (', length(theseSampledNodes), ')\n'))
                colonBool <- grepl(':', theseSampledNodes)
                lengthGToneBool <- sapply(theseSampledNodes, length) > 1
                multivariateBool <- colonBool | lengthGToneBool
                univariateList <- theseSampledNodes[!multivariateBool]
                multivariateList <- theseSampledNodes[multivariateBool]
                if(length(univariateList) > 0) {   ## univariate samplers:
                    theseUniVars <- model$getVarNames(nodes = univariateList)
                    uniNodesListByVar <- lapply(theseUniVars, function(var)
                        unlist(univariateList[(univariateList == var) |
                                                  grepl(paste0('^', gsub('\\.','\\\\\\.',var), '\\['), univariateList)]))
                    if(length(unlist(uniNodesListByVar)) != length(univariateList)) stop('something went wrong')
                    for(j in seq_along(uniNodesListByVar)) {
                        theseNodes <- uniNodesListByVar[[j]]
                        isIndexedVector <- grepl("\\[", theseNodes)

                        theseNodesIndexed <- theseNodes[isIndexedVector]
                        if(length(theseNodesIndexed)) {
                            numElements <- length(theseNodesIndexed)
                            sTag <- ifelse(numElements > 1, 's', '')
                            cat(paste0(indent, theseUniVars[j], '[]  (', numElements, ' element', sTag, ')\n'))
                        }
                        theseNodesNotIndexed <- theseNodes[!isIndexedVector]
                        if(length(theseNodesNotIndexed)) {
                            if(length(theseNodesNotIndexed) == 1) cat(paste0(indent, theseNodesNotIndexed))
                            if(length(theseNodesNotIndexed) >  1 && length(unique(theseNodesNotIndexed)) > 1) stop('internal error in printSamplersByType method', call. = FALSE)
                            if(length(theseNodesNotIndexed) >  1) cat(paste0(indent, theseNodesNotIndexed[1], '  (', length(theseNodesNotIndexed), ')'))
                            cat('\n')
                        }
                    }
                }
                if(length(multivariateList) > 0) {   ## multivariate samplers:
                    multiLengthGToneBool <- sapply(multivariateList, length) > 1
                    LGoneNodes <- multivariateList[multiLengthGToneBool]
                    LEoneNodes <- multivariateList[!multiLengthGToneBool]
                    if(length(LEoneNodes) > 0) {
                        theseMultiVars <- model$getVarNames(nodes = LEoneNodes)
                        multiNodesListByVar <- lapply(theseMultiVars, function(var)
                            unlist(LEoneNodes[ grepl(paste0('^', gsub('\\.','\\\\\\.',var), '\\['), LEoneNodes) ]))
                        if(length(unlist(multiNodesListByVar)) != length(LEoneNodes)) stop('something went wrong')
                        for(j in seq_along(multiNodesListByVar)) {
                            theseNodes <- multiNodesListByVar[[j]]
                            numElements <- length(theseNodes)
                            if(numElements > 4) {
                                sTag <- ifelse(numElements>1, 's', '')
                                cat(paste0(indent, theseMultiVars[j], '[]  (', numElements, ' multivariate element', sTag, ')'))
                                cat('\n')
                            } else { theseNodesIndent <- paste0(indent, theseNodes)
                                     cat(paste0(theseNodesIndent, collapse = '\n'), '\n') }
                        }
                    }
                    if(length(LGoneNodes) > 0) {
                        LGoneNodesCompressed <- sapply(LGoneNodes, function(nns) if(length(nns)==1) nns else paste0(nns, collapse = ', '))
                        LGoneNodesCompressedIndent <- paste0(indent, LGoneNodesCompressed)
                        cat(paste0(LGoneNodesCompressedIndent, collapse = '\n'), '\n')
                    }
                }
            }
        },

        getSamplers = function(ind) {
            '
Returns a list of samplerConf objects.

Arguments:

ind: A numeric vector or character vector.  A numeric vector may be used to specify the indices of the samplerConf objects to return, or a character vector may be used to indicate a set of target nodes and/or variables, for which all samplers acting on these nodes will be returned. For example, getSamplers(\'x\') will return all samplerConf objects whose target is model node \'x\', or whose targets are contained (entirely or in part) in the model variable \'x\'.  If omitted, then all samplerConf objects in this MCMC configuration object are returned.
'
            if(missing(ind))        ind <- seq_along(samplerConfs)
            if(is.character(ind))   ind <- findSamplersOnNodes(ind)
            if(length(ind) > 0 && max(ind) > length(samplerConfs)) stop('MCMC configuration doesn\'t have that many samplers')
            return(samplerConfs[ind])
        },

        findSamplersOnNodes = function(nodes) {
            if(length(samplerConfs) == 0) return(integer())
            nodes <- model$expandNodeNames(nodes, returnScalarComponents = TRUE, sort = TRUE)
            samplerConfNodesList <- lapply(samplerConfs, function(sc) sc$targetAsScalar)
            return(unique(unlist(lapply(nodes, function(n) which(unlist(lapply(samplerConfNodesList, function(scn) n %in% scn)))))))
        },

        getSamplerDefinition = function(ind, print = FALSE) {
            '
Returns the nimbleFunction definition of an MCMC sampler.

Arguments:

ind: A numeric vector or character vector.  A numeric vector may be used to specify the index of the sampler definition to return, or a character vector may be used to indicate a target node for which the sampler acting on this nodes will be printed. For example, getSamplerDefinition(\'x[2]\') will return the definition of the sampler whose target is model node \'x[2]\'.  If more than one sampler function is specified, only the first is returned.

Returns a list object, containing the setup function, run function, and additional member methods for the specified nimbleFunction sampler.
'
            if(is.character(ind))   ind <- findSamplersOnNodes(ind)
            if(length(ind) > 1) {
                messageIfVerbose('  [Note] More than one sampler specified, only returning the first.')
                ind <- ind[1]
            }
            if((ind <= 0) || (ind > length(samplerConfs))) stop('Invalid sampler specified')
            if(print) printSamplers(ind)
            def <- getDefinition(samplerConfs[[ind]]$samplerFunction)
            return(def)
        },

        setSamplerExecutionOrder = function(order, print = FALSE) {
            '
Sets the ordering in which sampler functions will execute.

This allows some samplers to be "turned off", or others to execute multiple times in a single MCMC iteration.  The ordering in which samplers execute can also be interleaved.

Arguments:

order: A numeric vector, specifying the ordering in which the sampler functions will execute.  The indices of execution specified in this numeric vector correspond to the enumeration of samplers printed by printSamplers(), or returned by getSamplers().  If this argument is omitted, the sampler execution ordering is reset so as to sequentially execute each sampler once.

print: A logical argument specifying whether to print the current list of samplers in the modified order of execution (default FALSE).
'
            if(missing(order)) order <- seq_along(samplerConfs)
            if(any(order < 1))                    stop('sampler execution ordering must all be positive integers')
            if(any(order != floor(order)))        stop('sampler execution ordering must all be positive integers')
            if(any(order > length(samplerConfs))) stop('sampler execution ordering contains indices larger than the number of sampler functions')
            samplerExecutionOrder <<- order
            if(print) printSamplers(executionOrder = TRUE)
            return(invisible(NULL))
        },

        getSamplerExecutionOrder = function() {
            '
Returns a numeric vector, specifying the ordering of sampler function execution.

The indices of execution specified in this numeric vector correspond to the enumeration of samplers printed by printSamplers(), or returned by getSamplers().
'
            return(samplerExecutionOrder)
        },
        
        addMonitors = function(..., ind = 1, print = TRUE) {
            '
Adds variables to the list of monitors.

Arguments:

...: One or more character vectors of indexed nodes, or variables, which are to be monitored.  These are added onto the current monitors list.

print: A logical argument specifying whether to print all current monitors (default TRUE).

Details: See the initialize() function
            '
            
            if(isMvSamplesReady(ind)){
            	messageIfVerbose('   [Note] Changing monitors, even though an MCMC has been built already. When compiling the MCMC, use resetFunctions = TRUE option.')
            	if(ind == 1)
                    mvSamples1Conf <<- NULL
            	if(ind == 2)
                    mvSamples2Conf <<- NULL
            }

            vars <- list(...)
            if(length(vars) == 1 && is.null(vars[[1]])) {
                if(getNimbleOption('MCMCmonitorAllSampledNodes')) {
                    vars <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
                } else {
                    vars <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE, topOnly = TRUE)
                }
            } else {
                vars <- unlist(vars)
            }
            vars <- unique(removeIndexing(vars))
            nl_checkVarNamesInModel(model, vars)
            if(ind == 1)     monitors  <<- sort(unique(c(monitors,  vars)))
            if(ind == 2)     monitors2 <<- sort(unique(c(monitors2, vars)))
            if(print && nimbleOptions('verbose')) printMonitors()
            return(invisible(NULL))
        },

        addMonitors2 = function(..., print = TRUE) {
            '
Adds variables to the list of monitors2.

Arguments:

...: One or more character vectors of indexed nodes, or variables, which are to be monitored.  These are added onto the current monitors2 list.

print: A logical argument specifying whether to print all current monitors (default TRUE).

Details: See the initialize() function
            '
            addMonitors(..., ind = 2, print = print)
        },

        setMonitors = function(..., ind = 1, print = TRUE) {
            '
Sets new variables to the list of monitors.

Arguments:

...: One or more character vectors of indexed nodes, or variables, which are to be monitored.  These replace the current monitors list.

print: A logical argument specifying whether to print all current monitors (default TRUE).

Details: See the initialize() function
            '
            if(ind == 1)   monitors  <<- character()
            if(ind == 2)   monitors2 <<- character()
            addMonitors(..., ind = ind, print = print)
        },

        setMonitors2 = function(..., print = TRUE) {
            '
Sets new variables to the list of monitors2.

Arguments:

...: One or more character vectors of indexed nodes, or variables, which are to be monitored.  These replace the current monitors2 list.

print: A logical argument specifying whether to print all current monitors (default TRUE).

Details: See the initialize() function
            '
            setMonitors(..., ind = 2, print = print)
        },
        
        resetMonitors = function() {
            '
Resets the current monitors and monitors2 lists to nothing.

Details: See the initialize() function
            '
            monitors  <<- character()
            monitors2 <<- character()
            
            if(isMvSamplesReady(1) || isMvSamplesReady(2)){
            	message('Changing monitors, even though an MCMC has been built already. When compiling the MCMC, use resetFunctions = TRUE option.')
            	mvSamples1Conf <<- NULL
            	mvSamples2Conf <<- NULL
            }

            
            return(invisible(NULL))
        },
        
        printMonitors = function() {
            '
Prints all current monitors and monitors2

Details: See the initialize() function
            '
            if(length(monitors)  > 0)   cat(paste0('thin = ',  thin,  ': ', paste0(monitors,  collapse = ', '), '\n'))
            if(length(monitors2) > 0)   cat(paste0('thin2 = ', thin2, ': ', paste0(monitors2, collapse = ', '), '\n'))
        },

        getMonitors = function() {
            '
Returns a character vector of the current monitors

Details: See the initialize() function
            '
            return(monitors)
        },

        getMonitors2 = function() {
            '
Returns a character vector of the current monitors2

Details: See the initialize() function
            '
            return(monitors2)
        },

        setThin  = function(thin, print = TRUE, ind = 1) {
            '
Sets the value of thin.

Arguments:

thin: The new value for the thinning interval \'thin\'.

print: A logical argument specifying whether to print all current monitors (default TRUE).

Details: See the initialize() function
            '
            if(thin < 1)              stop('cannot use thin < 1', call. = FALSE)
            if(thin != floor(thin))   stop('cannot use non-integer thin', call. = FALSE)
            if(ind == 1)   thin  <<- thin
            if(ind == 2)   thin2 <<- thin
            if(print && nimbleOptions('verbose')) printMonitors()
            return(invisible(NULL))
        },
        setThin2 = function(thin2, print = TRUE) {
            '
Sets the value of thin2.

Arguments:

thin2: The new value for the thinning interval \'thin2\'.

print: A logical argument specifying whether to print all current monitors (default TRUE).

Details: See the initialize() function
            '
            setThin(thin = thin2, print = print, ind = 2)
        },

        getMvSamplesConf  = function(ind = 1){
            
            if(isMvSamplesReady(ind) == TRUE) {
                if(ind == 1) return(mvSamples1Conf)
                return(mvSamples2Conf)
            }
            else{
                makeMvSamplesConf(ind)
                if(ind == 1)
                    output <- mvSamples1Conf
                if(ind == 2)
                    output <- mvSamples2Conf
                return(output)
            }
        },
        
        isMvSamplesReady = function(ind){
            if(ind == 1) return(is(mvSamples1Conf, 'function'))		#Probably really want to give mvConfs there own class...
            if(ind == 2) return(is(mvSamples2Conf, 'function'))
            stop('invalid indicator for isMvSsamplesReady')
        },
        
        makeMvSamplesConf = function(ind){
            modelSymbolObjects = model$getSymbolTable()$getSymbolObjects()
            if(ind == 1) monitorNames = monitors
            if(ind == 2) monitorNames = monitors2
            if(!all(monitorNames %in% names(modelSymbolObjects))) stop('some monitor names are not in the model symbol table; this should never occur')
            thisModelValuesConf = modelValuesConf(symbolTable(symbols = modelSymbolObjects[monitorNames]))
            if(ind == 1) mvSamples1Conf <<- thisModelValuesConf
            if(ind == 2) mvSamples2Conf <<- thisModelValuesConf     	
        },
        
        setUnsampledNodes = function() {
            samplerTargetNodes <- model$expandNodeNames(unlist(lapply(samplerConfs, `[[`, 'target')))
            additionalNodesBeingSampled <- character()
            samplerNames <- sapply(samplerConfs, `[[`, 'name')
            ## special case for posterior_predictive sampler:
            samplerInd <- which(samplerNames == 'posterior_predictive')
            postPredSamplerTargetNodes <- if(length(samplerInd)) model$expandNodeNames(unlist(lapply(samplerConfs[samplerInd], `[[`, 'target'))) else character()
            postPredSamplerNodesSampled <- model$getDependencies(postPredSamplerTargetNodes, stochOnly = TRUE, downstream = TRUE, includePredictive = TRUE)
            additionalNodesBeingSampled <- c(additionalNodesBeingSampled, postPredSamplerNodesSampled)
            postPredSamplerDownstreamNodes <<- setdiff(postPredSamplerNodesSampled, postPredSamplerTargetNodes)
            ## special case for RW_PF and RW_PF_block samplers:
            samplerInd <- which(samplerNames %in% c('RW_PF', 'RW_PF_block'))
            pfSamplerLatentNodes <- if(length(samplerInd)) model$expandNodeNames(unlist(lapply(samplerConfs[samplerInd], function(sconf) sconf$control$latents))) else character()
            additionalNodesBeingSampled <- c(additionalNodesBeingSampled, pfSamplerLatentNodes)
            ##
            allNodesBeingSampled <- unique(c(samplerTargetNodes, additionalNodesBeingSampled))
            unsampledNodes <<- setdiff(model$getNodeNames(stochOnly = TRUE, includeData = FALSE), allNodesBeingSampled)
        },
        
        getUnsampledNodes = function() {
            setUnsampledNodes()
            return(unsampledNodes)
        },
        
        warnUnsampledNodes = function(includeConfGetUnsampledNodes = TRUE) {
            if(length(unsampledNodes)) {
                numUnsampled <- length(unsampledNodes)
                sTag <- if(numUnsampled > 1) 's' else ''
                msg <- paste0('  [Warning] No samplers assigned for ', numUnsampled, ' node', sTag)
                if(includeConfGetUnsampledNodes)   msg <- paste0(msg, ', use conf$getUnsampledNodes() for node name', sTag)
                msg <- paste0(msg, '.')
                messageIfVerbose(msg)
            }
        },

        printComments = function(...) {
            setUnsampledNodes()
            anyComments <-
                length(postPredSamplerDownstreamNodes) ||
                (getNimbleOption('MCMCwarnUnsampledStochasticNodes') && length(unsampledNodes))
            if(anyComments) {
                cat('===== Comments =====\n')
                if(length(postPredSamplerDownstreamNodes))   message('  [Note] Additional downstream predictive nodes are also being sampled by posterior_predictive sampler.')
                if(getNimbleOption('MCMCwarnUnsampledStochasticNodes'))   warnUnsampledNodes(...)
            }
        },

        show = function(...) {
            cat('===== Monitors =====\n')
            printMonitors()
            cat('===== Samplers =====\n')
            if(length(samplerConfs)) printSamplers(byType = TRUE) else cat('(no samplers assigned)\n')
            printComments(...)
        }
    )
)







#' Build the MCMCconf object for construction of an MCMC object
#'
#' Creates a default MCMC configuration for a given model.
#'
#'@param model A NIMBLE model object, created from \code{\link{nimbleModel}}
#'@param nodes An optional character vector, specifying the nodes and/or variables for which samplers should be created.
#'Nodes may be specified in their indexed form, \code{y[1, 3]}.  Alternatively, nodes specified without indexing will be expanded fully, e.g., \code{x} will be expanded to \code{x[1]}, \code{x[2]}, etc.
#'If missing, the default value is all non-data stochastic nodes.
#'If NULL, then no samplers are added.
#'@param control An optional list of control arguments to sampler functions.  If a control list is provided, the elements will be provided to all sampler functions which utilize the named elements given.
#'For example, the standard Metropolis-Hastings random walk sampler (\link{sampler_RW}) utilizes control list elements \code{adaptive}, \code{adaptInterval}, and \code{scale}.
#' (Internally it also uses \code{targetNode}, but this should not generally be provided as a control list element).
#'The default values for control list arguments for samplers (if not otherwise provided as an argument to configureMCMC() ) are in the setup code of the sampling algorithms.
#'@param monitors A character vector of node names or variable names, to record during MCMC sampling.
#'This set of monitors will be recorded with thinning interval \code{thin}, and the samples will be stored into the \code{mvSamples} object.
#'The default value is all top-level stochastic nodes of the model -- those having no stochastic parent nodes.
#'@param monitors2 A character vector of node names or variable names, to record during MCMC sampling.
#'This set of monitors will be recorded with thinning interval \code{thin2}, and the samples will be stored into the \code{mvSamples2} object.
#'The default value is an empty character vector, i.e. no values will be recorded.
#'@param thin The thinning interval for \code{monitors}.  Default value is one.
#'@param thin2 The thinning interval for \code{monitors2}.  Default value is one.
#'@param useConjugacy A logical argument, with default value TRUE.  If specified as FALSE, then no conjugate samplers will be used, even when a node is determined to be in a conjugate relationship.
#'@param onlyRW A logical argument, with default value FALSE.  If specified as TRUE, then Metropolis-Hastings random walk samplers (\link{sampler_RW}) will be assigned for all non-terminal continuous-valued nodes nodes. Discrete-valued nodes are assigned a slice sampler (\link{sampler_slice}), and terminal nodes are assigned a posterior_predictive sampler (\link{sampler_posterior_predictive}).
#'@param onlySlice A logical argument, with default value FALSE.  If specified as TRUE, then a slice sampler is assigned for all non-terminal nodes. Terminal nodes are still assigned a posterior_predictive sampler.
#'@param multivariateNodesAsScalars A logical argument, with default value FALSE.  If specified as TRUE, then non-terminal multivariate stochastic nodes will have scalar samplers assigned to each of the scalar components of the multivariate node.  The default value of FALSE results in a single block sampler assigned to the entire multivariate node.  Note, multivariate nodes appearing in conjugate relationships will be assigned the corresponding conjugate sampler (provided \code{useConjugacy == TRUE}), regardless of the value of this argument.
#'@param enableWAIC A logical argument, specifying whether to enable WAIC calculations for the resulting MCMC algorithm.  Defaults to the value of \code{nimbleOptions('MCMCenableWAIC')}, which in turn defaults to FALSE.  Setting \code{nimbleOptions('enableWAIC' = TRUE)} will ensure that WAIC is enabled for all calls to \code{\link{configureMCMC}} and \code{\link{buildMCMC}}.
#'@param controlWAIC A named list of inputs that control the behavior of the WAIC calculation. See \code{help(waic)}.
#'@param print A logical argument, specifying whether to print the ordered list of default samplers.
#'@param autoBlock A logical argument specifying whether to use an automated blocking procedure to determine blocks of model nodes for joint sampling.  If TRUE, an MCMC configuration object will be created and returned corresponding to the results of the automated parameter blocking.  Default value is FALSE.
#'@param oldConf An optional MCMCconf object to modify rather than creating a new MCMCconf from scratch
#'@param ... Additional named control list elements for default samplers, or additional arguments to be passed to the \code{\link{autoBlock}} function when \code{autoBlock = TRUE}
#'@author Daniel Turek
#'@export 
#'@details See \code{\link{MCMCconf}} for details on how to manipulate the \code{MCMCconf} object
#'@seealso \code{\link{buildMCMC}} \code{\link{runMCMC}} \code{\link{nimbleMCMC}}
configureMCMC <- function(model, nodes, control = list(), 
                          monitors, thin = 1, monitors2 = character(), thin2 = 1,
                          useConjugacy = getNimbleOption('MCMCuseConjugacy'),
                          onlyRW = FALSE, onlySlice = FALSE,
                          multivariateNodesAsScalars = getNimbleOption('MCMCmultivariateNodesAsScalars'),
                          enableWAIC = getNimbleOption('MCMCenableWAIC'), controlWAIC = list(),
                          print = getNimbleOption('verbose'),
                          autoBlock = FALSE, oldConf,
                          ## samplerAssignmentRules system deprecated Nov 2020 -DT
                          ##rules = getNimbleOption('MCMCdefaultSamplerAssignmentRules'),
                          ...) {
    
    ## samplerAssignmentRules system deprecated Nov 2020 -DT
    ##if(!inherits(rules, 'samplerAssignmentRules')) stop('rules argument must be a samplerAssignmentRules object')

    if(!missing(oldConf)){
        if(!is(oldConf, 'MCMCconf'))
            stop('oldConf must be an MCMCconf object, as built by the configureMCMC function')
        return(makeNewConfFromOldConf(oldConf))	
    }
    
    if(missing(model))        stop('Either oldConf or model must be supplied')
    if(missing(monitors))     monitors <- NULL

    if(autoBlock) return(autoBlock(model, ...)$conf)

    thisConf <- MCMCconf(model = model, nodes = nodes, control = control, ##rules = rules,
                         monitors = monitors, thin = thin, monitors2 = monitors2, thin2 = thin2,
                         useConjugacy = useConjugacy,
                         onlyRW = onlyRW, onlySlice = onlySlice,
                         multivariateNodesAsScalars = multivariateNodesAsScalars,
                         enableWAIC = enableWAIC, controlWAIC = controlWAIC,
                         print = print, ...)
    return(invisible(thisConf))
}



# This is function which builds a new MCMCconf from an old MCMCconf
# This is required to be able to a new C-based MCMC without recompiling
makeNewConfFromOldConf <- function(oldMCMCconf){
    newMCMCconf <- configureMCMC(oldMCMCconf$model, nodes = NULL, print = FALSE)
    newMCMCconf$monitors <- oldMCMCconf$monitors
    newMCMCconf$monitors2 <- oldMCMCconf$monitors2
    newMCMCconf$thin <- oldMCMCconf$thin
    newMCMCconf$thin2 <- oldMCMCconf$thin2
    newMCMCconf$samplerConfs <- oldMCMCconf$samplerConfs
    newMCMCconf$samplerExecutionOrder <- oldMCMCconf$samplerExecutionOrder
    newMCMCconf$controlDefaults <- oldMCMCconf$controlDefaults
    ##newMCMCconf$namedSamplerLabelMaker <- oldMCMCconf$namedSamplerLabelMaker  ## usage long since deprecated (Dec 2020)
    newMCMCconf$mvSamples1Conf <- oldMCMCconf$mvSamples1Conf
    newMCMCconf$mvSamples2Conf <- oldMCMCconf$mvSamples2Conf
    return(newMCMCconf)	
}


newSpacesFunction <- function(m) {
    log10max <- floor(log10(m))
    function(i) paste0(rep(' ', log10max-floor(log10(i))), collapse = '')
}



