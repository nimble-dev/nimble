
samplerConf <- setRefClass(
    Class = 'samplerConf',
    fields = list(
        name            = 'ANY',
        samplerFunction = 'ANY',
        target          = 'ANY',
        control         = 'ANY',
        targetAsScalar  = 'ANY'
    ),
    methods = list(
        initialize = function(name, samplerFunction, target, control, model) {
            name <<- name
            samplerFunction <<- samplerFunction
            target <<- target
            control <<- control
            targetAsScalar <<- model$expandNodeNames(target, returnScalarComponents = TRUE)
        },
        setName = function(name) name <<- name,
        setSamplerFunction = function(fun) samplerFunction <<- fun,
        setTarget = function(target, model) {
            target <<- target
            targetAsScalar <<- model$expandNodeNames(target, returnScalarComponents = TRUE)
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
#' @aliases MCMCconf addSampler removeSamplers setSamplers printSamplers getSamplers setSamplerExecutionOrder getSamplerExecutionOrder addMonitors addMonitors2 resetMonitors getMonitors getMonitors2 printMonitors setThin setThin2
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
        samplerConfs        = 'ANY',
        samplerExecutionOrder = 'ANY',
        controlDefaults     = 'ANY',
        namedSamplerLabelMaker = 'ANY',
        mvSamples1Conf      = 'ANY',
        mvSamples2Conf      = 'ANY'
    ),
    
    methods = list(
        
        initialize = function(model, nodes, control = list(), rules,
            monitors,                thin  = 1,
            monitors2 = character(), thin2 = 1,
            useConjugacy = TRUE,
            onlyRW = FALSE,
            onlySlice = FALSE,
            multivariateNodesAsScalars = FALSE,
            enableWAIC = nimbleOptions('enableWAIC'),
            warnNoSamplerAssigned = TRUE,
            print = FALSE, ...) {
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

enableWAIC: A logical argument, specifying whether to enable WAIC calculations for the resulting MCMC algorithm.  Defaults to the value of nimbleOptions(\'enableWAIC\'), which in turn defaults to FALSE.  Setting nimbleOptions(\'enableWAIC\' = TRUE) will ensure that WAIC is enabled for all calls to configureMCMC and buildMCMC.

warnNoSamplerAssigned: A logical argument specifying whether to issue a warning when no sampler is assigned to a node, meaning there is no matching sampler assignment rule. Default is TRUE.

print: A logical argument specifying whether to print the ordered list of default samplers.  Default is FALSE.

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
            thin  <<- thin
            thin2 <<- thin2
            enableWAIC <<- enableWAIC
            samplerConfs <<- list()
            samplerExecutionOrder <<- numeric()
            controlDefaults <<- list(...)
            namedSamplerLabelMaker <<- labelFunctionCreator('namedSampler')
            for(i in seq_along(control))     controlDefaults[[names(control)[i]]] <<- control[[i]]
            if(identical(nodes, character())) { nodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
                                            } else             { if(is.null(nodes) || length(nodes)==0)     nodes <- character(0)
                                                                 nl_checkVarNamesInModel(model, removeIndexing(nodes))
                                                                 nodes <- model$expandNodeNames(nodes)            }
            
            nodes <- model$topologicallySortNodes(nodes)   ## topological sort
            if(!(all(model$isStoch(nodes)))) { stop('assigning samplers to non-stochastic nodes: ', paste0(nodes[!model$isStoch(nodes)], collapse=', ')) }    ## ensure all target node(s) are stochastic

            if(getNimbleOption('MCMCuseSamplerAssignmentRules')) {
                ## use new system of samplerAssignmentRules
                isEndNodeAll <- model$isEndNode(nodes)
                isMultivariateAll <- model$isMultivariate(nodes)
                isDiscreteAll <- model$isDiscrete(nodes)
                isBinaryAll <- model$isBinary(nodes)
                nodeDistributionsAll <- model$getDistribution(nodes)
                if(useConjugacy) {
                    conjugacyResultsAll <- model$checkConjugacy(nodes)
                    isConjugateAll <- nodes %in% names(conjugacyResultsAll)
                }
                
                ruleSelectFunction <- function() {}
                body(ruleSelectFunction) <- rules$makeRuleSelectionCodeBlock()
                
                for(i in seq_along(nodes)) {
                    node <- nodes[i]
                    isEndNode <- isEndNodeAll[i]
                    isMultivariate <- isMultivariateAll[i]
                    isDiscrete <- isDiscreteAll[i]
                    isBinary <- isBinaryAll[i]
                    nodeDistribution <- nodeDistributionsAll[i]
                    if(useConjugacy) isConjugate <- isConjugateAll[i]
                    ruleSelectFunction()
                }
            } else {
                ## use old (static) system for assigning default samplers
                isEndNode <- model$isEndNode(nodes)
                if(useConjugacy) conjugacyResultsAll <- model$checkConjugacy(nodes)
                
                for(i in seq_along(nodes)) {
                    node <- nodes[i]
                    discrete <- model$isDiscrete(node)
                    binary <- model$isBinary(node)
                    nodeDist <- model$getDistribution(node)
                    nodeScalarComponents <- model$expandNodeNames(node, returnScalarComponents = TRUE)
                    nodeLength <- length(nodeScalarComponents)
                    
                    ## if node has 0 stochastic dependents, assign 'posterior_predictive' sampler (e.g. for predictive nodes)
                    if(isEndNode[i]) { addSampler(target = node, type = 'posterior_predictive');     next }
                    
                    ## for multivariate nodes, either add a conjugate sampler, RW_multinomial, or RW_block sampler
                    if(nodeLength > 1) {
                        if(useConjugacy) {
                            conjugacyResult <- conjugacyResultsAll[[node]]
                            if(!is.null(conjugacyResult)) {
                                addConjugateSampler(conjugacyResult = conjugacyResult,
                                                    dynamicallyIndexed = model$modelDef$varInfo[[model$getVarNames(nodes=node)]]$anyDynamicallyIndexed);     next }
                        }
                        if(nodeDist == 'dmulti')       { addSampler(target = node, type = 'RW_multinomial');     next }
                        if(nodeDist == 'ddirch')       { addSampler(target = node, type = 'RW_dirichlet');       next }
                        if(nodeDist == 'dwish')        { addSampler(target = node, type = 'RW_wishart');         next }
                        if(nodeDist == 'dinvwish')     { addSampler(target = node, type = 'RW_wishart');         next }
                        if(nodeDist == 'dcar_normal')  { addSampler(target = node, type = 'CAR_normal');         next }
                        if(nodeDist == 'dcar_proper')  { addSampler(target = node, type = 'CAR_proper');         next }
                        if(nodeDist == 'dCRP')         { addSampler(target = node, type = 'CRP');                next }
                        if(multivariateNodesAsScalars) {
                            for(scalarNode in nodeScalarComponents) {
                                if(onlySlice) addSampler(target = scalarNode, type = 'slice')
                                else          addSampler(target = scalarNode, type = 'RW')    };     next }
                        addSampler(target = node, type = 'RW_block', silent = TRUE);     next }
                    
                    if(onlyRW && !discrete)   { addSampler(target = node, type = 'RW'   );     next }
                    if(onlySlice)             { addSampler(target = node, type = 'slice');     next }
                    
                    ## if node passes checkConjugacy(), assign 'conjugate_dxxx' sampler
                    if(useConjugacy) {
                        conjugacyResult <- conjugacyResultsAll[[node]]
                        if(!is.null(conjugacyResult)) {
                            addConjugateSampler(conjugacyResult = conjugacyResult,
                                                dynamicallyIndexed = model$modelDef$varInfo[[model$getVarNames(nodes=node)]]$anyDynamicallyIndexed);     next }
                    }
                    
                    ## if node is discrete 0/1 (binary), assign 'binary' sampler
                    if(binary) { addSampler(target = node, type = 'binary');     next }
                    
                    ## for categorical nodes, assign a 'categorical' sampler
                    if(nodeDist == 'dcat') { addSampler(target = node, type = 'categorical');     next }
                    
                    ## if node distribution is discrete, assign 'slice' sampler
                    if(discrete) { addSampler(target = node, type = 'slice');     next }
                    
                    ## if node distribution is dgamma and its dependency is dCRP, assign 'augmented_BetaGamma' sampler
                    if(nodeDist == 'dgamma'){
                      depNode <- model$getDependencies(node, self=FALSE)
                      depNodeDist <- model$getDistribution(depNode)
                      if(length(depNodeDist) == 1 && depNodeDist == 'dCRP'){
                        addSampler(target = node, type = 'CRP_concentration')
                        next
                      }
                    }
                    
                    ## default: 'RW' sampler
                    addSampler(target = node, type = 'RW');     next
                }
            }
            
            if(print)   printSamplers()
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
        
        addSampler = function(target, type = 'RW', control = list(), print = FALSE, name, silent = FALSE, ...) {
            '
Adds a sampler to the list of samplers contained in the MCMCconf object.

Arguments:

target: The target node or nodes to be sampled.  This may be specified as a character vector of model node and/or variable names.  This argument is required.

type: The type of sampler to add, specified as either a character string or a nimbleFunction object.  If the character argument type=\'newSamplerType\', then either newSamplerType or sampler_newSamplertype must correspond to a nimbleFunction (i.e. a function returned by nimbleFunction, not a specialized nimbleFunction).  Alternatively, the type argument may be provided as a nimbleFunction itself rather than its name.  In that case, the \'name\' argument may also be supplied to provide a meaningful name for this sampler.  The default value is \'RW\' which specifies scalar adaptive Metropolis-Hastings sampling with a normal proposal distribution. This default will result in an error if \'target\' specifies more than one target node.

control: A list of control arguments specific to the sampler function. These will override those specified in the control list argument to configureMCMC().

print: Logical argument, specifying whether to print the details of the newly added sampler, as well as its position in the list of MCMC samplers.

name: Optional character string name for the sampler, which is used by the printSamplers method.  If \'name\' is not provided, the \'type\' argument is used to generate the sampler name.

silent: Logical argument, specifying whether to print warning messages when assigning samplers.

...: Additional named arguments passed through ... will be used as additional control list elements.

Details: A single instance of the newly configured sampler is added to the end of the list of samplers for this MCMCconf object.

Invisibly returns a list of the current sampler configurations, which are samplerConf reference class objects.
'

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
                    message('Note: Assigning an RW_block sampler to nodes with very different scales can result in low MCMC efficiency.  If all nodes assigned to RW_block are not on a similar scale, we recommend providing an informed value for the \"propCov\" control list argument, or using the AFSS sampler instead.')
                }
                if(exists(type) && is.nfGenerator(eval(as.name(type)))) {   ## try to find sampler function 'type'
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

            ##libraryTag <- if(nameProvided) namedSamplerLabelMaker() else thisSamplerName   ## unique tag for each 'named' sampler, internal use only
            ##if(is.null(controlNamesLibrary[[libraryTag]]))   controlNamesLibrary[[libraryTag]] <<- mcmc_findControlListNamesInCode(samplerFunction)   ## populate control names library
            ##requiredControlNames <- controlNamesLibrary[[libraryTag]]
            controlArgs <- c(control, list(...))
            thisControlList <- mcmc_generateControlListArgument(control=controlArgs, controlDefaults=controlDefaults)  ## should name arguments
            
            newSamplerInd <- length(samplerConfs) + 1
            samplerConfs[[newSamplerInd]] <<- samplerConf(name=thisSamplerName, samplerFunction=samplerFunction, target=target, control=thisControlList, model=model)
            samplerExecutionOrder <<- c(samplerExecutionOrder, newSamplerInd)
            
            if(print) printSamplers(newSamplerInd)
            return(invisible(samplerConfs))
        },
        
        removeSamplers = function(ind, print = FALSE) {
            '
Removes one or more samplers from an MCMCconf object.

Arguments:

This function also has the side effect of resetting the sampler execution ordering so as to iterate over the remaining set of samplers, sequentially, executing each sampler once.

ind: A numeric vector or character vector specifying the samplers to remove.  A numeric vector may specify the indices of the samplers to be removed.  Alternatively, a character vector may be used to specify a set of model nodes and/or variables, and all samplers whose \'target\' is among these nodes will be removed.  If omitted, then all samplers are removed.

print: A logical argument specifying whether to print the current list of samplers once the removal has been done (default FALSE).
'      
            if(missing(ind))        ind <- seq_along(samplerConfs)
            if(is.character(ind))   ind <- findSamplersOnNodes(ind)
            if(length(ind) > 0 && max(ind) > length(samplerConfs)) stop('MCMC configuration doesn\'t have that many samplers')
            samplerConfs[ind] <<- NULL
            samplerExecutionOrder <<- seq_along(samplerConfs)
            if(print) printSamplers()
            return(invisible(NULL))
        },

        removeSampler = function(...){
            '
Alias for removeSamplers method
'
            removeSamplers(...)
        },
        
        setSamplers = function(ind, print = FALSE) {
            '
Sets the ordering of the list of MCMC samplers.

This function also has the side effect of resetting the sampler execution ordering so as to iterate over the specified set of samplers, sequentially, executing each sampler once.

Arguments:

ind: A numeric vector or character vector.  A numeric vector may be used to specify the indicies for the new list of MCMC samplers, in terms of the current ordered list of samplers.
For example, if the MCMCconf object currently has 3 samplers, then the ordering may be reversed by calling MCMCconf$setSamplers(3:1), or all samplers may be removed by calling MCMCconf$setSamplers(numeric(0)).

Alternatively, a character vector may be used to specify a set of model nodes and/or variables, and the sampler list will modified to only those samplers acting on these target nodes.

As another alternative, a list of samplerConf objects may be used as the argument, in which case this ordered list of samplerConf objects will define the samplers in this MCMC configuration object, completely over-writing the current list of samplers.  No checking is done to ensure the validity of the contents of these samplerConf objects; only that all elements of the list argument are, in fact, samplerConf objects.

print: A logical argument specifying whether to print the new list of samplers (default FALSE).
'   
            if(missing(ind))        ind <- numeric(0)
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
        
        printSamplers = function(ind, type, displayControlDefaults = FALSE, displayNonScalars = FALSE, displayConjugateDependencies = FALSE, executionOrder = FALSE) {
            '
Prints details of the MCMC samplers.

Arguments:

ind: A numeric vector or character vector.  A numeric vector may be used to specify the indices of the samplers to print, or a character vector may be used to indicate a set of target nodes and/or variables, for which all samplers acting on these nodes will be printed. For example, printSamplers(\'x\') will print all samplers whose target is model node \'x\', or whose targets are contained (entirely or in part) in the model variable \'x\'.  If omitted, then all samplers are printed.

type: a character vector containing sampler type names.  Only samplers with one of these specified types, as printed by this printSamplers method, will be displayed.  Standard regular expression mathing using is also applied.

displayConjugateDependencies: A logical argument, specifying whether to display the dependency lists of conjugate samplers (default FALSE).

displayNonScalars: A logical argument, specifying whether to display the values of non-scalar control list elements (default FALSE).

executionOrder: A logical argument, specifying whether to print the sampler functions in the (possibly modified) order of execution (default FALSE).
'
            if(missing(ind))        ind <- seq_along(samplerConfs)
            if(is.character(ind))   ind <- findSamplersOnNodes(ind)
            if(length(ind) > 0 && max(ind) > length(samplerConfs)) stop('MCMC configuration doesn\'t have that many samplers')
            if(!missing(type)) {
                if(!is.character(type)) stop('type argument must have type character')
                ## find sampler indices with 'name' matching anything in 'type' argument:
                typeInd <- unique(unname(unlist(lapply(type, grep, x = lapply(conf$samplerConfs, `[[`, 'name')))))
                ind <- intersect(ind, typeInd)
            }
            makeSpaces <- if(length(ind) > 0) newSpacesFunction(max(ind)) else NULL
            if(executionOrder)      ind <- samplerExecutionOrder[samplerExecutionOrder %in% ind]
            for(i in ind)
                cat(paste0('[', i, '] ', makeSpaces(i), samplerConfs[[i]]$toStr(displayControlDefaults, displayNonScalars, displayConjugateDependencies), '\n'))
            if(!executionOrder && !identical(as.numeric(samplerExecutionOrder), as.numeric(seq_along(samplerConfs)))) {
                cat('These sampler functions have a modified order of execution.\n')
                cat('To print samplers in the modified order of execution, use printSamplers(executionOrder = TRUE).\n')
            }
            return(invisible(NULL))
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
            nodes <- model$expandNodeNames(nodes, returnScalarComponents = TRUE)
            which(unlist(lapply(samplerConfs, function(ss) any(nodes %in% ss$targetAsScalar))))
        },

        getSamplerDefinition = function(ind) {
            '
Returns the nimbleFunction definition of an MCMC sampler.

Arguments:

ind: A numeric vector or character vector.  A numeric vector may be used to specify the index of the sampler definition to return, or a character vector may be used to indicate a target node for which the sampler acting on this nodes will be printed. For example, getSamplerDefinition(\'x[2]\') will return the definition of the sampler whose target is model node \'x[2]\'.  If more than one sampler function is specified, only the first is returned.

Returns a list object, containing the setup function, run function, and additional member methods for the specified nimbleFunction sampler.
'
            if(is.character(ind))   ind <- findSamplersOnNodes(ind)
            if(length(ind) > 1) {
                message('More than one sampler specified, only returning the first')
                ind <- ind[1]
            }
            if((ind <= 0) || (ind > length(samplerConfs))) stop('Invalid sampler specified')
            printSamplers(ind)
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
            	cat('Changing monitors, even though an MCMC has been built already. When compiling the MCMC, use resetFunctions = TRUE option\n')
            	if(ind == 1)
                    mvSamples1Conf <<- NULL
            	if(ind == 2)
                    mvSamples2Conf <<- NULL
            }

            vars <- list(...)
            if(length(vars) == 1 && is.null(vars[[1]])) {
                vars <- model$getNodeNames(topOnly = TRUE, stochOnly = TRUE)
            } else {
                vars <- unlist(vars)
            }
            vars <- unique(removeIndexing(vars))
            nl_checkVarNamesInModel(model, vars)
            if(ind == 1)     monitors  <<- unique(c(monitors,  vars))
            if(ind == 2)     monitors2 <<- unique(c(monitors2, vars))
            if(print) printMonitors()
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
        
        resetMonitors = function() {
            '
Resets the current monitors and monitors2 lists to nothing.

Details: See the initialize() function
            '
            monitors  <<- character()
            monitors2 <<- character()
            
            if(isMvSamplesReady(1) || isMvSamplesReady(2)){
            	cat('Changing monitors, even though an MCMC has been built already. When compiling the MCMC, use resetFunctions = TRUE option\n')
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

        setThin  = function(thin, print = TRUE) {
            '
Sets the value of thin.

Arguments:

thin: The new value for the thinning interval \'thin\'.

print: A logical argument specifying whether to print all current monitors (default TRUE).

Details: See the initialize() function
            '
            thin <<- thin
            if(print) printMonitors()
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
            thin2 <<- thin2
            if(print) printMonitors()
            return(invisible(NULL))
        },

        setEnableWAIC = function(waic = TRUE) {
            '
Sets the value of enableWAIC.

Arguments:

waic: A logical argument, indicating whether to enable WAIC calculations in the resulting MCMC algorithm (default TRUE).
'
            enableWAIC <<- as.logical(waic)
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

        show = function() {
            cat('MCMC configuration object\n')
        }
    )
)

checkCRPconjugacy <- function(model, target) {
    ## Checks if can use conjugacy in drawing new components for dCRP node updating.
    ## Should find conjugacy if there are no deterministic nodes between the "observations" and
    ## the cluster parameters or one intermediate node.
    ## For now, no conjugacy when multiple parameters in the observation distribution are
    ## being clustered.
    
    conjugate <- FALSE 
    
    targetElementExample <- model$expandNodeNames(target, returnScalarComponents=TRUE)[1]
    
    ## Find the cluster variables ("tilde" variables)
    
    tildeVars <- NULL
    idx <- 1
    ## Check that cluster ID appears in dependency expression and that
    ## potential tilde variable is a vector simply indexed by the cluster ID.
    ## For now we do not detect conjugacy for cases like mu[xi[i], 2],
    ## as these are not handled by our CRP conjugate samplers.
    deps <- model$getDependencies(targetElementExample, self=FALSE)
    for(i in seq_along(deps)) {
      expr <- cc_getNodesInExpr(model$getValueExpr(deps[i]))
      expr <- parse(text = expr)[[1]]
      if(is.call(expr) && length(expr) == 3 && expr[[1]] == '[' &&
         expr[[3]] == targetElementExample) {
          tildeVars[idx] <- deparse(expr[[2]])
          idx <- idx + 1
      }
      ## Allow for exact model structure where conjugate sampler works for ddirch/multi.:
      ## p[xi[i], 1:k]
      if(is.call(expr) && length(expr) == 4 && expr[[1]] == '[' &&
         length(expr[[4]]) == 3 && expr[[4]][[2]] == 1 && 
         expr[[3]] == targetElementExample) {
          tildeVars[idx] <- deparse(expr[[2]])
          idx <- idx + 1
      }
    }
    ## Note that models with thetaTilde[n-xi[i]+1] or thetaTilde[xi[i]+1]
    ## won't have xi[i] detected above and so no conjugacy will be found,
    ## which is what we want because conjugate samplers assume xi[i] is the thetaTilde element being used.

    ## New checking for conjugacy using tilde variables: check conjugacy for one tilde node and
    ## then make sure all tilde nodes and all of their dependent nodes are exchangeable
    if(length(tildeVars) == 1){  ## for now avoid case of mixing over multiple parameters
        clusterNodes <- model$expandNodeNames(tildeVars[1])  # e.g., 'thetatilde[1]',...,
        ## Avoid non-nodes from truncated clustering,
        ## e.g., avoid 'thetaTilde[3:10]' if only have 2 thetaTilde nodes but 10 obs.
        clusterNodes <- clusterNodes[clusterNodes %in% model$getNodeNames(stochOnly = TRUE, includeData = FALSE)]
        conjugacy <- model$checkConjugacy(clusterNodes[1], restrictLink = 'identity')
        if(length(conjugacy)) {
            if(length(unique(model$getDeclID(clusterNodes))) == 1) { ## make sure all tilde nodes from same declaration (i.e., exchangeable)
                depNodes <- model$getDependencies(clusterNodes[1], stochOnly = TRUE, self=FALSE)
                if(length(unique(model$getDeclID(depNodes))) == 1) { ## make sure all dependent nodes from same declaration (i.e., exchangeable)

                    conjugacyType <- paste0(conjugacy[[1]]$type, '_', sub('dep_', '', names(conjugacy[[1]]$control))) 
                    conjugate <- TRUE
                }
            }
        }
        ## Check that prior for tilde nodes are truly exchangeable.
        valueExprs <- sapply(clusterNodes, function(x) model$getValueExpr(x))
        names(valueExprs) <- NULL
        if(length(unique(valueExprs)) != 1)
            conjugate <- FALSE
    }
    
    if(conjugate) return(conjugacyType) else return(NULL)
}


rule <- setRefClass(
    Class = 'rule',
    fields = list(
        condition = 'ANY',
        sampler = 'ANY',
        name = 'character'
    ),
    methods = list(
        initialize = function(condition, sampler, name, nameProvided) {
            condition <<- condition
            sampler <<- sampler      ## must be of class: character, function, or call
            name <<- name
        },
        getCondition = function()   return(condition),
        getSampler = function()     return(sampler),
        getName = function()        return(name),
        show = function() {
            if(is.character(sampler) || is.function(sampler)) sampPrint <- as.name(name) else sampPrint <- sampler    ## ifelse() fails here, seems to be R bug w/ quoted expressions
            cat('condition: ', paste0(deparse(condition),sep='\n'), sep='')
            cat('sampler:   ', paste0(deparse(sampPrint),sep='\n'), sep='')
        }
    )
)


addRuleToCodeBlock <- function(oldCode, rule) {
    condition <- rule$getCondition()
    sampler <- rule$getSampler()
    name <- rule$getName()
    if(is.call(sampler)) {
        substitute(
            if(CONDITION) { SAMPLER } else OLDCODE,
            list(CONDITION = condition,
                 SAMPLER = sampler,
                 NAME = name,
                 OLDCODE = oldCode))
    } else {
        substitute(
            if(CONDITION) { addSampler(target = node, type = SAMPLER, name = NAME, silent = TRUE) } else OLDCODE,
            list(CONDITION = condition,
                 SAMPLER = sampler,
                 NAME = name,
                 OLDCODE = oldCode))
    }
}

#' Class \code{samplerAssignmentRules}
#' @aliases samplerAssignmentRules addRule reorder printRules
#' @export
#' @description
#' Objects of this class specify an ordered set of rules for assigning MCMC sampling algorithms to the stochastic nodes in a BUGS model.
#' This feature can be enabled by setting the NIMBLE option \code{MCMCuseSamplerAssignmentRules} to \code{TRUE}.
#' The rules can be modified to alter under what circumstances various samplers are assigned, and with what precedence.
#' When assigning samplers to each stochastic node, the set of rules is traversed beginning with the first, until a matching rule is found.
#' When a matching rule is found, the sampler specified by that rule is assigned (or general code for sampler assignment is executed),
#' and the assignment process proceeds to the next stochastic node.  That is, a maximum of one rule can be invoked for each stochastic node.
#' If no matching rule is found, an (optional) warning is issued and no sampler is assigned.
#' Objects of this class may be passed using the \code{rules} argument to \code{\link{configureMCMC}} to customize the sampler assignment process.
#' See documentation below for method \code{initialize()} for details of creating a samplerAssignmentRules object, 
#' and methods \code{addRule()} and \code{reorder()} for adding and modifying the sampler assignment rules.
#' The default behaviour of \code{configureMCMC} can be modified by setting the nimble option \'MCMCsamplerAssignmentRules\' to a customized samplerAssignmentRules object.
#' The default behaviour of \code{configureMCMC} can be restored using \code{nimbleOptions(MCMCdefaultSamplerAssignmentRules = samplerAssignmentRules())}.
#' @author Daniel Turek
#' @seealso \code{\link{configureMCMC}}
#' @examples
#' \dontrun{
#' ## enable the use of samplerAssignmentRules:
#' nimbleOptions(MCMCuseSamplerAssignmentRules = TRUE)
#' 
#' ## omitting empty=TRUE creates a copy of nimble's default rules
#' my_rules <- samplerAssignmentRules(empty = TRUE)
#' 
#' my_rules$addRule(quote(model$isEndNode(node)), "posterior_predictive")
#' my_rules$addRule(quote(model$isDiscrete(node)), "my_new_discrete_sampler")
#' my_rules$addRule(TRUE, "RW")   ## default catch-all sampler assignment
#'
#' ## print the ordered set of sampler assignment rules
#' my_rules$printRules()
#'
#' ## assign samplers according to my_rules object
#' conf <- configureMCMC(Rmodel, rules = my_rules)
#' conf$printSamplers()
#'
#' ## view the current (default) assignment rules used by configureMCMC()
#' nimbleOptions(MCMCdefaultSamplerAssignmentRules)
#'
#' ## change default behaviour of configureMCMC() to use my_rules
#' nimbleOptions(MCMCdefaultSamplerAssignmentRules = my_rules)
#' 
#' ## reset configureMCMC() to use default rules
#' nimbleOptions(MCMCdefaultSamplerAssignmentRules = samplerAssignmentRules())
#' }
samplerAssignmentRules <- setRefClass(
    Class = 'samplerAssignmentRules',
    fields = list(
        ruleList = 'list'       ## list of rule objects
    ),
    methods = list(
        initialize = function(empty = FALSE, print = FALSE) {
            '
Creates a new samplerAssignmentRules object, which is a container for an ordered set of rules for MCMC sampler assignments.  Objects of this class may be passed using the \'rules\' argument to configureMCMC(), to customize the process of assigning samplers to stochastic model nodes.  By default, new samplerAssignmentRules objects are initialized having an exact copy of the default sampler assignment rules used by NIMBLE, and can thereafter be modified using the addRule() and reorder() methods.

Arguments:

empty: Logical argument (default = FALSE).  If TRUE, then a new samplerAssignmentRules object is created containing no rules.  The default behaviour creates new objects containing an exact copy of the default sampler assignment rules used by NIMBLE.

print: Logical argument specifying whether to print the ordered list of sampler assignment rules (default FALSE).
'
            ruleList <<- list()
            if(!empty) addDefaultSamplerAssignmentRules()
            if(print) printRules()
        },
        makeRuleSelectionCodeBlock = function() {
            code <- quote(if(warnNoSamplerAssigned) {    ## no matching rule was found
                warning(paste0('No matching rule found, and no sampler assigned to node: ', node))
            })
            for(i in rev(seq_along(ruleList)))   ## important to add rules in *reverse* order
                code <- addRuleToCodeBlock(code, ruleList[[i]])
            return(code)
        },
        addRule = function(condition, sampler, position, name, print = FALSE) {
            '
Add a new rule for assigning sampler(s) to the samplerAssignmentRules object.  A rule consists of two parts: (1) a \'condition\' which determines when the rule is invoked, and (2) a \'sampler\' which governs the assignment of sampler(s) when the rule is invoked.  New rules can be inserted at an arbitrary position in the ordered set of rules.

Arguments:

condition: The \'condition\' argument must be a quoted R expression object, which will be evaluated and interpreted as a logical to control whether or not the rule is invoked.  The condition will be evaluated in an environment which contains the BUGS \'model\' object, the \'node\' name to which the rules (and hence the sampler assignment process) are being applied, and other sampler assignment related arguments of configureMCMC() (e.g., \'useConjugacy\' and \'multivariateNodesAsScalars\').  Thus, the condition expression may involve these names, as well as methods of BUGS model objects.  Creating an R expression object will generally use the function quote(...).  For example: addRule(condition = quote(model$isBinary(node)), ...).  Model-specific rules for particular nodes could be specified as: addRule(condition = quote(node == \'x\' || node == \'y\'), ...), or addRule(condition = quote(grepl(\'^sigma\', node)), ...).  Rules for specific distributions can be created as: addRule(condition = quote(model$getDistribution(node) == \'dpois\'), ...).  More complex, multi-line conditions can be specified as: addRule(condition = quote({ expression1;  expression2;  expression3 }), ...), the final expression of which will be evaluated and interpreted as a logical to control rule invocation.  A default, catch-all rule can specified as: addRule(condition = TRUE, ...), which will always be invoked whenever the sampler assignment process reaches this rule, and therefore subsequent rules will never be invoked.

sampler: The \'sampler\' argument controls the sampler assignment process, once a rule is invoked (i.e., the \'condition\' evaluated to TRUE).  The \'sampler\' argument must take one of three different forms: (1) a character string giving the name of an MCMC nimbleFunction sampler, (2) an unspecialized nimbleFunction object which is a valid MCMC sampler, or (3) an arbitrary quoted R expression object, which will be executed to perform the sampler assignment process, and should generally make use of the method addSampler().  Example (1): addRule(..., sampler = \'slice\'), for assigning a \'slice\' sampler when the rule is invoked. Example (2): addRule(..., sampler = my_sampler_nimbleFunction), for assigning the sampling algorithm defined in the object my_sampler_nimbleFunction.  Note the same behaviour will result from: addRule(..., sampler = \'my_sampler_nimbleFunction\'), which will be also more informative when the list of assignment rules is printed.  Example (3): addRule(..., sampler = quote( addSampler(target=node, type=\'RW\', control=list(adaptive=FALSE)) )), for assigning a non-adaptive RW sampler.  Or, perhaps, addRule(..., sampler = quote( { addSampler(target=node, type=\'RW\');  addSampler(target=node, type=\'slice\') } )), for assigning two samplers.  This third option for specifying \'sampler\' as an arbitrary R expression allows for completely general and flexible behaviour of sampler assignment.

position: Index of the position to add the new rule.  By default, new rules are added at the end of the current ordered set of rules (giving it the lowest priority in the sampler assignment process).  Specifying a position inserts the new rule at that position, and does not over-write an existing rule.

name: Optional character string name for the sampler to be added, which is used by subsequent print methods.  If \'name\' is not provided, the \'sampler\' argument is used to generate the name.  Note, if the \'sampler\' argument is provided as an R expression making use of the addSampler method, then the \'name\' argument will not be passed on to the MCMC configuration object, and instead any call(s) to addSampler can explicitly make use of its own \'name\' argument.

print: Logical argument specifying whether to print the newly-added sampler assignment rule (default FALSE).
'
            numRules <- length(ruleList)
            if(missing(position)) position <- numRules + 1  ## default adds new rules at the end
            if(position < 1) stop('cannot add new rules before position = 1')
            if(position > numRules + 1)  stop('position specified is too high')
            if(position < numRules + 1)
                ruleList[(position+1):(numRules+1)] <<- ruleList[position:numRules]
            nameProvided <- !missing(name)
            if(!nameProvided) {
                if(is.character(sampler)) {
                    name <- gsub('^sampler_', '', sampler)   ## removes 'sampler_' from beginning, if present
                } else if(is.function(sampler)) {
                    samplerArg <- substitute(sampler)
                    if(is.name(samplerArg)) {
                        name <- gsub('^sampler_', '', deparse(samplerArg))
                    } else { name <- 'custom_function' }
                } else { name <- 'if this ever is printed, then something went wrong' }  ## sampler argument is a quoted expression, such as quote({ ... })
            }
            ruleList[[position]] <<- rule(condition, sampler, name, nameProvided)
            if(print) printRules(position)
        },
        reorder = function(ind, print = FALSE) {
            '
Reorder the current ordered list of sampler assignment rules.  This method can be used to reorder the existing rules, as well as delete one or more rules.

Arguments:

ind: The indices of the current set of rules to keep.  Assuming there are 10 rules, reorder(1:5) will remove the final five rules, reorder(c(10,1:9)) will move the last (lowest priority) rule to the first position (highest priority), and reorder(8) deletes all rules except the eighth, making it the only (and hence first, highest priority) rule.

print: Logical argument specifying whether to print the resulting ordered list of sampler assignment rules (default FALSE).
'
            if(min(ind) < 1) stop('index is below one')
            if(max(ind) > length(ruleList)) stop('index is greater than number of rules')
            ruleList <<- ruleList[ind]
            if(print) printRules()
        },
        addDefaultSamplerAssignmentRules = function() {
            
	    ## posterior predictive nodes
            addRule(quote(isEndNode), 'posterior_predictive')
            
	    ## conjugate nodes
            addRule(quote(useConjugacy && isConjugate),
                    quote(addConjugateSampler(conjugacyResult = conjugacyResultsAll[[node]],
                                              dynamicallyIndexed = model$modelDef$varInfo[[model$getVarNames(nodes=node)]]$anyDynamicallyIndexed)))
            
	    ## multinomial
            addRule(quote(nodeDistribution == 'dmulti'), 'RW_multinomial')
            
	    ## dirichlet
            addRule(quote(nodeDistribution == 'ddirch'), 'RW_dirichlet')
            
	    ## wishart
            addRule(quote(nodeDistribution == 'dwish'), 'RW_wishart')
            
            ## inverse-wishart
            addRule(quote(nodeDistribution == 'dinvwish'), 'RW_wishart')
            
            ## CAR models
            addRule(quote(model$getDistribution(node) == 'dcar_normal'), 'CAR_normal')
            addRule(quote(model$getDistribution(node) == 'dcar_proper'), 'CAR_proper')

            ## multivariate & multivariateNodesAsScalars: univariate RW
            addRule(quote(isMultivariate && multivariateNodesAsScalars),
                    quote(for(scalarNode in model$expandNodeNames(node, returnScalarComponents = TRUE)) {
                        if(onlySlice) addSampler(target = scalarNode, type = 'slice')
                        else          addSampler(target = scalarNode, type = 'RW')
                    }))
            
            ## multivariate: RW_block
            addRule(quote(isMultivariate), 'RW_block')

	    ## onlyRW argument
            addRule(quote(onlyRW && !isDiscrete), 'RW')
            
	    ## onlySlice argument
            addRule(quote(onlySlice), 'slice')
            
	    ## binary-valued nodes
            addRule(quote(isBinary), 'binary')
            
	    ## categorical
            addRule(quote(nodeDistribution == 'dcat'), 'categorical')
            
	    ## discrete-valued nodes
            addRule(quote(isDiscrete), 'slice')
            
	    ## default for continuous-valued nodes: RW
            addRule(TRUE, 'RW')
        },
        printRules = function(ind) {
            '
Prints the ordered set of sampler assignment rules.

Arguments:

ind: A set of indicies, specifying which sampler assignment rules to print.  If omitted, all rules are printed.
'
            if(length(ruleList) == 0) {
                cat('Empty list of sampler assignment rules\n')
                return()
            }
            if(missing(ind)) ind <- 1:length(ruleList)   ## default is to print entire rule list
            for(i in ind) {
                cat('[[', i, ']]\n', sep='')
                print(ruleList[[i]])
                cat('\n')
            }
        },
        show = function() {
            printRules()
        }
    )
)



## set nimbleOption for configureMCMC() to default behaviour
nimbleOptions(MCMCdefaultSamplerAssignmentRules = samplerAssignmentRules())




#' Build the MCMCconf object for construction of an MCMC object
#'
#' Creates a defaut MCMC configuration for a given model.  The resulting object is suitable as an argument to \code{\link{buildMCMC}}. The assignment of sampling algorithms may be controlled using the \code{rules} argument, if provided.
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
#' @param enableWAIC A logical argument, specifying whether to enable WAIC calculations for the resulting MCMC algorithm.  Defaults to the value of \code{nimbleOptions('enableWAIC')}, which in turn defaults to FALSE.  Setting \code{nimbleOptions('enableWAIC' = TRUE)} will ensure that WAIC is enabled for all calls to \code{configureMCMC} and \code{buildMCMC}.
#'@param rules An object of class samplerAssignmentRules, which governs the assigment of MCMC sampling algorithms to stochastic model nodes.  The default set of sampler assignment rules is specified by the nimble option \'MCMCdefaultSamplerAssignmentRules\'.
#'@param warnNoSamplerAssigned A logical argument, with default value TRUE.  This specifies whether to issue a warning when no sampler is assigned to a node, meaning there is no matching sampler assignment rule.
#'@param print A logical argument, specifying whether to print the ordered list of default samplers.
#'@param autoBlock A logical argument specifying whether to use an automated blocking procedure to determine blocks of model nodes for joint sampling.  If TRUE, an MCMC configuration object will be created and returned corresponding to the results of the automated parameter blocking.  Default value is FALSE.
#'@param oldConf An optional MCMCconf object to modify rather than creating a new MCMCconf from scratch
#'@param ... Additional named control list elements for default samplers, or additional arguments to be passed to the \code{autoBlock()} function when \code{autoBlock = TRUE}
#'@author Daniel Turek
#'@export 
#'@details See \code{MCMCconf} for details on how to manipulate the \code{MCMCconf} object
#'@seealso \code{\link{samplerAssignmentRules}} \code{\link{buildMCMC}} \code{\link{runMCMC}} \code{\link{nimbleMCMC}}
configureMCMC <- function(model, nodes, control = list(), 
                          monitors, thin = 1, monitors2 = character(), thin2 = 1,
                          useConjugacy = TRUE, onlyRW = FALSE, onlySlice = FALSE, multivariateNodesAsScalars = FALSE,
                          enableWAIC = nimbleOptions('enableWAIC'),
                          print = FALSE, autoBlock = FALSE, oldConf,
                          rules = getNimbleOption('MCMCdefaultSamplerAssignmentRules'),
                          warnNoSamplerAssigned = TRUE, ...) {
    
    if(class(rules) != 'samplerAssignmentRules') stop('rules argument must be a samplerAssignmentRules object')

    if(!missing(oldConf)){
        if(!is(oldConf, 'MCMCconf'))
            stop('oldConf must be an MCMCconf object, as built by the configureMCMC function')
        return(makeNewConfFromOldConf(oldConf))	
    }
    
    if(missing(model))        stop('Either oldConf or model must be supplied')
    if(missing(nodes))        nodes <- character()
    if(missing(monitors))     monitors <- NULL

    if(autoBlock) return(autoBlock(model, ...)$conf)

    thisConf <- MCMCconf(model = model, nodes = nodes, control = control, rules = rules,
                         monitors = monitors, thin = thin, monitors2 = monitors2, thin2 = thin2,
                         useConjugacy = useConjugacy,
                         onlyRW = onlyRW, onlySlice = onlySlice,
                         multivariateNodesAsScalars = multivariateNodesAsScalars,
                         enableWAIC = enableWAIC,
                         warnNoSamplerAssigned = warnNoSamplerAssigned,
                         print = print, ...)
    return(thisConf)	
}



# This is function which builds a new MCMCconf from an old MCMCconf
# This is required to be able to a new C-based MCMC without recompiling
makeNewConfFromOldConf <- function(oldMCMCconf){
    newMCMCconf <- configureMCMC(oldMCMCconf$model, nodes = NULL)
    newMCMCconf$monitors <- oldMCMCconf$monitors
    newMCMCconf$monitors2 <- oldMCMCconf$monitors2
    newMCMCconf$thin <- oldMCMCconf$thin
    newMCMCconf$thin2 <- oldMCMCconf$thin2
    newMCMCconf$samplerConfs <- oldMCMCconf$samplerConfs
    newMCMCconf$samplerExecutionOrder <- oldMCMCconf$samplerExecutionOrder
    newMCMCconf$controlDefaults <- oldMCMCconf$controlDefaults
    newMCMCconf$namedSamplerLabelMaker <- oldMCMCconf$namedSamplerLabelMaker
    newMCMCconf$mvSamples1Conf <- oldMCMCconf$mvSamples1Conf
    newMCMCconf$mvSamples2Conf <- oldMCMCconf$mvSamples2Conf
    return(newMCMCconf)	
}


newSpacesFunction <- function(m) {
    log10max <- floor(log10(m))
    function(i) paste0(rep(' ', log10max-floor(log10(i))), collapse = '')
}



