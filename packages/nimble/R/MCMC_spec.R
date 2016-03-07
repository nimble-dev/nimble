

controlDefaultList <- list(
    adaptive = TRUE,
    adaptScaleOnly = FALSE,
    adaptInterval = 200,
    scale = 1,
    propCov = 'identity',
    sliceWidth = 1,
    sliceMaxSteps = 100,
    m = 1000
)



samplerSpec <- setRefClass(
    Class = 'samplerSpec',
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
        buildSampler = function(model, mvSaved) {
            samplerFunction(model=model, mvSaved=mvSaved, target=target, control=control)
        },
        toStr = function() {
            tempList <- list()
            tempList[[paste0(name, ' sampler')]] <- paste0(target, collapse = ', ')
            mcmc_listContentsToStr(c(tempList, control))
        }
    )
)


## NOTE: methods are documented as a "docstring" with each method - see 'removeSamplers' below. roxygen will automatically grab info from these docstrings and inject into the Rd in the Methods Section
## NOTE: including the name of the class in @aliases is important because by default we only get help("MCMCspec-class") and not help(MCMCspec)
## NOTE: the empty lines are important in the final formatting, so please don't remove any of them in your own help info

#' Class \code{MCMCspec}
#' @aliases MCMCspec addSampler removeSamplers setSamplers getSamplers addMonitors addMonitors2 resetMonitors getMonitors setThin setThin2
#' @export
#' @description
#' Objects of this class specify an MCMC algorithm, specific to a particular model.  Objects are normally created by calling \link{configureMCMC}.
#' Given an MCMCspec object, the actual MCMC function can be built by calling \link{buildMCMC}\code{(spec)}.
#' See documentation below for method initialize() for details of creating an MCMCspec object.
#' @author Daniel Turek
#' @seealso \link{configureMCMC}
#' @examples
#' code <- nimbleCode({
#'  mu ~ dnorm(0, 1)
#'  x ~ dnorm(mu, 1)
#' })
#' Rmodel <- nimbleModel(code)
#' spec <- configureMCMC(Rmodel)
#' spec$setSamplers(1)
#' spec$addSampler(target = 'x', type = 'slice', control = list(adaptInterval = 100))
#' spec$addMonitors('mu')
#' spec$addMonitors2('x')
#' spec$setThin(5)
#' spec$setThin2(10)
#' spec$getMonitors()
#' spec$getSamplers()
MCMCspec <- setRefClass(
    
    Class = 'MCMCspec',                           
    
    fields = list(
        model               = 'ANY',
        monitors            = 'ANY',
        monitors2           = 'ANY',
        thin                = 'ANY',
        thin2               = 'ANY',
        samplerSpecs        = 'ANY',
        controlDefaults     = 'ANY',
        controlNamesLibrary = 'ANY',
        mvSamples1Spec      = 'ANY',
        mvSamples2Spec      = 'ANY'
    ),
    
    methods = list(
        
        initialize = function(model, nodes, control = list(),
            monitors,                thin  = 1,
            monitors2 = character(), thin2 = 1,
            useConjugacy = TRUE, onlyRW = FALSE, onlySlice = FALSE, multivariateNodesAsScalars = FALSE,
            print = FALSE) {	
            '
Creates a MCMC specification for a given model.  The resulting object is suitable as an argument to buildMCMC.

Arguments:

model: A NIMBLE model object, created from nimbleModel(...)

nodes: An optional character vector, specifying the nodes for which samplers should be created.
Nodes may be specified in their indexed form, \'y[1, 3]\', or nodes specified without indexing will be expanded fully, e.g., \'x\' will be expanded to \'x[1]\', \'x[2]\', etc.
If missing, the default value is all non-data stochastic nodes.
If NULL, then no samplers are added.

control: An optional list of control arguments to sampler functions.  If a control list is provided, the elements will be provided to all sampler functions which utilize the named elements given.
For example, the standard Metropolis-Hastings random walk sampler (sampler_RW) utilizes control list elements \'adaptive\', \'adaptInterval\', \'scale\', 
and also \'targetNode\' however this should not generally be provided as a control list element to configureMCMC().
The default values for control list arguments for samplers (if not otherwise provided as an argument to configureMCMC) are in the NIMBLE system option \'MCMCcontrolDefaultList\'.

monitors: A character vector of node names or variable names, to record during MCMC sampling.
This set of monitors will be recorded with thinning interval \'thin\', and the samples will be stored into the \'mvSamples\' object.
The default value is all top-level stochastic nodes of the model -- those having no stochastic parent nodes.

monitors2: A character vector of node names or variable names, to record during MCMC sampling.
This set of monitors will be recorded with thinning interval \'thin2\', and the samples will be stored into the \'mvSamples2\' object.
The default value is an empty character vector, i.e. no values will be recorded.

thin: The thinning interval for \'monitors\'.  Default value is one.

thin2: The thinning interval for \'monitors2\'.  Default value is one.

useConjugacy: A logical argument, with default value TRUE.  If specified as FALSE, then no conjugate samplers will be used, even when a node is determined to be in a conjugate relationship.

onlyRW: A logical argument, with default value FALSE.  If specified as TRUE, then Metropolis-Hastings random walk samplers (sampler_RW) will be assigned for all non-terminal continuous-valued nodes nodes.
Discrete-valued nodes are assigned a slice sampler (sampler_slice), and terminal (predictive) nodes are assigned an end sampler (sampler_end).

onlySlice: A logical argument, with default value FALSE.  If specified as TRUE, then a slice sampler is assigned for all non-terminal nodes.
Terminal (predictive) nodes are still assigned an end sampler (sampler_end).

multivariateNodesAsScalars: A logical argument, with default value FALSE.  If specified as TRUE, then non-terminal multivariate stochastic nodes will have scalar samplers assigned to each of the scalar components of the multivariate node.  The default value of FALSE results in a single block sampler assigned to the entire multivariate node.  Note, multivariate nodes appearing in conjugate relationships will be assigned the corresponding conjugate sampler (provided useConjugacy == TRUE), regardless of the value of this argument.

print: A logical argument, specifying whether to print the ordered list of default samplers.
'
            
            samplerSpecs <<- list(); controlDefaults <<- list(); controlNamesLibrary <<- list(); monitors <<- character(); monitors2 <<- character();
            model <<- model
            addMonitors( monitors,  print = FALSE)
            addMonitors2(monitors2, print = FALSE)
            thin  <<- thin
            thin2 <<- thin2
            samplerSpecs    <<- list()
            ## moved controlDefaultList to be a NIMBLE system option (as a single list: MCMCcontrolDefaultList)
            controlDefaults <<- getNimbleOption('MCMCcontrolDefaultList')
            for(i in seq_along(control))     controlDefaults[[names(control)[i]]] <<- control[[i]]
            controlNamesLibrary <<- list()
            if(identical(nodes, character())) { nodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
                                            } else             { if(is.null(nodes) || length(nodes)==0)     nodes <- character(0)
                                                                 nl_checkVarNamesInModel(model, removeIndexing(nodes))
                                                                 nodes <- model$expandNodeNames(nodes)            }
            
            nodes <- model$topologicallySortNodes(nodes)   ## topological sort
            isNodeEnd <- model$isNodeEnd(nodes)


            if(useConjugacy) conjugacyResultsAll <- model$checkConjugacy2(nodes)

            for(i in seq_along(nodes)) {
            	node <- nodes[i]
                discrete <- model$isDiscrete(node)
                nodeScalarComponents <- model$expandNodeNames(node, returnScalarComponents = TRUE)
                nodeLength <- length(nodeScalarComponents)
                
                ## if node has 0 stochastic dependents, assign 'end' sampler (e.g. for predictive nodes)
             	if(isNodeEnd[i]) { addSampler(target = node, type = 'end', print = print);     next }
                
                ## for multivariate nodes, either add a conjugate sampler, or RW_block sampler
                if(nodeLength > 1) {
                    if(useConjugacy) {
                        conjugacyResult <- conjugacyResultsAll[[node]]
                        if(!is.null(conjugacyResult)) {
                            addSampler(target = conjugacyResult$target, type = conjugacyResult$type, control = conjugacyResult$control, print = print);     next }
                    }
                    if(multivariateNodesAsScalars) {
                        for(scalarNode in nodeScalarComponents) {
                            addSampler(target = scalarNode, type = 'RW', print = print) };     next }
                    addSampler(target = node, type = 'RW_block', print = print);     next }

                ## node is scalar, non-end node
                if(onlyRW && !discrete)   { addSampler(target = node, type = 'RW',    print = print);     next }
                if(onlySlice)             { addSampler(target = node, type = 'slice', print = print);     next }
                
                ## if node passes checkConjugacy(), assign 'conjugate_dxxx' sampler
                if(useConjugacy) {
                    conjugacyResult <- conjugacyResultsAll[[node]]
                    if(!is.null(conjugacyResult)) {
                        addSampler(target = conjugacyResult$target, type = conjugacyResult$type, control = conjugacyResult$control, print = print);     next }
                }
                
                ## if node distribution is discrete, assign 'slice' sampler
                if(discrete) { addSampler(target = node, type = 'slice', print = print);     next }
                
                ## default: 'RW' sampler
                addSampler(target = node, type = 'RW', print = print);     next
            }
        },
        
        addSampler = function(target, type = 'RW', control = list(), print = TRUE, name) {
            '
Adds a sampler to the list of samplers contained in the MCMCspec object.

Arguments:

target: The target node or nodes to be sampled.  This may be specified as a character vector of model node and/or variable names.  This argument is required.

type: The type of sampler to add, specified as either a character string or a nimbleFunction object.  If the character argument type=\'newSamplerType\', then either newSamplerType or sampler_newSamplertype must correspond to a nimbleFunction (i.e. a function returned by nimbleFunction, not a specialized nimbleFunction).  Alternatively, the type argument may be provided as a nimbleFunction itself rather than its name.  In that case, the \'name\' argument may also be supplied to provide a meaningful name for this sampler.  The default value is \'RW\' which specifies scalar adaptive Metropolis-Hastings sampling with a normal proposal distribution. This default will result in an error if \'target\' specifies more than one target node.

control: A list of control arguments specific to the sampler function.
These will override the defaults provided in the NIMBLE system option \'MCMCcontrolDefaultList\', and any specified in the control list argument to configureMCMC().
An error results if the sampler function requires any control elements which are 
not present in this argument, the control list argument to configureMCMC(), or in the NIMBLE system option \'MCMCcontrolDefaultList\'.

print: Logical argument, specifying whether to print the details of the newly added sampler, as well as its position in the list of MCMC samplers.

name: A character string name for the sampler, which is only used when the \'type\' argument is provided as a nimbleFunction generator object.  If \'name\' is not provided, then the text of the \'type\' argument (i.e. the result of \'deparse(substitute(type))\') is used as the default sampler name.

Details: A single instance of the newly specified sampler is added to the end of the list of samplers for this MCMCspec object.

Invisibly returns a list of the current sampler specifications, which are samplerSpec reference class objects.
'

            if(is.character(type)) {
                if(type == 'end') type <- 'sampler_end'  ## because 'end' is an R function
                thisSamplerName <- gsub('^sampler_', '', type)   ## removes 'sampler_' from beginning of name, if present
                if(exists(type)) {   ## try to find sampler function 'type'
                    samplerFunction <- eval(as.name(type))
                } else {
                    sampler_type <- paste0('sampler_', type)   ## next, try to find sampler function 'sampler_type'
                    if(exists(sampler_type)) {
                        samplerFunction <- eval(as.name(sampler_type))
                    } else stop(paste0('cannot find sampler type \'', type, '\''))
                }
            } else if(is.function(type)) {
                thisSamplerName <- if(!missing(name)) name else gsub('^sampler_', '', deparse(substitute(type)))
                samplerFunction <- type
            } else stop('sampler type must be character name or function')
            if(!is.character(thisSamplerName)) stop('Sampler name should be a character string')
            if(!is.function(samplerFunction)) stop('Sampler type does not specify a function')

            if(is.null(controlNamesLibrary[[thisSamplerName]]))   controlNamesLibrary[[thisSamplerName]] <<- mcmc_findControlListNamesInCode(samplerFunction)   ## populate control names library
            controlListNames <- controlNamesLibrary[[thisSamplerName]]
            thisControlList <- controlDefaults           ## start with all the defaults
            thisControlList[names(control)] <- control   ## add in any controls provided as an argument
            missingControlNames <- setdiff(controlListNames, names(thisControlList))
            missingControlNames <- missingControlNames[!grepl('^dependents_', missingControlNames)]   ## dependents for conjugate samplers are exempted from this check
            if(length(missingControlNames) != 0)  stop(paste0('Required control names are missing for ', thisSamplerName, ' sampler: ', paste0(missingControlNames, collapse=', ')))
            if(!all(names(control) %in% controlListNames))   warning(paste0('Superfluous control names were provided for ', thisSamplerName, ' sampler: ', paste0(setdiff(names(control), controlListNames), collapse=', ')))
            thisControlList <- thisControlList[controlListNames]
            
            newSamplerInd <- length(samplerSpecs) + 1
            samplerSpecs[[newSamplerInd]] <<- samplerSpec(name=thisSamplerName, samplerFunction=samplerFunction, target=target, control=thisControlList, model=model)
            
            if(print) getSamplers(newSamplerInd)
            return(invisible(samplerSpecs))
        },
        
        removeSamplers = function(ind, print = TRUE) {
            '
Removes one or more samplers from an MCMCspec object.

Arguments:

ind: A numeric vector or character vector specifying the samplers to remove.  A numeric vector may specify the indices of the samplers to be removed.  Alternatively, a character vector may be used to specify a set of model nodes and/or variables, and all samplers whose \'target\' is among these nodes will be removed.  If omitted, then all samplers are removed.

print: A logical argument, default value TRUE, specifying whether to print the current list of samplers once the removal has been done.
'      
            if(missing(ind))        ind <- seq_along(samplerSpecs)
            if(is.character(ind))   ind <- findSamplersOnNodes(ind)
            if(length(ind) > 0 && max(ind) > length(samplerSpecs)) stop('MCMC specification doesn\'t have that many samplers')
            samplerSpecs[ind] <<- NULL
            if(print) getSamplers()
            return(invisible(NULL))
        },
        
        setSamplers = function(ind, print = TRUE) {
            '
Sets the ordering of the list of MCMC samplers.

Arguments:

ind: A numeric vector or character vector.  A numeric vector may be used to specify the indicies for the new list of MCMC samplers, in terms of the current ordered list of samplers.
For example, if the MCMCspec object currently has 3 samplers, then the ordering may be reversed by calling mcmcspec$setSamplers(3:1), or all samplers may be removed by calling mcmcspec$setSamplers(numeric(0)).  Alternatively, a character vector may be used to specify a set of model nodes and/or variables, and the sampler list will modified to only those samplers acting on these target nodes.

print: A logical argument, default value TRUE, specifying whether to print the new list of samplers.
'   
            if(missing(ind))        ind <- numeric(0)
            if(is.character(ind))   ind <- findSamplersOnNodes(ind)
            if(length(ind) > 0 && max(ind) > length(samplerSpecs)) stop('MCMC specification doesn\'t have that many samplers')
            samplerSpecs <<- samplerSpecs[ind]
            if(print) getSamplers()
            return(invisible(NULL))
        },
        
        getSamplers = function(ind) {
            '
Prints details of the MCMC samplers.

Arguments:

ind: A numeric vector or character vector.  A numeric vector may be used to specify the indices of the samplers to print, or a character vector may be used to indicate a set of target nodes and/or variables, for which all samplers acting on these nodes will be printed. For example, getSamplers(\'x\') will print all samplers whose target is model node \'x\', or whose targets are contained (entirely or in part) in the model variable \'x\'.  If omitted, then all samplers are printed.

Invisibly returns a list of the current sampler specifications for the specified samplers, which are samplerSpec reference class objects.
'
            if(missing(ind))        ind <- seq_along(samplerSpecs)
            if(is.character(ind))   ind <- findSamplersOnNodes(ind)
            if(length(ind) > 0 && max(ind) > length(samplerSpecs)) stop('MCMC specification doesn\'t have that many samplers')
            makeSpaces <- if(length(ind) > 0) newSpacesFunction(max(ind)) else NULL
            for(i in ind)
                cat(paste0('[', i, '] ', makeSpaces(i), samplerSpecs[[i]]$toStr(), '\n'))
            if(length(ind) == 1) return(invisible(samplerSpecs[[ind]]))
            return(invisible(samplerSpecs[ind]))
        },

        findSamplersOnNodes = function(nodes) {
            if(length(samplerSpecs) == 0) return(integer())
            nodes <- model$expandNodeNames(nodes, returnScalarComponents = TRUE)
            which(unlist(lapply(samplerSpecs, function(ss) any(nodes %in% ss$targetAsScalar))))
        },
        
        addMonitors = function(vars, ind = 1, print = TRUE) {
            '
Adds variables to the list of monitors.

Arguments:

vars: A character vector of indexed nodes, or variables, which are to be monitored.  These are added onto the current monitors list.

print: A logical argument, specifying whether to print all current monitors.

Details: See the initialize() function
            '
            
            if(isMvSamplesReady(ind)){
            	cat('Changing monitors, even though an MCMC has been built already. When compiling the MCMC, use resetFunctions = TRUE option\n')
            	if(ind == 1)
                    mvSamples1Spec <<- NULL
            	if(ind == 2)
                    mvSamples2Spec <<- NULL
            }
            
            
            
            if(is.null(vars))  vars <- model$getNodeNames(topOnly = TRUE, stochOnly = TRUE)
            vars <- unique(removeIndexing(vars))
            nl_checkVarNamesInModel(model, vars)
            if(ind == 1)     monitors  <<- unique(c(monitors,  vars))
            if(ind == 2)     monitors2 <<- unique(c(monitors2, vars))
            if(print) getMonitors()
            return(invisible(NULL))
        },

        addMonitors2 = function(vars, print = TRUE) {
            '
Adds variables to the list of monitors2.

Arguments:

vars: A character vector of indexed nodes, or variables, which are to be monitored.  These are added onto the current monitors2 list.

print: A logical argument, specifying whether to print all current monitors.

Details: See the initialize() function
            '
            addMonitors(vars, ind = 2, print = print)
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
            	mvSamples1Spec <<- NULL
            	mvSamples2Spec <<- NULL
            }

            
            return(invisible(NULL))
        },
        
        getMonitors = function() {
            '
Prints all current monitors and monitors2

Details: See the initialize() function
            '
            if(length(monitors)  > 0)   cat(paste0('thin = ', thin,  ': ', paste0(monitors,  collapse = ', '), '\n'))
            if(length(monitors2) > 0)   cat(paste0('thin2 = ', thin2, ': ', paste0(monitors2, collapse = ', '), '\n'))
        },
        
        setThin  = function(thin, print = TRUE) {
            '
Sets the value of thin.

Arguments:

thin: The new value for the thinning interval \'thin\'.

print: A logical argument, specifying whether to print all current monitors.

Details: See the initialize() function
            '
            thin  <<- thin
            if(print) getMonitors()
            return(invisible(NULL))
        },
        setThin2 = function(thin2, print = TRUE) {
            '
Sets the value of thin2.

Arguments:

thin2: The new value for the thinning interval \'thin2\'.

print: A logical argument, specifying whether to print all current monitors.

Details: See the initialize() function
            '
            thin2 <<- thin2
            if(print) getMonitors()
            return(invisible(NULL))
        },
        
        getMvSamplesSpec  = function(ind = 1){
            
            if(isMvSamplesReady(ind) == TRUE) {
                if(ind == 1) return(mvSamples1Spec)
                return(mvSamples2Spec)
            }
            else{
                makeMvSamplesSpec(ind)
                if(ind == 1)
                    output <- mvSamples1Spec
                if(ind == 2)
                    output <- mvSamples2Spec
                return(output)
            }
        },
        
        isMvSamplesReady = function(ind){
            if(ind == 1) return(is(mvSamples1Spec, 'function'))		#Probably really want to give mvSpecs there own class...
            if(ind == 2) return(is(mvSamples2Spec, 'function'))
            stop('invalid indicator for isMvSsamplesReady')
        },
        
        makeMvSamplesSpec = function(ind){
            modelSymbolObjects = model$getSymbolTable()$getSymbolObjects()
            if(ind == 1) monitorNames = monitors
            if(ind == 2) monitorNames = monitors2
            if(!all(monitorNames %in% names(modelSymbolObjects))) stop('some monitor names are not in the model symbol table; this should never occur')
            thisModelValuesSpec = modelValuesSpec(symbolTable(symbols = modelSymbolObjects[monitorNames]))
            if(ind == 1) mvSamples1Spec <<- thisModelValuesSpec
            if(ind == 2) mvSamples2Spec <<- thisModelValuesSpec     	
        },

        show = function() {
            cat('MCMC specification object\n')
        }
    )
)



## This appeared to be roxygen content that didn't belong here.  It is in BUGS_readBUGS.R.  -Perry 8/1/15 
## # Turn BUGS model code into an object for use in nimbleModel or readBUGSmodel
## #
## # Simply keeps model code as an R call object, the form needed by \code{nimbleModel} and optionally usable by \code{readBUGSmodel}
## # 
## # @param code expression providing the code for the model 
## # @author Daniel Turek
## # @export
## # @details It is equivalent to use the R function \code{quote}.  \code{nimbleCode} is simply provided as a more readable alternative for NIMBLE users not familiar with \code{quote}.
## # @examples
## # code <- nimbleCode({
## #     x ~ dnorm(mu, sd = 1)
## #     mu ~ dnorm(0, sd = prior_sd)
## # })




#' Build the MCMCspec object for construction of an MCMC object
#'
#' Creates a defaut MCMC specification for a given model.  The resulting object is suitable as an argument to \link{buildMCMC}. 
#'
#'@param model A NIMBLE model object, created from \link{nimbleModel}
#'@param nodes An optional character vector, specifying the nodes and/or variables for which samplers should be created.
#'Nodes may be specified in their indexed form, \code{y[1, 3]}.  Alternatively, nodes specified without indexing will be expanded fully, e.g., \code{x} will be expanded to \code{x[1]}, \code{x[2]}, etc.
#'If missing, the default value is all non-data stochastic nodes.
#'If NULL, then no samplers are added.
#'@param control An optional list of control arguments to sampler functions.  If a control list is provided, the elements will be provided to all sampler functions which utilize the named elements given.
#'For example, the standard Metropolis-Hastings random walk sampler (\link{sampler_RW}) utilizes control list elements \code{adaptive}, \code{adaptInterval}, and \code{scale}.
#' (Internally it also uses \code{targetNode}, but this should not generally be provided as a control list element).
#'The default values for control list arguments for samplers (if not otherwise provided as an argument to configureMCMC() ) are in the NIMBLE system option \code{MCMCcontrolDefaultList}.
#'@param monitors A character vector of node names or variable names, to record during MCMC sampling.
#'This set of monitors will be recorded with thinning interval \code{thin}, and the samples will be stored into the \code{mvSamples} object.
#'The default value is all top-level stochastic nodes of the model -- those having no stochastic parent nodes.
#'@param monitors2 A character vector of node names or variable names, to record during MCMC sampling.
#'This set of monitors will be recorded with thinning interval \code{thin2}, and the samples will be stored into the \code{mvSamples2} object.
#'The default value is an empty character vector, i.e. no values will be recorded.
#'@param thin The thinning interval for \code{monitors}.  Default value is one.
#'@param thin2 The thinning interval for \code{monitors2}.  Default value is one.
#'@param useConjugacy A logical argument, with default value TRUE.  If specified as FALSE, then no conjugate samplers will be used, even when a node is determined to be in a conjugate relationship.
#'@param onlyRW A logical argument, with default value FALSE.  If specified as TRUE, then Metropolis-Hastings random walk samplers (\link{sampler_RW}) will be assigned for all non-terminal continuous-valued nodes nodes.
#'Discrete-valued nodes are assigned a slice sampler (\link{sampler_slice}), and terminal (predictive) nodes are assigned an end sampler (\link{sampler_end}).
#'@param onlySlice A logical argument, with default value FALSE.  If specified as TRUE, then a slice sampler is assigned for all non-terminal nodes.
#'Terminal (predIctive) nodes are still assigned an end sampler (sampler_end).
#'@param multivariateNodesAsScalars A logical argument, with default value FALSE.  If specified as TRUE, then non-terminal multivariate stochastic nodes will have scalar samplers assigned to each of the scalar components of the multivariate node.  The default value of FALSE results in a single block sampler assigned to the entire multivariate node.  Note, multivariate nodes appearing in conjugate relationships will be assigned the corresponding conjugate sampler (provided \code{useConjugacy == TRUE}), regardless of the value of this argument.
#'@param print A logical argument, specifying whether to print the ordered list of default samplers.
#'@param autoBlock A logical argument specifying whether to use an automated blocking procedure to determine blocks of model nodes for joint sampling.  If TRUE, an MCMC specification object will be created and returned corresponding to the results of the automated parameter blocking.  Default value is FALSE.
#'@param oldSpec An optional MCMCspec object to modify rather than creating a new MCMCspec from scratch
#'@param ... Additional arguments to be passed to the \code{autoBlock()} function when \code{autoBlock = TRUE}
#' @author Daniel Turek
#'@details See \code{MCMCspec} for details on how to manipulate the \code{MCMCspec} object
configureMCMC <- function(model, nodes, control = list(), 
                          monitors, thin = 1, monitors2 = character(), thin2 = 1,
                          useConjugacy = TRUE, onlyRW = FALSE, onlySlice = FALSE, multivariateNodesAsScalars = FALSE,
                          print = FALSE, autoBlock = FALSE, oldSpec, ...) {
    
    if(!missing(oldSpec)){
        if(!is(oldSpec, 'MCMCspec'))
            stop('oldSpec must be an MCMCspec object, as built by the configureMCMC function')
        return(makeNewSpecFromOldSpec(oldSpec))	
    }
    
    if(missing(model))        stop('Either oldSpec or model must be supplied')
    if(missing(nodes))        nodes <- character()
    if(missing(monitors))     monitors <- NULL

    if(autoBlock) return(autoBlock(model, ...)$spec)
    
    thisSpec <- MCMCspec(model = model, nodes = nodes, control = control, 
                         monitors = monitors, thin = thin, monitors2 = monitors2, thin2 = thin2,
                         useConjugacy = useConjugacy, onlyRW = onlyRW, onlySlice = onlySlice,
                         multivariateNodesAsScalars = multivariateNodesAsScalars, print = print)
    return(thisSpec)	
}



# This is function which builds a new MCMCspec from an old MCMCspec
# This is required to be able to a new C-based MCMC without recompiling
makeNewSpecFromOldSpec <- function(oldMCMCspec){
    newMCMCspec <- configureMCMC(oldMCMCspec$model, nodes = NULL)
    newMCMCspec$monitors <- oldMCMCspec$monitors
    newMCMCspec$monitors2 <- oldMCMCspec$monitors2
    newMCMCspec$thin <- oldMCMCspec$thin
    newMCMCspec$thin2 <- oldMCMCspec$thin2
    newMCMCspec$samplerSpecs <- oldMCMCspec$samplerSpecs
    newMCMCspec$controlDefaults <- oldMCMCspec$controlDefaults
    newMCMCspec$controlNamesLibrary <- oldMCMCspec$controlNamesLibrary
    newMCMCspec$mvSamples1Spec <- oldMCMCspec$mvSamples1Spec
    newMCMCspec$mvSamples2Spec <- oldMCMCspec$mvSamples2Spec
    return(newMCMCspec)	
}


newSpacesFunction <- function(m) {
    log10max <- floor(log10(m))
    function(i) paste0(rep(' ', log10max-floor(log10(i))), collapse = '')
}




