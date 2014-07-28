



conjugacyRelationshipsInputList <- list(
    
    list(prior = 'beta',
         dependents = list(
             bern   = list(param = 'prob', link = 'identity', contribution_shape1 = 'value', contribution_shape2 = '1 - value'   ),
             bin    = list(param = 'prob', link = 'identity', contribution_shape1 = 'value', contribution_shape2 = 'size - value'),
             negbin = list(param = 'prob', link = 'identity', contribution_shape1 = 'size',  contribution_shape2 = 'value'       )
         ),
         posterior = 'beta(shape1 = prior_shape1 + contribution_shape1,
                           shape2 = prior_shape2 + contribution_shape2)'),
    
    list(prior = 'gamma',
         dependents = list(
             pois  = list(param = 'lambda', link = 'multiplicative', contribution_shape = 'value', contribution_rate = 'coeff'                           ),
             norm  = list(param = 'tau',    link = 'multiplicative', contribution_shape = '1/2',   contribution_rate = 'coeff/2 * (value-mean)^2'        ),
             lnorm = list(param = 'tau',    link = 'multiplicative', contribution_shape = '1/2',   contribution_rate = 'coeff/2 * (log(value)-meanlog)^2'),
             gamma = list(param = 'rate',   link = 'multiplicative', contribution_shape = 'shape', contribution_rate = 'coeff   * value'                 ),
             exp   = list(param = 'rate',   link = 'multiplicative', contribution_shape = '1',     contribution_rate = 'coeff   * value'                 )
             ## dexp  = list(param = 'rate',   link = 'multiplicative', contribution_shape = '1',     contribution_rate = 'coeff   * abs(value-location)'   )
             ## par = list(...)    ## need to figure this out
         ),
         posterior = 'gamma(shape = prior_shape + contribution_shape,
                            scale = 1 / (prior_rate + contribution_rate))'),
    
    list(prior = 'norm',
         dependents = list(
             norm  = list(param = 'mean',    link = 'linear', contribution_mean = 'coeff * (value-offset) * tau',      contribution_tau = 'coeff^2 * tau'),
             lnorm = list(param = 'meanlog', link = 'linear', contribution_mean = 'coeff * (log(value)-offset) * tau', contribution_tau = 'coeff^2 * tau')
         ),
         posterior = 'norm(mean = (prior_mean*prior_tau + contribution_mean) / (prior_tau + contribution_tau),
                           sd   = (prior_tau + contribution_tau)^(-0.5))')
    
    ##### waiting for dpar() distribution
    # list(prior = 'par',
    #      dependents = list(
    #          unif = list(param = 'max', link = 'multiplicative', contribution_alpha = '1', contribution_not_used = 'coeff'),
    #          par  = list(param = 'c',   link = 'multiplicative', contribution_alpha = '-alpha')
    #      ),
    #      posterior = 'par(alpha = prior_alpha + contribution_alpha,
    #                       c     = max(prior_c, max(dependents_dunif_values/dependents_dunif_coeff)))'),
    #####
)




##############################################################################################
##############################################################################################
## reference class definitions:
## conjugacyRelationshipsClass
## conjugacyClass
## dependentClass
## posteriorClass
##############################################################################################
##############################################################################################

conjugacyRelationshipsClass <- setRefClass(
    Class = 'conjugacyRelationshipsClass',
    fields = list(
        conjugacys = 'ANY'  ## a (named) list of conjugacyClass objects, each describes the conjugacies for a particular prior distribution (name is prior distribution name)
    ),
    methods = list(
        initialize = function(crl) {
        	conjugacys <<- list()
            for(i in seq_along(crl)) {
                conjugacys[[i]] <<- conjugacyClass(crl[[i]])
            }
            names(conjugacys) <<- unlist(lapply(conjugacys, function(cr) cr$prior))
        },
        checkConjugacy = function(model, targetNode) {
            if(!(targetNode %in% model$getNodeNames()))       stop('checking conjugacy of a node not in model')
            if(model$getNodeInfo()[[targetNode]]$type != 'stoch')  stop('checking conjugacy of non-stochastic node')
            depNodes <- model$getDependencies(targetNode, stochOnly = TRUE, self = FALSE)
            if(length(depNodes) == 0)  return(NULL)   # no dependent stochastic nodes: not conjugate, return NULL
            
            for(conjugacyObj in conjugacys) {  # conjugacyObj is a conjugacyClass object
                conjugacyResult <- conjugacyObj$checkConjugacy(model, targetNode, depNodes)    ## workhorse for checking conjugacy
                if(is.null(conjugacyResult))     next
                return(conjugacyResult)
            }
            return(NULL)  # didn't find a matching conjugacy class: not conjugate, return NULL
        },
        generateConjugateSamplerDefinitions = function() {
            conjugateSamplerDefinitions <- list()
            for(conjugacyObj in conjugacys) {  # conjugacyObj is a conjugacyClass object
                samplerName <- cc_makeConjugateSamplerName(conjugacyObj$samplerType)
                conjugateSamplerDefinitions[[samplerName]] <- conjugacyObj$generateConjugateSamplerDef()    ## workhorse for creating conjugate sampler nimble functions
            }
            return(conjugateSamplerDefinitions)
        }
    )
)

conjugacyClass <- setRefClass(
    Class = 'conjugacyClass',
    fields = list(
        samplerType = 			'ANY', 		## name of the sampler for this conjugacy class, e.g. 'conjugate_norm'
        prior =					'ANY', 		## name of the prior distribution, e.g. 'dnorm'
        dependents = 			'ANY', 		## (named) list of dependentClass objects, each contains conjugacy information specific to a particular sampling distribution (name is sampling distribution name)
        dependentDistNames = 	'ANY', 		## character vector of the names of all allowable dependent sampling distributions.  same as: names(dependents)
        posteriorObject = 		'ANY',   	## an object of posteriorClass
        needsLinearityCheck = 	'ANY', 		## logical specifying whether we need to do the linearity check; if any dependents require either 'coeff' or 'offset'
        model = 				'ANY' 	    ## ONLY EXISTS TO PREVENT A WARNING for '<<-', in the code for generating the conjugate sampler function
    ),
    methods = list(
        initialize = function(cr) {
        	dependents <<- list()
            samplerType <<- cc_makeSamplerTypeName(cr$prior)
            prior <<- cc_makeDDistributionName(cr$prior)
            initialize_addDependents(cr$dependents)
            needsLinearityCheck <<- any(c(unlist(lapply(dependents, function(dep) dep$needsCoeff)), unlist(lapply(dependents, function(dep) dep$needsOffset))))
            posteriorObject <<- posteriorClass(cr$posterior)
        },
        initialize_addDependents = function(depList) {
            for(i in seq_along(depList)) {
                dependents[[i]] <<- dependentClass(depList[[i]], names(depList)[i])
            }
            names(dependents) <<- cc_makeDDistributionName(names(depList))
            dependentDistNames <<- names(dependents)
        },
        
        ## workhorse for checking conjugacy
        checkConjugacy = function(model, targetNode, depNodes) {
            if(cc_getNodeDistributionText(model, targetNode) != prior)     return(NULL)    # check prior distribution of targetNode
            control <- initControl(model, targetNode)
            
            for(depNode in depNodes) {
                depNodeDist <- cc_getNodeDistributionText(model, depNode)
                if(!(depNodeDist %in% dependentDistNames))     return(NULL)    # check sampling distribution of depNode
                dependentObj <- dependents[[depNodeDist]]
                linearityCheckExpr <- cc_getNodeParamExpr(model, depNode, dependentObj$param)   # extracts the expression for 'param' from 'depNode'
                linearityCheckExpr <- cc_expandDetermNodesInExpr(linearityCheckExpr, model)
                ## next line is a NEW ADDITION, prevents a minor bug in conjugacy checking:
                ## when targetNode doesn't appear in 'param' expr (hence passes the linearlity check),
                ## and targetNode appears in *exactly one* other parameter expr (hence passing cc_otherParamsCheck()),
                ## which also explains why depNode is identified as a dependent node in the first place.
                ## we simply ensure that targetNode actually does appear in the conjugate parameter expression,
                ## thus the conjugacy check will fail if targetNode appears in any other parameter expressions (failing in cc_otherParamsCheck())
                if(!cc_nodeInExpr(targetNode, linearityCheckExpr))       return(NULL)
                linearityCheck <- cc_checkLinearity(linearityCheckExpr, targetNode)   # determines whether paramExpr is linear in targetNode
                if(!cc_linkCheck(linearityCheck, dependentObj$link))     return(NULL)
                if(!cc_otherParamsCheck(model, depNode, targetNode))     return(NULL)   # ensure targetNode appears in only *one* depNode parameter expression
                control[[paste0('dependents_', depNodeDist)]] <- c(control[[paste0('dependents_', depNodeDist)]], depNode)
            }
            return(list(samplerType=samplerType, control=control))   # all dependent nodes passed the conjugacy check
        },
        
        initControl = function(model, targetNode) {
            control <- list()
            control$targetNode <- targetNode
            for(depDist in dependentDistNames)     control[[paste0('dependents_', depDist)]] <- character()
            return(control)
        },
        
        ## workhorse for creating conjugate sampler nimble functions
        generateConjugateSamplerDef = function() {
            substitute(
                nimbleFunction(contains = sampler_BASE,
                               setup    = SETUPFUNCTION,
                               run      = RUNFUNTION,
                               methods  = list(getPosteriorLogDensity = GETPOSTERIORLOGDENSITYFUNCTION,
                                               reset                  = function() {}),
                               where    = getLoadingNamespace()
                ),
                list(SETUPFUNCTION                  = genSetupFunction(),
                     RUNFUNTION                     = genRunFunction(),
                     GETPOSTERIORLOGDENSITYFUNCTION = genGetPosteriorLogDensityFunction()
                )
            )
        },
        
        genSetupFunction = function() {
            functionBody <- codeBlockClass()
            
            ## preliminaries
            functionBody$addCode({
                targetNode <- control$targetNode
                calcNodes       <- model$getDependencies(targetNode)
                calcNodesDeterm <- model$getDependencies(targetNode, determOnly = TRUE)
                my_calcCoeffAndOffset <- calcCoeffAndOffset()
            })
            
            ## make a nodeFunctionList of length=1, to hold the targetNode nodeFunction
            functionBody$addCode({
                targetNode_nodeFunctionList <- nimbleFunctionList(NF_VIRTUAL)
                targetNode_nodeFunctionList[[1]] <- model$nodeFunctions[[targetNode]]
            }, list(NF_VIRTUAL = as.name(paste0('node_stoch_', prior))))
            
            ## create lists of dependent node names, and nodeFunctions
            for(distName in dependentDistNames) {
                functionBody$addCode({
                    DEP_NODENAMES <- control$DEP_CONTROL_NAME
                    DEP_NODEFUNCTIONS <- nimbleFunctionList(NF_VIRTUAL)
                    for(i in seq_along(DEP_NODENAMES)) {
                        DEP_NODEFUNCTIONS[[i]] <- model$nodeFunctions[[DEP_NODENAMES[i]]]
                    }
                },
                list(DEP_CONTROL_NAME  = as.name(paste0('dependents_', distName)),
                     DEP_NODENAMES     = as.name(paste0('dependents_', distName, '_nodeNames')),
                     DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')),
                     NF_VIRTUAL        = as.name(paste0('node_stoch_', distName)))
                )
            }
            
            functionDef <- quote(function(model, mvSaved, control) {})
            functionDef[[3]] <- functionBody$getCode()
            functionDef[[4]] <- NULL   ## removes the 'scrref' attribute
            return(functionDef)
        },
        
        genRunFunction = function() {
            functionBody <- codeBlockClass()
            
            ## only if we're verifying conjugate posterior distributions: get initial targetValue, and modelLogProb -- getLogProb(model, calcNodes)
            if(nimbleOptions$verifyConjugatePosteriors) {
                functionBody$addCode({ modelLogProb0 <- getLogProb(model, calcNodes)
                                       origValue <- model[[targetNode]] })
            }
            
            addPosteriorQuantitiesGenerationCode(functionBody)    ## adds code to generate the quantities prior_xxx, and contribution_xxx
            
            ## generate new value, store, calculate, copy, etc...
            functionBody$addCode({
                newValue <- RPOSTERIORCALL
                model[[targetNode]] <<- newValue
                calculate(model, calcNodes)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
            }, list(RPOSTERIORCALL = posteriorObject$rCallExpr))
            
            ## only if we're verifying conjugate posterior distributions: figure out if conjugate posterior distribution is correct
            if(nimbleOptions$verifyConjugatePosteriors) {
                functionBody$addCode({modelLogProb1 <- getLogProb(model, calcNodes)
                                      posteriorLogDensity0 <- DPOSTERIORCALL_ORIG
                                      posteriorLogDensity1 <- DPOSTERIORCALL_NEW
                                      posteriorVerification <- modelLogProb0 - posteriorLogDensity0 - modelLogProb1 + posteriorLogDensity1
                                      if(abs(posteriorVerification) > 1e-10)     { nimPrint('conjugate posterior density appears to be wrong') }
                }, list(DPOSTERIORCALL_ORIG = eval(substitute(substitute(expr, list(VALUE=quote(origValue))), list(expr=posteriorObject$dCallExpr))),
                        DPOSTERIORCALL_NEW  = eval(substitute(substitute(expr, list(VALUE=quote(newValue))),  list(expr=posteriorObject$dCallExpr)))))
            }
            
            functionDef <- quote(function() {})
            functionDef[[3]] <- functionBody$getCode()
            functionDef[[4]] <- NULL   ## removes the 'scrref' attribute
            return(functionDef)
        },
        
        genGetPosteriorLogDensityFunction = function() {
            functionBody <- codeBlockClass()
            
            addPosteriorQuantitiesGenerationCode(functionBody)    ## adds code to generate the quantities prior_xxx, and contribution_xxx
            
            ## calculate and return the (log)density for the current value or targetNode
            functionBody$addCode({targetValue <- model[[targetNode]]
                                  posteriorLogDensity <- DPOSTERIORCALL
                                  returnType(double())
                                  return(posteriorLogDensity)
            }, list(DPOSTERIORCALL = eval(substitute(substitute(expr, list(VALUE=quote(targetValue))), list(expr=posteriorObject$dCallExpr)))))
            
            functionDef <- quote(function() {})
            functionDef[[3]] <- functionBody$getCode()
            functionDef[[4]] <- NULL   ## removes the 'scrref' attribute
            return(functionDef)
        },
        
        addPosteriorQuantitiesGenerationCode = function(functionBody) {
            
            ## get current value of prior parameters which appear in the posterior expression
            for(priorParam in posteriorObject$neededPriorParams) {
                functionBody$addCode(PRIOR_PARAM_VAR <- nfMethod(targetNode_nodeFunctionList[[1]], GET_PARAM_NAME)(),
                                     list(PRIOR_PARAM_VAR = as.name(paste0('prior_', priorParam)),
                                          GET_PARAM_NAME  =         paste0('get_', priorParam)))
            }
            
            ## get current values of ALL dependent nodes.  stated without proof, ALL values seem to be needed in the posterior
            for(distName in dependentDistNames) {
                functionBody$addCode({ declare(DEP_VALUES_VAR, double(1, length(DEP_NODEFUNCTIONS)))      ## DECLARE() statement
                                       getValues(DEP_VALUES_VAR, model, DEP_NODENAMES) },
                                     list(DEP_VALUES_VAR    = as.name(paste0('dependents_', distName, '_values')),
                                          DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')),
                                          DEP_NODENAMES     = as.name(paste0('dependents_', distName, '_nodeNames'))))
            }
            
            ## get values of dependent node parameters, only which we need for the posterior expression
            for(distName in dependentDistNames) {
                neededParams <- dependents[[distName]]$neededParamsForPosterior
                if(length(neededParams) == 0)     next
                forLoopBody <- codeBlockClass()
                for(param in neededParams) {
                    functionBody$addCode(declare(DEP_PARAM_VAR, double(1, length(DEP_NODEFUNCTIONS))),                               ## DECLARE() statement
                                         list(DEP_PARAM_VAR     = as.name(paste0('dependents_', distName, '_', param)),              ## DECLARE() statement
                                              DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions'))))       ## DECLARE() statement
                    forLoopBody$addCode(DEP_PARAM_VAR[i] <- nfMethod(DEP_NODEFUNCTIONS[[i]], GET_PARAM_NAME)(),
                                        list(DEP_PARAM_VAR     = as.name(paste0('dependents_', distName, '_', param)),
                                             DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')),
                                             GET_PARAM_NAME    =         paste0('get_', param)))
                }
                functionBody$addCode(for(i in seq_along(DEP_NODEFUNCTIONS)) FORLOOPBODY,
                                     list(DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')),
                                          FORLOOPBODY       = forLoopBody$getCode()))
            }
            
            ## only do the following if we need the linearity values 'coeff' or 'offset'
            if(needsLinearityCheck) {
                ## set new value into targetNode, and calculate all deterministic dependents
                functionBody$addCode({
                    x1 <- model[[targetNode]]
                    delta <- 1
                    x2 <- x1 + delta
                    model[[targetNode]] <<- x2
                    lp <- calculate(model, targetNode)
                    while(is.nan(lp)) {
                        delta <- delta * -1/2
                        x2 <- x1 + delta
                        model[[targetNode]] <<- x2
                        lp <- calculate(model, targetNode)
                    }
                    calculate(model, calcNodesDeterm)    ## flow the new value of targetNode down through all dependent deterministic nodes
                })
                
                for(distName in dependentDistNames) {
                    if(!dependents[[distName]]$needsCoeff && !dependents[[distName]]$needsOffset)     next
                    forLoopBody <- codeBlockClass()
                    forLoopBody$addCode({
                        y1 <- DEP_PARAM_VAR[i]
                        y2 <- nfMethod(DEP_NODEFUNCTIONS[[i]], GET_PARAM_NAME)()
                        coeffAndOffset <- my_calcCoeffAndOffset(x1, x2, y1, y2)
                    },
                    list(DEP_PARAM_VAR     = as.name(paste0('dependents_', distName, '_', dependents[[distName]]$param)),
                         DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')),
                         GET_PARAM_NAME    =         paste0('get_', dependents[[distName]]$param))
                    )
                    if(dependents[[distName]]$needsCoeff) {
                        functionBody$addCode(declare(DEP_COEFF_VAR, double(1, length(DEP_NODEFUNCTIONS))),                               ## DECLARE() statement
                                             list(DEP_COEFF_VAR     = as.name(paste0('dependents_', distName, '_coeff')),                ## DECLARE() statement
                                                  DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions'))))       ## DECLARE() statement
                        forLoopBody$addCode(DEP_COEFF_VAR[i] <- coeffAndOffset[1],
                                            list(DEP_COEFF_VAR = as.name(paste0('dependents_', distName, '_coeff'))))
                    }
                    if(dependents[[distName]]$needsOffset) {
                        functionBody$addCode(declare(DEP_OFFSET_VAR, double(1, length(DEP_NODEFUNCTIONS))),                              ## DECLARE() statement
                                             list(DEP_OFFSET_VAR    = as.name(paste0('dependents_', distName, '_offset')),               ## DECLARE() statement
                                                  DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions'))))       ## DECLARE() statement
                        forLoopBody$addCode(DEP_OFFSET_VAR[i] <- coeffAndOffset[2],
                                            list(DEP_OFFSET_VAR = as.name(paste0('dependents_', distName, '_offset'))))
                    }
                    functionBody$addCode(for(i in seq_along(DEP_NODEFUNCTIONS)) FORLOOPBODY,
                                         list(DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')),
                                              FORLOOPBODY       = forLoopBody$getCode()))
                }
            }
            
            ## set values for each contribution to the posterior expression, e.g., 'contribution_scale'
            for(contributionName in posteriorObject$neededContributionNames) {
                functionBody$addCode(CONTRIB_NAME <- 0,
                                     list(CONTRIB_NAME = as.name(contributionName)))
            }
            for(distName in dependentDistNames) {
                if(!any(posteriorObject$neededContributionNames %in% dependents[[distName]]$contributionNames))     next
                depParamsAvailable <- dependents[[distName]]$neededParamsForPosterior
                subList <- lapply(depParamsAvailable, function(param) parse(text=paste0('dependents_', distName, '_', param, '[i]'))[[1]])
                names(subList) <- depParamsAvailable
                subList$value  <- parse(text=paste0('dependents_', distName, '_values[i]'))[[1]]
                subList$coeff  <- parse(text=paste0('dependents_', distName, '_coeff[i]'))[[1]]
                subList$offset <- parse(text=paste0('dependents_', distName, '_offset[i]'))[[1]]
                forLoopBody <- codeBlockClass()
                for(contributionName in posteriorObject$neededContributionNames) {
                    if(!(contributionName %in% dependents[[distName]]$contributionNames))     next
                    contributionExpr <- eval(substitute(substitute(EXPR, subList), list(EXPR=dependents[[distName]]$contributionExprs[[contributionName]])))
                    forLoopBody$addCode(CONTRIB_NAME <- CONTRIB_NAME + CONTRIB_EXPR,
                                        list(CONTRIB_NAME = as.name(contributionName),
                                             CONTRIB_EXPR = contributionExpr)
                    )
                }
                functionBody$addCode(for(i in seq_along(DEP_NODEFUNCTIONS)) FORLOOPBODY,
                                     list(DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')),
                                          FORLOOPBODY       = forLoopBody$getCode()))
            }
        }
    )
)

dependentClass <- setRefClass(
    Class = 'dependentClass',
    fields = list(
        distribution = 				'ANY',   ## the name of the (dependent) sampling distribution, e.g. 'dnorm'
        param = 					'ANY', 	 ## the name of the sampling distribution parameter in which targetNode must appear
        link = 						'ANY',   ## the link: 'linear', 'multiplicative', or 'identity'
        contributionExprs = 		'ANY', 	 ## a (named) list of expressions, giving the (additive) contribution to any parameters of the posterior. names correspond to variables in the posterior expressions
        contributionNames = 		'ANY', 	 ## names of the contributions to the parameters of the posterior distribution.  same as names(posteriorExprs)
        neededParamsForPosterior = 	'ANY', 	 ## names of all parameters appearing in the posteriorExprs
        needsCoeff = 				'ANY', 	 ## logical, whether 'coeff' appears anywhere in the posteriorExprs
        needsOffset = 				'ANY' 	 ## logical, whether 'offset' appears anywhere in the posteriorExprs
    ),
    methods = list(
        initialize = function(depInfoList, depDistName) {
        	contributionExprs <<- list()
            distribution <<- cc_makeDDistributionName(depDistName)
            param <<- depInfoList$param
            link <<- depInfoList$link
            initialize_contributionExprs(depInfoList)
            initialize_neededParamsForPosterior()
        },
        initialize_contributionExprs = function(depInfoList) {
            depInfoList[c('param', 'link')] <- NULL
            depInfoList <- lapply(depInfoList, function(di) parse(text=di)[[1]])
            contributionExprs <<- depInfoList
            contributionNames <<- names(depInfoList)
        },
        initialize_neededParamsForPosterior = function() {
            posteriorVars <- unique(unlist(lapply(contributionExprs, all.vars)))
            neededParamsForPosterior <<- unique(c(param, posteriorVars[!(posteriorVars %in% c('value', 'coeff', 'offset'))]))
            needsCoeff  <<- 'coeff'  %in% posteriorVars
            needsOffset <<- 'offset' %in% posteriorVars
        }
    )
)

posteriorClass <- setRefClass(
    Class = 'posteriorClass',
    fields = list(
        posteriorExpr = 			'ANY',   ## the full, parsed, posterior distribution expression, e.g. dnorm(mean = prior_mean + ..., sd = ...)
        rDistribution = 			'ANY', 	 ## the *R* name of the posterior distribution, e.g. 'rnorm'
        dDistribution = 			'ANY', 	 ## the *R* name of the posterior density distribution, e.g. 'dnorm'
        argumentExprs = 			'ANY', 	 ## (named) list of expressions for each argument to the posterior distribution. names are the posterior distribution argument names
        argumentNames = 			'ANY',   ## character vector of the argument names to the posterior distribution.  same as: names(argumentExprs)
        rCallExpr = 				'ANY',   ## the actual 'rnorm(1, ...)' call, which will be substituted into the conjugate sampler function
        dCallExpr = 				'ANY',   ## the 'dnorm(value, ...)' call, which can be used to get values of the posterior density
        neededPriorParams = 		'ANY',   ## the names of any prior parameters (e.g., 'mean') which appear in the posterior expression as 'prior_mean'
        neededContributionNames = 	'ANY' 	 ## the names of contributions from dependent nodes, such as 'contribution_scale'
    ),
    methods = list(
        initialize = function(posteriorText) {
            posteriorExpr <<- parse(text = posteriorText)[[1]]
            rDistribution <<- cc_makeRDistributionName(as.character(posteriorExpr[[1]]))
            dDistribution <<- cc_makeDDistributionName(as.character(posteriorExpr[[1]]))
            argumentExprs <<- as.list(posteriorExpr)[-1]
            argumentNames <<- names(argumentExprs)
            rCallExpr <<- as.call(c(as.name(rDistribution), 1, argumentExprs))
            dCallExpr <<- as.call(c(as.name(dDistribution), quote(VALUE), argumentExprs, log = 1))
            posteriorVars <- all.vars(rCallExpr)
            neededPriorParams <<- gsub('^prior_', '', posteriorVars[grepl('^prior_', posteriorVars)])
            neededContributionNames <<- posteriorVars[grepl('^contribution_', posteriorVars)]
        }
    )
)

##############################################################################################
##############################################################################################
## utility functions
##############################################################################################
##############################################################################################


cc_makeSamplerTypeName       <- function(distName)     return(paste0('conjugate_', distName))    ## 'norm' --> 'conjugate_norm'
cc_makeConjugateSamplerName  <- function(samplerType)  return(paste0('sampler_', samplerType))   ## 'conjugate_norm' --> 'sampler_conjugate_norm'
cc_makeDDistributionName     <- function(distName)     return(paste0('d', distName))             ## 'norm' --> 'dnorm'
cc_makeRDistributionName     <- function(distName)     return(paste0('r', distName))             ## 'norm' --> 'rnorm'

## returns the text for the distribution of a stochastic node, e.g., 'dnorm'
cc_getNodeDistributionText <- function(model, node)     return(model$getNodeInfo()[[node]]$getDistribution())

## returns the expr corresponding to 'param' in the distribution of 'node'
cc_getNodeParamExpr <- function(model, node, param)     return(model$getNodeInfo()[[node]]$getParamExpr(param))

## returns NULL if param is not a parameter of the distribution for node, even after checking for re-parametizations
## returns list(expr = ..., paramFound = ...), giving the expression for the parameter, and the name of the actual parameter in which is was found.
# cc_findParamExpr <- function(model, node, param) {     -- obsolete?? (DT)
#     ## lord help anyone who tries to make sense of this code
#     paramsList <- as.list(cc_getNodeValueExpr(model, node)[-1])       # extracts the list of all parameters, for the distribution of node
#     if(param %in% names(paramsList))  return(list(expr = paramsList[[param]], paramFound = param))
#     
#     ## now, check for a possible re-parameterization
#     if(is.null(reparameterizationsList[[cc_getNodeDistributionText(model, node)]]))     return(NULL)   ## no possible re-parameterizations
#     reparamInfo <- reparameterizationsList[[cc_getNodeDistributionText(model, node)]]
#     for(i in seq_along(reparamInfo)) {
#         reparamExpr <- reparamInfo[[i]]
#         if(!cc_nodeInExpr(param, reparamExpr))   next   ## this reparameterization doesn't have 'param' in it
#         distParam <- names(reparamInfo)[i]
#         if(!(distParam %in% names(paramsList)))  next   ## our distribution for node doesn't have the right parameter
#         depNodeParamExpr <- paramsList[[distParam]]
#         depNodeParamExpr <- cc_expandDetermNodesInExpr(depNodeParamExpr, model)
#         parseTreeResult <- cc_comparePTexpressions(reparamExpr, depNodeParamExpr, param)
#         if(is.logical(parseTreeResult) && parseTreeResult == TRUE)   stop('this case should never occur....')
#         if(is.logical(parseTreeResult) && parseTreeResult == FALSE)   next
#         return(list(expr = parseTreeResult, paramFound = distParam))
#     }
#     return(NULL)   ## no suitable reparameterization found; return NULL
# }

# cc_comparePTexpressions <- function(templateExpr, actualExpr, param, foundExpr = NULL) {     -- obsolete?? (DT)
#     ## likewise, lord help anyone who tries to make sense of this code, either
#     if(class(templateExpr) == 'name') {
#         if(templateExpr == param) {
#             if(is.null(foundExpr))     return(actualExpr)
#             return(FALSE)
#         }
#         if(identical(templateExpr, actualExpr))  return(TRUE) else return(FALSE)
#     }
#     if(class(templateExpr) %in% c('numeric', 'integer', 'logical')) {
#         if(identical(templateExpr, actualExpr))  return(TRUE) else return(FALSE)
#     }
#     if(length(templateExpr) != length(actualExpr)) return(FALSE)
#     for(i in seq_along(templateExpr)) {
#         result <- cc_comparePTexpressions(templateExpr[[i]], actualExpr[[i]], param, foundExpr)
#         if(is.logical(result) && result == FALSE) return(FALSE)
#         if(is.logical(result) && result == TRUE) next
#         foundExpr <- result
#     }
#     if(is.null(foundExpr))  return(TRUE)
#     return(foundExpr)
# }

##  returns the entire RHS valueExpr for 'node'
cc_getNodeValueExpr <- function(model, node) {
    if(!any(node == model$getNodeNames()))   stop(paste0('node not present in model: ', node))
    return(model$getNodeInfo()[[node]]$getValueExpr())
}

## expands all deterministic nodes in expr, to create a single expression with only stochastic nodes
cc_expandDetermNodesInExpr <- function(expr, model) {
    if(is.numeric(expr))     return(expr)     # return numeric
    if(is.name(expr)    ||    (is.call(expr) && (expr[[1]] == '['))) {    # expr is a name, or an indexed name
        exprText <- deparse(expr)
        if(!any(exprText == model$getMaps('nodeNames') ) )          stop('name found which isn\'t a model node')
        if( any(exprText == model$getMaps('nodeNamesStoch') ) )        return(expr)      # return stochastic nodes
        if( any(exprText == model$getMaps('nodeNamesLHSinferred') ) )  return(expr)      # return LHS nodes inferred from a multivariate stochastic distribution
        if( any(exprText == model$getMaps('nodeNamesRHSonly') ) )   stop('something wrong with model; possible failure to specify constants = ..., for a RHS-only node')
        if(!any(exprText == model$getMaps('nodeNamesDeterm') ) )    stop('something went wrong, possibly bad model specification')
        return(cc_expandDetermNodesInExpr(expr=cc_getNodeValueExpr(model,node=exprText), model))   # precess and return the value expression for this deterministic node
    }
    if(is.call(expr)) {
        for(i in seq_along(expr)[-1])    expr[[i]] <- cc_expandDetermNodesInExpr(expr[[i]], model)
        return(expr) }
    stop(paste0('something went wrong processing: ', deparse(expr)))
}

## verifies that 'link' is satisfied by the results of linearityCheck
cc_linkCheck <- function(linearityCheck, link) {
    if(!(link %in% c('identity', 'multiplicative', 'linear')))    stop(paste0('unknown link: \'', link, '\''))
    if(is.null(linearityCheck))    return(FALSE)
    offset <- linearityCheck$offset
    scale  <- linearityCheck$scale
    if(link == 'identity'       && offset == 0 && scale == 1)     return(TRUE)
    if(link == 'multiplicative' && offset == 0)                   return(TRUE)
    if(link == 'linear')                                          return(TRUE)
    return(FALSE)
}

## checks the parameter expressions in the stochastic distribution of depNode
## returns FALSE if we find 'targetNode' in ***more than one*** of these expressions
cc_otherParamsCheck <- function(model, depNode, targetNode) {
    paramsList <- as.list(cc_getNodeValueExpr(model, depNode)[-1])       # extracts the list of all parameters, for the distribution of depNode
    timesFound <- 0   ## for success, we'll find targetNode in only *one* parameter expression
    for(i in seq_along(paramsList)) {
        expr <- cc_expandDetermNodesInExpr(paramsList[[i]], model)
        if(cc_nodeInExpr(targetNode, expr))     { timesFound <- timesFound + 1 }    ## we found 'targetNode'
    }
    if(timesFound == 0)     stop('something went wrong; targetNode not found in any parameter expressions')
    if(timesFound == 1)     return(TRUE)
    if(timesFound  > 1)     return(FALSE)
}

## determines whether node appears anywhere in expr
cc_nodeInExpr <- function(node, expr) { return(node %in% cc_getNodesInExpr(expr)) }

## determines which nodes apppear in an expression
cc_getNodesInExpr <- function(expr) {
    if(is.numeric(expr)) return(character(0))   ## expr is numeric
    if(is.logical(expr)) return(character(0))   ## expr is logical
    if(is.name(expr) || (is.call(expr) && (expr[[1]] == '['))) return(deparse(expr))   ## expr is a node name
    if(is.call(expr)) return(unlist(lapply(expr[-1], cc_getNodesInExpr)))   ## expr is some general call
    stop(paste0('something went wrong processing: ', deparse(expr)))
}



##############################################################################################
##############################################################################################
## general, flexible, check for linear relationship
## returns list(offset, scale), or NULL
##############################################################################################
##############################################################################################

cc_checkLinearity <- function(expr, targetNode) {
    
    ## targetNode doesn't appear in expr
    if(!cc_nodeInExpr(targetNode, expr)) {
        if(is.call(expr) && expr[[1]] == '(') return(cc_checkLinearity(expr[[2]], targetNode))
        return(list(offset = expr, scale = 0))
    }
    
    ## expr is exactly the targetNode
    if(identical(targetNode, deparse(expr)))
        return(list(offset = 0, scale = 1))
    
    if(!is.call(expr))   stop('something went wrong')
    
    ## process the expression contents of the parentheses
    if(expr[[1]] == '(')
        return(cc_checkLinearity(expr[[2]], targetNode))
    
    ## minus sign: change to a plus sign, and invert the sign of the RHS
    if(expr[[1]] == '-') {
        if(length(expr) == 3) {
            checkLinearityLHS <- cc_checkLinearity(expr[[2]], targetNode)
            checkLinearityRHS <- cc_checkLinearity(expr[[3]], targetNode)
            if(is.null(checkLinearityLHS) || is.null(checkLinearityRHS)) return(NULL)
            return(list(offset = cc_combineExprsSubtraction(checkLinearityLHS$offset, checkLinearityRHS$offset),
                        scale  = cc_combineExprsSubtraction(checkLinearityLHS$scale,  checkLinearityRHS$scale)))
        }
        if(length(expr) == 2) {
            checkLin <- cc_checkLinearity(expr[[2]], targetNode)
            if(is.null(checkLin)) return(NULL)
            return(list(offset = cc_negateExpr(checkLin$offset),
                        scale  = cc_negateExpr(checkLin$scale)))
        }
        stop('something went wrong')
    }
    
    if(expr[[1]] == '+') {
        checkLinearityLHS <- cc_checkLinearity(expr[[2]], targetNode)
        checkLinearityRHS <- cc_checkLinearity(expr[[3]], targetNode)
        if(is.null(checkLinearityLHS) || is.null(checkLinearityRHS)) return(NULL)
        return(list(offset = cc_combineExprsAddition(checkLinearityLHS$offset, checkLinearityRHS$offset),
                    scale  = cc_combineExprsAddition(checkLinearityLHS$scale,  checkLinearityRHS$scale)))
    }
    
    if(expr[[1]] == '*') {
        if(cc_nodeInExpr(targetNode, expr[[2]]) && cc_nodeInExpr(targetNode, expr[[3]])) return(NULL)
        checkLinearityLHS <- cc_checkLinearity(expr[[2]], targetNode)
        checkLinearityRHS <- cc_checkLinearity(expr[[3]], targetNode)
        if(is.null(checkLinearityLHS) || is.null(checkLinearityRHS)) return(NULL)
        if((checkLinearityLHS$scale != 0) && (checkLinearityRHS$scale != 0)) stop('something went wrong')
        if(checkLinearityLHS$scale != 0) {
            return(list(offset = cc_combineExprsMultiplication(checkLinearityLHS$offset, checkLinearityRHS$offset),
                        scale  = cc_combineExprsMultiplication(checkLinearityLHS$scale,  checkLinearityRHS$offset)))
        }
        if(checkLinearityRHS$scale != 0) {
            return(list(offset = cc_combineExprsMultiplication(checkLinearityLHS$offset, checkLinearityRHS$offset),
                        scale  = cc_combineExprsMultiplication(checkLinearityLHS$offset, checkLinearityRHS$scale)))
        }
        stop('something went wrong')
    }
    
    if(expr[[1]] == '/') {
        if(cc_nodeInExpr(targetNode, expr[[3]])) return(NULL)
        checkLinearityLHS <- cc_checkLinearity(expr[[2]], targetNode)
        checkLinearityRHS <- cc_checkLinearity(expr[[3]], targetNode)
        if(is.null(checkLinearityLHS) || is.null(checkLinearityRHS)) return(NULL)
        if(checkLinearityLHS$scale == 0) stop('something went wrong')
        if(checkLinearityRHS$scale != 0) stop('something went wrong')
        return(list(offset = cc_combineExprsDivision(checkLinearityLHS$offset, checkLinearityRHS$offset),
                    scale  = cc_combineExprsDivision(checkLinearityLHS$scale,  checkLinearityRHS$offset)))
    }
    
    ## returns not conjugate (NULL) if we don't recognize the call (expr[[1]])
    return(NULL)
}

cc_negateExpr <- function(expr) {
    if(expr == 0) return(0)
    if(is.numeric(expr)) return(-1*expr)
    return(substitute(-EXPR, list(EXPR = expr)))
}

cc_combineExprsAddition <- function(expr1, expr2) {
    if((expr1 == 0) && (expr2 == 0)) return(0)
    if((expr1 != 0) && (expr2 == 0)) return(expr1)
    if((expr1 == 0) && (expr2 != 0)) return(expr2)
    if(is.numeric(expr1) && is.numeric(expr2)) return(expr1 + expr2)
    return(substitute(EXPR1 + EXPR2, list(EXPR1=expr1, EXPR2=expr2)))
}

cc_combineExprsSubtraction <- function(expr1, expr2) {
    if((expr1 == 0) && (expr2 == 0)) return(0)
    if((expr1 != 0) && (expr2 == 0)) return(expr1)
    if((expr1 == 0) && (expr2 != 0)) return(cc_negateExpr(expr2))
    if(is.numeric(expr1) && is.numeric(expr2)) return(expr1 - expr2)
    return(substitute(EXPR1 - EXPR2, list(EXPR1=expr1, EXPR2=expr2)))
}

cc_combineExprsMultiplication <- function(expr1, expr2) {
    if((expr1 == 0) || (expr2 == 0)) return(0)
    if((expr1 == 1) && (expr2 == 1)) return(1)
    if((expr1 != 1) && (expr2 == 1)) return(expr1)
    if((expr1 == 1) && (expr2 != 1)) return(expr2)
    if(is.numeric(expr1) && is.numeric(expr2)) return(expr1 * expr2)
    return(substitute(EXPR1 * EXPR2, list(EXPR1=expr1, EXPR2=expr2)))
}

cc_combineExprsDivision <- function(expr1, expr2) {
    if(                (expr2 == 0)) stop('error, division by 0')
    if(                (expr2 == 1)) return(expr1)
    if((expr1 == 0)                ) return(0)
    if(is.numeric(expr1) && is.numeric(expr2)) return(expr1 / expr2)
    return(substitute(EXPR1 / EXPR2, list(EXPR1=expr1, EXPR2=expr2)))
}




##############################################################################################
##############################################################################################
## create object: conjugacyRelationshipsObject
## also, generate all conjugate sampler nimbleFunctions
##############################################################################################
##############################################################################################


conjugacyRelationshipsObject <- conjugacyRelationshipsClass(conjugacyRelationshipsInputList)

conjugateSamplerDefinitions <- conjugacyRelationshipsObject$generateConjugateSamplerDefinitions()
createNamedObjectsFromList(conjugateSamplerDefinitions)
# createNamedObjectsFromList(conjugateSamplerDefinitions, writeToFile = 'TEMP_conjugateSamplerDefinitions.R')











