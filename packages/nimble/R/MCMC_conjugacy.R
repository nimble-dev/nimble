



conjugacyRelationshipsInputList <- list(
    
    ## beta
    list(prior = 'dbeta',
         link = 'identity',
         dependents = list(
             dbern   = list(param = 'prob', contribution_shape1 = 'value', contribution_shape2 = '1 - value'   ),
             dbin    = list(param = 'prob', contribution_shape1 = 'value', contribution_shape2 = 'size - value'),
             dnegbin = list(param = 'prob', contribution_shape1 = 'size',  contribution_shape2 = 'value'       )),
         posterior = 'dbeta(shape1 = prior_shape1 + contribution_shape1,
                            shape2 = prior_shape2 + contribution_shape2)'),
    
    ## Dirichlet - added by CJP 1/14/15
    ## at moment can't do dcat because of limitations of conjugacy processing relying on nimble code
    list(prior = 'ddirch',
         link = 'identity',
         dependents = list(
             dmulti    = list(param = 'prob', contribution_alpha = 'value')),
        #     dcat      = list(param = 'prob', contribution_alpha = as.numeric((1:length(prob)) == value)'),
        #     dcat      = list(param = 'prob', contribution_alpha = {tmp = rep(0,length(prob)); tmp[value]=1; tmp}')),
          posterior = 'ddirch(alpha = prior_alpha + contribution_alpha)'), 
    
    ## gamma
    list(prior = 'dgamma',
         link = 'multiplicative',
         dependents = list(
             dpois  = list(param = 'lambda', contribution_shape = 'value', contribution_rate = 'coeff'                           ),
             dnorm  = list(param = 'tau',    contribution_shape = '1/2',   contribution_rate = 'coeff/2 * (value-mean)^2'        ),
             dlnorm = list(param = 'tau',    contribution_shape = '1/2',   contribution_rate = 'coeff/2 * (log(value)-meanlog)^2'),
             dgamma = list(param = 'rate',   contribution_shape = 'shape', contribution_rate = 'coeff   * value'                 ),
             dexp   = list(param = 'rate',   contribution_shape = '1',     contribution_rate = 'coeff   * value'                 )),
             ## ddexp  = list(param = 'rate',   contribution_shape = '1',     contribution_rate = 'coeff   * abs(value-location)'   )
             ## dpar = list(...)    ## need to figure this out
         posterior = 'dgamma(shape = prior_shape + contribution_shape,
                             scale = 1 / (prior_rate + contribution_rate))'),
    
    ## normal
    list(prior = 'dnorm',
         link = 'linear',
         dependents = list(
             dnorm  = list(param = 'mean',    contribution_mean = 'coeff * (value-offset) * tau',      contribution_tau = 'coeff^2 * tau'),
             dlnorm = list(param = 'meanlog', contribution_mean = 'coeff * (log(value)-offset) * tau', contribution_tau = 'coeff^2 * tau')),
         posterior = 'dnorm(mean = (prior_mean*prior_tau + contribution_mean) / (prior_tau + contribution_tau),
                            sd   = (prior_tau + contribution_tau)^(-0.5))'),
    
    ## pareto
    # list(prior = 'dpar',      ##### waiting for dpar() distribution
    #      link = 'multiplicative',
    #      dependents = list(
    #          dunif = list(param = 'max', contribution_alpha = '1', contribution_not_used = 'coeff'),
    #          dpar  = list(param = 'c',   contribution_alpha = '-alpha')),
    #      posterior = 'dpar(alpha = prior_alpha + contribution_alpha,
    #                        c     = max(prior_c, max(dependents_dunif_values/dependents_dunif_coeff)))'),
    #####
    
    ## multivariate-normal
    list(prior = 'dmnorm',
         link = 'linear',
         dependents = list(
             dmnorm = list(param = 'mean', contribution_mean = 't(coeff) %*% prec %*% asCol(value-offset)', contribution_prec = 't(coeff) %*% prec %*% coeff')),
         posterior = 'dmnorm_chol(mean       = (inverse( (prior_prec + contribution_prec) ) %*% (prior_prec %*% asCol(prior_mean) + contribution_mean))[,1],
                                  chol       = chol(prior_prec + contribution_prec),
                                  prec_param = 1)'),

    ## wishart
    list(prior = 'dwish',
         link = 'linear',
         dependents = list(
             dmnorm = list(param = 'prec', contribution_R = 'asCol(value-mean) %*% asRow(value-mean) %*% coeff', contribution_df = '1')),
         posterior = 'dwish_chol(chol        = chol(prior_R + contribution_R),
                                 df          = prior_df + contribution_df,
                                 scale_param = 0)')

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
            gIDs_4_checking <- numeric(0)
            try(gIDs_4_checking <- model$modelDef$nodeName2GraphIDs(targetNode), silent = TRUE)
            if(length(gIDs_4_checking) == 0)       stop('checking conjugacy of a node not in model')
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

setMethod('[[',   'conjugacyRelationshipsClass',
          function(x, i) {
              return(x$conjugacys[[i]])
          }
)

conjugacyClass <- setRefClass(
    Class = 'conjugacyClass',
    fields = list(
        samplerType = 			'ANY', 		## name of the sampler for this conjugacy class, e.g. 'conjugate_dnorm'
        prior =					'ANY', 		## name of the prior distribution, e.g. 'dnorm'
        link =     				'ANY',      ## the link ('linear', 'multiplicative', or 'identity')
        dependents = 			'ANY', 		## (named) list of dependentClass objects, each contains conjugacy information specific to a particular sampling distribution (name is sampling distribution name)
        dependentDistNames = 	'ANY', 		## character vector of the names of all allowable dependent sampling distributions.  same as: names(dependents)
        posteriorObject = 		'ANY',   	## an object of posteriorClass
        needsLinearityCheck = 	'ANY', 		## logical specifying whether we need to do the linearity check; if the link is 'multiplicative' or 'linear'
        model = 				'ANY' 	    ## ONLY EXISTS TO PREVENT A WARNING for '<<-', in the code for generating the conjugate sampler function
    ),
    methods = list(
        initialize = function(cr) {
        	dependents <<- list()
            samplerType <<- cc_makeSamplerTypeName(cr$prior)
            prior <<- cr$prior
            link <<- cr$link
            initialize_addDependents(cr$dependents)
            needsLinearityCheck <<- link %in% c('multiplicative', 'linear')
            posteriorObject <<- posteriorClass(cr$posterior)
            model <<- NA
        },
        initialize_addDependents = function(depList) {
            for(i in seq_along(depList)) {
                dependents[[i]] <<- dependentClass(depList[[i]], names(depList)[i])
            }
            names(dependents) <<- names(depList)
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
                if(!cc_nodeInExpr(targetNode, linearityCheckExpr))                return(NULL)
                if(cc_vectorizedComponentCheck(targetNode, linearityCheckExpr))   return(NULL)   # if targetNode is vectorized, make sure non of it's components appear in expr
                linearityCheck <- cc_checkLinearity(linearityCheckExpr, targetNode)   # determines whether paramExpr is linear in targetNode
                if(!cc_linkCheck(linearityCheck, link))                           return(NULL)
                if(!cc_otherParamsCheck(model, depNode, targetNode))              return(NULL)   # ensure targetNode appears in only *one* depNode parameter expression
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
                targetNode      <- control$targetNode
                calcNodes       <- model$getDependencies(targetNode)
                calcNodesDeterm <- model$getDependencies(targetNode, determOnly = TRUE)
                ######### my_calcCoeffAndOffset <- calcCoeffAndOffset()   # no longer needed -DT
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
            
            ## if this conjugate sampler is for a multivariate node (i.e., nDim > 0), then we need to determine the size (d)
            if(distributions[[prior]]$types$value$nDim > 0) {
                functionBody$addCode(d <- model$getNodeInfo()[[targetNode]]$targetNodeIndexSizes[1])
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
                functionBody$addCode({ 
                					   modelLogProb0 <- getLogProb(model, calcNodes)
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
                                      if(abs(posteriorVerification) > 1e-8)     {
                                      nimPrint('conjugate posterior density appears to be wrong, off by ', posteriorVerification) }
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
            
            ## get values of all dependent nodes, and values of dependent node parameters (those needed for posterior expression)
            for(distName in dependentDistNames) {
                forLoopBody <- codeBlockClass()
                
                depNodeValueNdim <- distributions[[distName]]$types$value$nDim
                functionBody$addCode(declare(DEP_VALUES_VAR, double(DEP_VALUES_VAR_NDIM, DECLARE_SIZE)),                              ## DECLARE() statement
                                     list(DEP_VALUES_VAR         = as.name(paste0('dependents_', distName, '_values')),               ## DECLARE() statement
                                          DEP_VALUES_VAR_NDIM    = 1 + depNodeValueNdim,                                              ## DECLARE() statement
                                          DECLARE_SIZE           = makeDeclareSizeField(substitute(length(DEP_NODEFUNCTIONS), list(DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')))), depNodeValueNdim)))
                forLoopBody$addCode(DEP_VALUES_VAR_INDEXED <- nfMethod(DEP_NODEFUNCTIONS[[i]], 'get_value')(),
                                    list(DEP_VALUES_VAR_INDEXED = makeIndexedVariable(as.name(paste0('dependents_', distName, '_values')), depNodeValueNdim),
                                         DEP_NODEFUNCTIONS     = as.name(paste0('dependents_', distName, '_nodeFunctions'))))
                
                neededParams <- dependents[[distName]]$neededParamsForPosterior
                for(param in neededParams) {
                    depNodeParamNdim <- distributions[[distName]]$types[[param]]$nDim
                    functionBody$addCode(declare(DEP_PARAM_VAR, double(DEP_PARAM_VAR_NDIM, DECLARE_SIZE)),                            ## DECLARE() statement
                                         list(DEP_PARAM_VAR      = as.name(paste0('dependents_', distName, '_', param)),              ## DECLARE() statement
                                              DEP_PARAM_VAR_NDIM = 1 + depNodeParamNdim,                                              ## DECLARE() statement
                                              DECLARE_SIZE       = makeDeclareSizeField(substitute(length(DEP_NODEFUNCTIONS), list(DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')))), depNodeParamNdim)))
                    forLoopBody$addCode(DEP_PARAM_VAR_INDEXED <- nfMethod(DEP_NODEFUNCTIONS[[i]], GET_PARAM_NAME)(),
                                        list(DEP_PARAM_VAR_INDEXED = makeIndexedVariable(as.name(paste0('dependents_', distName, '_', param)), depNodeParamNdim),
                                             DEP_NODEFUNCTIONS     = as.name(paste0('dependents_', distName, '_nodeFunctions')),
                                             GET_PARAM_NAME        =         paste0('get_', param)))
                }
                
                functionBody$addCode(for(i in seq_along(DEP_NODEFUNCTIONS)) FORLOOPBODY,
                                     list(DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')),
                                          FORLOOPBODY       = forLoopBody$getCode()))
            }
            
            ## if we need to determine 'coeff' and/or 'offset'
            if(needsLinearityCheck) {
                targetNodeNdim <- distributions[[prior]]$types$value$nDim
                targetCoeffNdim <- switch(as.character(targetNodeNdim), `0`=0, `1`=2, `2`=2, stop())
                
                ## all the declare statements
                for(distName in dependentDistNames) {                                                           ## DECLARE() statement
                    functionBody$addCode({                                                                      ## DECLARE() statement
                        declare(DEP_OFFSET_VAR, double(DEP_OFFSET_VAR_NDIM, DECLARE_SIZE_OFFSET))               ## DECLARE() statement
                        declare(DEP_COEFF_VAR,  double(DEP_COEFF_VAR_NDIM,  DECLARE_SIZE_COEFF))                ## DECLARE() statement
                    }, list(DEP_OFFSET_VAR      = as.name(paste0('dependents_', distName, '_offset')),          ## DECLARE() statement
                            DEP_COEFF_VAR       = as.name(paste0('dependents_', distName, '_coeff')),           ## DECLARE() statement
                            DEP_OFFSET_VAR_NDIM = 1 + targetNodeNdim,                                           ## DECLARE() statement
                            DEP_COEFF_VAR_NDIM  = 1 + targetCoeffNdim,                                          ## DECLARE() statement
                            DECLARE_SIZE_OFFSET = makeDeclareSizeField(substitute(length(DEP_NODEFUNCTIONS), list(DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')))), targetNodeNdim),
                            DECLARE_SIZE_COEFF  = makeDeclareSizeField(substitute(length(DEP_NODEFUNCTIONS), list(DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')))), targetCoeffNdim)))
                }
                
                switch(as.character(targetNodeNdim),
                       `0` = {
                           functionBody$addCode({
                               model[[targetNode]] <<- 0
                               calculate(model, calcNodesDeterm)
                           })
                           for(distName in dependentDistNames) {
                               functionBody$addCode(
                                   for(i in seq_along(DEP_NODEFUNCTIONS)) {
                                       DEP_OFFSET_VAR[i] <- nfMethod(DEP_NODEFUNCTIONS[[i]], GET_PARAM_NAME)()
                                   }, list(DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')),
                                           DEP_OFFSET_VAR    = as.name(paste0('dependents_', distName, '_offset')),
                                           GET_PARAM_NAME    =         paste0('get_', dependents[[distName]]$param)))
                           }
                           functionBody$addCode({
                               model[[targetNode]] <<- 1
                               calculate(model, calcNodesDeterm)
                           })
                           for(distName in dependentDistNames) {
                               functionBody$addCode(
                                   for(i in seq_along(DEP_NODEFUNCTIONS)) {
                                       DEP_COEFF_VAR[i] <- nfMethod(DEP_NODEFUNCTIONS[[i]], GET_PARAM_NAME)() - DEP_OFFSET_VAR[i]
                                   }, list(DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')),
                                           DEP_COEFF_VAR     = as.name(paste0('dependents_', distName, '_coeff')),
                                           GET_PARAM_NAME    =         paste0('get_', dependents[[distName]]$param),
                                           DEP_OFFSET_VAR    = as.name(paste0('dependents_', distName, '_offset'))))
                           }
                       },
                       `1` = {
                           functionBody$addCode({
                               model[[targetNode]] <<- model[[targetNode]] * 0
                               calculate(model, calcNodesDeterm)
                           })
                           for(distName in dependentDistNames) {
                               functionBody$addCode(
                                   for(i in seq_along(DEP_NODEFUNCTIONS)) {
                                       DEP_OFFSET_VAR[i, 1:d] <- nfMethod(DEP_NODEFUNCTIONS[[i]], GET_PARAM_NAME)()
                                   }, list(DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')),
                                           DEP_OFFSET_VAR    = as.name(paste0('dependents_', distName, '_offset')),
                                           GET_PARAM_NAME    =         paste0('get_', dependents[[distName]]$param)))
                           }
                           forLoopBody <- codeBlockClass()
                           forLoopBody$addCode({
                               unitVector <- model[[targetNode]] * 0
                               unitVector[sizeIndex] <- 1
                               model[[targetNode]] <<- unitVector
                               calculate(model, calcNodesDeterm)
                           })
                           for(distName in dependentDistNames) {
                               forLoopBody$addCode(
                                   for(i in seq_along(DEP_NODEFUNCTIONS)) {
                                       DEP_COEFF_VAR[i, 1:d, sizeIndex] <- nfMethod(DEP_NODEFUNCTIONS[[i]], GET_PARAM_NAME)() - DEP_OFFSET_VAR[i, 1:d]
                                   }, list(DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')),
                                           DEP_COEFF_VAR     = as.name(paste0('dependents_', distName, '_coeff')),
                                           GET_PARAM_NAME    =         paste0('get_', dependents[[distName]]$param),
                                           DEP_OFFSET_VAR    = as.name(paste0('dependents_', distName, '_offset'))))
                           }
                           functionBody$addCode(for(sizeIndex in 1:d) FORLOOPBODY,
                                                list(FORLOOPBODY = forLoopBody$getCode()))
                       },
                       `2` = {
                           functionBody$addCode({
                               identityMatrix <- model[[targetNode]] * 0
                               for(sizeIndex in 1:d)   { identityMatrix[sizeIndex, sizeIndex] <- 1 }
                               model[[targetNode]] <<- identityMatrix   ## initially, propogate through X = I
                               calculate(model, calcNodesDeterm)
                           })
                           for(distName in dependentDistNames) {
                               functionBody$addCode(
                                   for(i in seq_along(DEP_NODEFUNCTIONS)) {
                                       DEP_OFFSET_VAR[i, 1:d, 1:d] <- nfMethod(DEP_NODEFUNCTIONS[[i]], GET_PARAM_NAME)()   ## DEP_OFFSET_VAR = A+B(I) = A+B
                                   }, list(DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')),
                                           DEP_OFFSET_VAR    = as.name(paste0('dependents_', distName, '_offset')),
                                           GET_PARAM_NAME    =         paste0('get_', dependents[[distName]]$param)))
                           }
                           functionBody$addCode({
                               model[[targetNode]] <<- identityMatrix * 2   ## now, propogate through X = 2I
                               calculate(model, calcNodesDeterm)
                           })
                           for(distName in dependentDistNames) {
                               functionBody$addCode(
                                   for(i in seq_along(DEP_NODEFUNCTIONS)) {
                                       DEP_COEFF_VAR[i, 1:d, 1:d] <- nfMethod(DEP_NODEFUNCTIONS[[i]], GET_PARAM_NAME)()   ## DEP_COEFF_VAR = A+B(2I) = A+2B
                                   }, list(DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')),
                                           DEP_COEFF_VAR     = as.name(paste0('dependents_', distName, '_coeff')),
                                           GET_PARAM_NAME    =         paste0('get_', dependents[[distName]]$param))
                               )
                           }
                           for(distName in dependentDistNames) {
                               functionBody$addCode(
                                   for(i in seq_along(DEP_NODEFUNCTIONS)) {
                                       DEP_COEFF_VAR[i, 1:d, 1:d] <- DEP_COEFF_VAR[i, 1:d, 1:d] - DEP_OFFSET_VAR[i, 1:d, 1:d]   ## now, DEP_COEFF_VAR = (A+2B)-(A+B) = B
                                       DEP_OFFSET_VAR[i, 1:d, 1:d] <- DEP_OFFSET_VAR[i, 1:d, 1:d] - DEP_COEFF_VAR[i, 1:d, 1:d]   ## now, DEP_OFFSET_VAR = (A+B)-(B) = A
                                   }, list(DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')),
                                           DEP_COEFF_VAR     = as.name(paste0('dependents_', distName, '_coeff')),
                                           DEP_OFFSET_VAR    = as.name(paste0('dependents_', distName, '_offset')))
                               )
                           }

                           ####### starting over above here
                           ## functionBody$addCode({
                           ##     model[[targetNode]] <<- model[[targetNode]] * 0
                           ##     calculate(model, calcNodesDeterm)
                           ## })
                           ## for(distName in dependentDistNames) {
                           ##     functionBody$addCode(
                           ##         for(i in seq_along(DEP_NODEFUNCTIONS)) {
                           ##             DEP_OFFSET_VAR[i, 1:d, 1:d] <- nfMethod(DEP_NODEFUNCTIONS[[i]], GET_PARAM_NAME)()
                           ##         }, list(DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')),
                           ##                 DEP_OFFSET_VAR    = as.name(paste0('dependents_', distName, '_offset')),
                           ##                 GET_PARAM_NAME    =         paste0('get_', dependents[[distName]]$param)))
                           ## }
                           ## functionBody$addCode({
                           ##     identityMatrix <- model[[targetNode]] * 0
                           ##     for(sizeIndex in 1:d)   { identityMatrix[sizeIndex, sizeIndex] <- 1 }
                           ##     model[[targetNode]] <<- identityMatrix
                           ##     calculate(model, calcNodesDeterm)
                           ## })
                           ## for(distName in dependentDistNames) {
                           ##     functionBody$addCode(
                           ##         for(i in seq_along(DEP_NODEFUNCTIONS)) {
                           ##             DEP_COEFF_VAR[i, 1:d, 1:d] <- nfMethod(DEP_NODEFUNCTIONS[[i]], GET_PARAM_NAME)() - DEP_OFFSET_VAR[i, 1:d, 1:d]
                           ##         }, list(DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')),
                           ##                 DEP_COEFF_VAR     = as.name(paste0('dependents_', distName, '_coeff')),
                           ##                 GET_PARAM_NAME    =         paste0('get_', dependents[[distName]]$param),
                           ##                 DEP_OFFSET_VAR    = as.name(paste0('dependents_', distName, '_offset')))
                           ##     )
                           ## }
                       },
                       stop()
                )
                
            } # end if(needsLinearityCheck)
            
            targetNodeNdim <- distributions[[prior]]$types$value$nDim
            targetCoeffNdim <- switch(as.character(targetNodeNdim), `0`=0, `1`=2, `2`=2, stop())
            
            functionBody$addCode(firstTime <- 1)
            
            for(distName in dependentDistNames) {
                if(!any(posteriorObject$neededContributionNames %in% dependents[[distName]]$contributionNames))     next
                depParamsAvailable <- dependents[[distName]]$neededParamsForPosterior
                subList <- lapply(depParamsAvailable, function(param) makeIndexedVariable(as.name(paste0('dependents_', distName, '_', param)), distributions[[distName]]$types[[param]]$nDim))
                names(subList) <- depParamsAvailable
                subList$value  <- makeIndexedVariable(as.name(paste0('dependents_', distName, '_values')), distributions[[distName]]$types$value$nDim)
                subList$offset <- makeIndexedVariable(as.name(paste0('dependents_', distName, '_offset')), targetNodeNdim)
                subList$coeff  <- makeIndexedVariable(as.name(paste0('dependents_', distName, '_coeff')),  targetCoeffNdim)
                forLoopBodyFirst <- codeBlockClass()
                forLoopBody      <- codeBlockClass()
                forLoopBodyFirst$addCode(firstTime <- 0)
                for(contributionName in posteriorObject$neededContributionNames) {
                    if(!(contributionName %in% dependents[[distName]]$contributionNames))     next
                    contributionExpr <- eval(substitute(substitute(EXPR, subList), list(EXPR=dependents[[distName]]$contributionExprs[[contributionName]])))
                    forLoopBodyFirst$addCode(CONTRIB_NAME <- CONTRIB_EXPR,
                                             list(CONTRIB_NAME = as.name(contributionName), CONTRIB_EXPR = contributionExpr))
                    forLoopBody$addCode(     CONTRIB_NAME <- CONTRIB_NAME + CONTRIB_EXPR,
                                             list(CONTRIB_NAME = as.name(contributionName), CONTRIB_EXPR = contributionExpr))
                }
                ifStatementBody <- codeBlockClass()
                ifStatementBody$addCode(
                    if(firstTime == 1) FORLOOPBODYFIRST else FORLOOPBODY,
                    list(FORLOOPBODYFIRST = forLoopBodyFirst$getCode(),
                         FORLOOPBODY      = forLoopBody$getCode()))
                functionBody$addCode(for(i in seq_along(DEP_NODEFUNCTIONS)) IFSTATEMENTBODY,
                                     list(DEP_NODEFUNCTIONS = as.name(paste0('dependents_', distName, '_nodeFunctions')),
                                          IFSTATEMENTBODY   = ifStatementBody$getCode()))
            }
        },
        
        makeDeclareSizeField = function(firstSize, nDim) {
            eval(substitute(switch(as.character(nDim),
                                   `0` = quote(FIRSTSIZE),
                                   `1` = quote(c(FIRSTSIZE, d)),
                                   `2` = quote(c(FIRSTSIZE, d, d)),
                                   stop()),
                            list(FIRSTSIZE = firstSize)))
        },
        
        makeIndexedVariable = function(varName, nDim) {
            eval(substitute(switch(as.character(nDim),
                                   `0` = quote(VARNAME[i]),
                                   `1` = quote(VARNAME[i, 1:d]),
                                   `2` = quote(VARNAME[i, 1:d, 1:d]),
                                   stop()),
                            list(VARNAME = varName)))
        }
    )
)

dependentClass <- setRefClass(
    Class = 'dependentClass',
    fields = list(
        distribution = 				'ANY',   ## the name of the (dependent) sampling distribution, e.g. 'dnorm'
        param = 					'ANY', 	 ## the name of the sampling distribution parameter in which targetNode must appear
        contributionExprs = 		'ANY', 	 ## a (named) list of expressions, giving the (additive) contribution to any parameters of the posterior. names correspond to variables in the posterior expressions
        contributionNames = 		'ANY', 	 ## names of the contributions to the parameters of the posterior distribution.  same as names(posteriorExprs)
        neededParamsForPosterior = 	'ANY'  	 ## names of all parameters appearing in the posteriorExprs
    ),
    methods = list(
        initialize = function(depInfoList, depDistName) {
        	contributionExprs <<- list()
            distribution <<- depDistName
            param <<- depInfoList$param
            initialize_contributionExprs(depInfoList)
            initialize_neededParamsForPosterior()
        },
        initialize_contributionExprs = function(depInfoList) {
            depInfoList['param'] <- NULL
            depInfoList <- lapply(depInfoList, function(di) parse(text=di)[[1]])
            contributionExprs <<- depInfoList
            contributionNames <<- names(depInfoList)
        },
        initialize_neededParamsForPosterior = function() {
            posteriorVars <- unique(unlist(lapply(contributionExprs, all.vars)))
            neededParamsForPosterior <<- unique(c(param, posteriorVars[!(posteriorVars %in% c('value', 'coeff', 'offset'))]))
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
            dDistribution <<- as.character(posteriorExpr[[1]])
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


cc_makeSamplerTypeName       <- function(distName)     return(paste0('conjugate_', distName))        ## 'dnorm' --> 'conjugate_dnorm'
cc_makeConjugateSamplerName  <- function(samplerType)  return(paste0('sampler_', samplerType))       ## 'conjugate_dnorm' --> 'sampler_conjugate_dnorm'
cc_makeRDistributionName     <- function(distName)     return(paste0('r', substring(distName, 2)))   ## 'dnorm' --> 'rnorm'

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
#    if(!any(node == model$getNodeNames()))   stop(paste0('node not present in model: ', node))			#Faster to just check if call fails
	output <- NULL
    try(output <- model$getNodeInfo()[[node]]$getValueExpr(), silent = TRUE)
    if(is.null(output))
    	stop(paste0('node not present in model: ', node) )
    return(output)
}

## special name used to represent vectors / arrays defined in terms of other stoch/determ nodes
cc_structureExprName <- quote(`_structureExpr`)

## expands all deterministic nodes in expr, to create a single expression with only stochastic nodes
cc_expandDetermNodesInExpr <- function(expr, model) {
    if(is.numeric(expr))     return(expr)     # return numeric
  

    if(is.name(expr)    ||    (is.call(expr) && (expr[[1]] == '['))) {    # expr is a name, or an indexed name
        exprText <- deparse(expr)
        
        
        graphID = NULL
        try(graphID <- model$modelDef$nodeName2GraphIDs(exprText), silent = TRUE)
        if(is.numeric(graphID)){
        	thisType <- model$modelDef$maps$types[graphID]
        	if(any(thisType == 'stoch') || any(thisType == 'LHSinferred') )
        		return(expr)
        	if(any(thisType == 'determ') ){
            if(length(model$expandNodeNames(exprText)) != 1){
              newExpr <- cc_createStructureExpr_fromModel(expr, model)
              for(i in seq_along(newExpr)[-1])    newExpr[[i]] <- cc_expandDetermNodesInExpr(newExpr[[i]], model)
              return(newExpr)
              }      
#            }
#        	  if(is.vectorized(exprText)) {
#        	    newExpr <- cc_createStructureExpr(expr)
#        	    for(i in seq_along(newExpr)[-1])    newExpr[[i]] <- cc_expandDetermNodesInExpr(newExpr[[i]], model)
#        	    return(newExpr)
#        	  }      
        	  
            return(cc_expandDetermNodesInExpr(expr=cc_getNodeValueExpr(model,node=exprText), model))

        }
        else
				stop(paste0('something went wrong processing: ', deparse(expr)))
			}
        
#        if(any(exprText == model$getMaps('nodeNamesStoch')))        return(expr)      # return stochastic nodes
#        if(any(exprText == model$getMaps('nodeNamesLHSinferred')))  return(expr)      # return LHS nodes inferred from a multivariate stochastic distribution
#        if(any(exprText == model$getMaps('nodeNamesDeterm')))
#            return(cc_expandDetermNodesInExpr(expr=cc_getNodeValueExpr(model,node=exprText), model))   # precess and return the value expression for this deterministic node
#        if(any(exprText == model$getMaps('nodeNamesRHSonly')))      stop('something wrong with model; possible failure to specify constants = ..., for a RHS-only node')
        return(expr)   # rather than throw an error, return expr; for the case where expr is the name of an array memberData object
    }
    if(is.call(expr)) {
        for(i in seq_along(expr)[-1])    expr[[i]] <- cc_expandDetermNodesInExpr(expr[[i]], model)
        return(expr) }
    stop(paste0('something went wrong processing: ', deparse(expr)))
}

## creates an expression of the form [cc_structureExprName](element11, element12, etc...) to represent vectors / arrays defined in terms of other stoch/determ nodes
cc_createStructureExpr <- function(expr) {
    expandedNodeNamesVector <- nl_expandNodeIndexExpr(expr)
    expandedNodeExprList <- lapply(expandedNodeNamesVector, function(x) parse(text=x)[[1]])
    structureExpr <- c(cc_structureExprName, expandedNodeExprList)
    structureExprCall <- as.call(structureExpr)
    return(structureExprCall)
}

## Same as above, but uses model to expandNodeNames
cc_createStructureExpr_fromModel <- function(expr, model) {
  expandedNodeNamesVector <- model$expandNodeNames(deparse(expr))
  expandedNodeExprList <- lapply(expandedNodeNamesVector, function(x) parse(text=x)[[1]])
  structureExpr <- c(cc_structureExprName, expandedNodeExprList)
  structureExprCall <- as.call(structureExpr)
  return(structureExprCall)
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
        if(cc_vectorizedComponentCheck(targetNode, expr))   return(FALSE)
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

## if targetNode is vectorized: determines if any components of targetNode appear in expr
cc_vectorizedComponentCheck <- function(targetNode, expr) {
    if(!is.vectorized(targetNode))     return(FALSE)
    targetNodesExpanded <- nl_expandNodeIndex(targetNode)
    if(any(unlist(lapply(targetNodesExpanded, function(node) cc_nodeInExpr(node, expr))))) return(TRUE)  ## any components of *vectorized* targetNode appear in expr
    return(FALSE)
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
    
    if(!is.call(expr))   stop('expression is not a call object')
    
    ## process the expression contents of the parentheses
    if(expr[[1]] == '(')
        return(cc_checkLinearity(expr[[2]], targetNode))
    
    ## we'll just have to skip over asRow() and asCol(), so they don't mess up the linearity check
    if(expr[[1]] == 'asRow' || expr[[1]] == 'asCol') {
        return(cc_checkLinearity(expr[[2]], targetNode))
    }
    
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
        stop('problem with negation expression')
    }
    
    if(expr[[1]] == '+') {
        checkLinearityLHS <- cc_checkLinearity(expr[[2]], targetNode)
        checkLinearityRHS <- cc_checkLinearity(expr[[3]], targetNode)
        if(is.null(checkLinearityLHS) || is.null(checkLinearityRHS)) return(NULL)
        return(list(offset = cc_combineExprsAddition(checkLinearityLHS$offset, checkLinearityRHS$offset),
                    scale  = cc_combineExprsAddition(checkLinearityLHS$scale,  checkLinearityRHS$scale)))
    }
    
    if(expr[[1]] == '*' || expr[[1]] == '%*%') {
        isMatrixMult <- ifelse(expr[[1]] == '%*%', TRUE, FALSE)
        if(cc_nodeInExpr(targetNode, expr[[2]]) && cc_nodeInExpr(targetNode, expr[[3]])) return(NULL)
        checkLinearityLHS <- cc_checkLinearity(expr[[2]], targetNode)
        checkLinearityRHS <- cc_checkLinearity(expr[[3]], targetNode)
        if(is.null(checkLinearityLHS) || is.null(checkLinearityRHS)) return(NULL)
        if((checkLinearityLHS$scale != 0) && (checkLinearityRHS$scale != 0)) stop('incompatible scales in * operation')
        if(checkLinearityLHS$scale != 0) {
            return(list(offset = cc_combineExprsMultiplication(checkLinearityLHS$offset, checkLinearityRHS$offset, isMatrixMult),
                        scale  = cc_combineExprsMultiplication(checkLinearityLHS$scale,  checkLinearityRHS$offset, isMatrixMult))) }
        if(checkLinearityRHS$scale != 0) {
            return(list(offset = cc_combineExprsMultiplication(checkLinearityLHS$offset, checkLinearityRHS$offset, isMatrixMult),
                        scale  = cc_combineExprsMultiplication(checkLinearityLHS$offset, checkLinearityRHS$scale , isMatrixMult))) }
        stop('something went wrong')
    }
    
    if(expr[[1]] == '/') {
        if(cc_nodeInExpr(targetNode, expr[[3]])) return(NULL)
        checkLinearityLHS <- cc_checkLinearity(expr[[2]], targetNode)
        checkLinearityRHS <- cc_checkLinearity(expr[[3]], targetNode)
        if(is.null(checkLinearityLHS) || is.null(checkLinearityRHS)) return(NULL)
        if(checkLinearityLHS$scale == 0) stop('left hand side scale == 0')
        if(checkLinearityRHS$scale != 0) stop('right hand side scale is non-zero')
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

cc_combineExprsMultiplication <- function(expr1, expr2, isMatrixMult) {
    if((expr1 == 0) || (expr2 == 0)) return(0)
    if((expr1 == 1) && (expr2 == 1)) return(1)
    if((expr1 != 1) && (expr2 == 1)) return(expr1)
    if((expr1 == 1) && (expr2 != 1)) return(expr2)
    if(is.numeric(expr1) && is.numeric(expr2)) return(expr1 * expr2)
    if(!isMatrixMult) return(substitute(EXPR1  *  EXPR2, list(EXPR1=expr1, EXPR2=expr2)))
    if( isMatrixMult) return(substitute(EXPR1 %*% EXPR2, list(EXPR1=expr1, EXPR2=expr2)))
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
##createNamedObjectsFromList(conjugateSamplerDefinitions)
createNamedObjectsFromList(conjugateSamplerDefinitions, writeToFile = 'TEMP_conjugateSamplerDefinitions.R')











