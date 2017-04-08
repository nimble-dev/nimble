
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
             ## dcat      = list(param = 'prob', contribution_alpha = as.numeric((1:length(prob)) == value)'),
             ## dcat      = list(param = 'prob', contribution_alpha = {tmp = rep(0,length(prob)); tmp[value]=1; tmp}')),
         posterior = 'ddirch(alpha = prior_alpha + contribution_alpha)'),

    ## gamma
    list(prior = 'dgamma',
         link = 'multiplicative',
         dependents = list(
             dpois  = list(param = 'lambda', contribution_shape = 'value', contribution_rate = 'coeff'                           ),
             dnorm  = list(param = 'tau',    contribution_shape = '1/2',   contribution_rate = 'coeff/2 * (value-mean)^2'        ),
             dlnorm = list(param = 'taulog', contribution_shape = '1/2',   contribution_rate = 'coeff/2 * (log(value)-meanlog)^2'),
             dgamma = list(param = 'rate',   contribution_shape = 'shape', contribution_rate = 'coeff   * value'                 ),
             dinvgamma = list(param = 'scale',   contribution_shape = 'shape', contribution_rate = 'coeff / value'                 ),
             dexp   = list(param = 'rate',   contribution_shape = '1',     contribution_rate = 'coeff   * value'                 ),
             dweib   = list(param = 'lambda',   contribution_shape = '1',     contribution_rate = 'coeff   * value^shape' )),
             ## ddexp  = list(param = 'rate',   contribution_shape = '1',     contribution_rate = 'coeff   * abs(value-location)'   )
             ## dpar = list(...)    ## contribution_shape=1; contribution_rate=coeff*log(value/c) 'c is 2nd param of pareto'
         posterior = 'dgamma(shape = prior_shape + contribution_shape,
                             scale = 1 / (prior_rate + contribution_rate))'),

    ## invgamma
    list(prior = 'dinvgamma',
         link = 'multiplicative',
         dependents = list(
             dnorm  = list(param = 'var',    contribution_shape = '1/2',   contribution_scale = '(value-mean)^2 / (coeff * 2)'),
             dlnorm = list(param = 'varlog', contribution_shape = '1/2',   contribution_scale = '(log(value)-meanlog)^2 / (coeff*2)'),
             dgamma = list(param = 'scale',   contribution_shape = 'shape', contribution_scale = 'value / coeff'                 ),
             dinvgamma = list(param = 'rate',   contribution_shape = 'shape', contribution_scale = '1 / (coeff * value)'     ),
             dexp   = list(param = 'scale',   contribution_shape = '1',     contribution_scale = 'value / coeff'            )),
             ## add ddexp
         posterior = 'dinvgamma(shape = prior_shape + contribution_shape,
                             rate = 1 / (prior_scale + contribution_scale))'),
    
    ## normal
    list(prior = 'dnorm',
         link = 'linear',
         dependents = list(
             dnorm  = list(param = 'mean',    contribution_mean = 'coeff * (value-offset) * tau',      contribution_tau = 'coeff^2 * tau'),
             dlnorm = list(param = 'meanlog', contribution_mean = 'coeff * (log(value)-offset) * taulog', contribution_tau = 'coeff^2 * taulog')),
         posterior = 'dnorm(mean = (prior_mean*prior_tau + contribution_mean) / (prior_tau + contribution_tau),
                            sd   = (prior_tau + contribution_tau)^(-0.5))'),

    #####
    ## pareto
    ## these are idiosyncratic enough that we probably want to skip them
    ## list(prior = 'dpar',      ##### waiting for dpar() distribution
    ##      link = 'multiplicative',
    ##      dependents = list(
    ##          dunif = list(param = 'max', contribution_alpha = '1', contribution_c = 'value/coeff'),  # only works if 0 is min of the unif so this will be hard to do
    ## careful with next one - data seem to impose upper bound on posterior, so not clear this goes through
    ##          dpar  = list(param = 'c',   contribution_alpha = '-alpha'), contribution_c = '0'),
    ##      posterior = 'dpar(alpha = prior_alpha + contribution_alpha,
    ##                        c     = max(prior_c, contribution_c)'),
    #####

    ## multivariate-normal
    list(prior = 'dmnorm',
         link = 'linear',
         dependents = list(
           ##dmnorm = list(param = 'mean', contribution_mean = '(t(coeff) %*% prec %*% asCol(value-offset))[,1]', contribution_prec = 't(coeff) %*% prec %*% coeff')),
             dmnorm = list(param = 'mean', contribution_mean = '(calc_dmnormConjugacyContributions(coeff, prec, value-offset, 1))[,1]', contribution_prec = 'calc_dmnormConjugacyContributions(coeff, prec, value-offset, 2)')),
         ## original less efficient posterior definition:
         ## posterior = 'dmnorm_chol(mean       = (inverse(prior_prec + contribution_prec) %*% (prior_prec %*% asCol(prior_mean) + asCol(contribution_mean)))[,1],
         ##                          cholesky   = chol(prior_prec + contribution_prec),
         ##                          prec_param = 1)'),
         posterior = '{ R <- chol(prior_prec + contribution_prec)
                        A <- prior_prec %*% asCol(prior_mean) + asCol(contribution_mean)
                        mu <- backsolve(R, forwardsolve(t(R), A))[,1]
                        dmnorm_chol(mean = mu, cholesky = R, prec_param = 1) }'),


    ## wishart
    list(prior = 'dwish',
         link = 'linear',
         dependents = list(
             ## parentheses added to the contribution_R calculation:
             ## colVec * (rowVec * matrix)
             ## Chris is checking to see whether this makes a difference for Eigen
             ## -DT April 2016
             dmnorm = list(param = 'prec', contribution_R = 'asCol(value-mean) %*% (asRow(value-mean) %*% coeff)', contribution_df = '1')),
         posterior = 'dwish_chol(cholesky    = chol(prior_R + contribution_R),
                                 df          = prior_df + contribution_df,
                                 scale_param = 0)')

    ## inverse wishart
    list(prior = 'dinvwish',
         link = 'identity',
         dependents = list(
             dmnorm = list(param = 'var', contribution_S = 'asCol(value-mean) %*% (asRow(value-mean)', contribution_df = '1')),
         posterior = 'dwish_chol(cholesky    = chol(prior_S + contribution_S),
                                 df          = prior_df + contribution_df,
                                 scale_param = 1)')

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
        checkConjugacy = function(model, nodeIDs) {
            maps <- model$modelDef$maps
            nodeDeclIDs <- maps$graphID_2_declID[nodeIDs] ## declaration IDs of the nodeIDs
            declID2nodeIDs <- split(nodeIDs, nodeDeclIDs) ## nodeIDs grouped by declarationID
            ansList <- ansList2 <- list()
            for(i in seq_along(declID2nodeIDs)) {         ## For each group of nodeIDs from the same declarationID
                nodeIDsFromOneDecl <- declID2nodeIDs[[i]]
                firstNodeName <- maps$graphID_2_nodeName[nodeIDsFromOneDecl[1]]
                if(model$isTruncated(firstNodeName)) next   ## we say non-conjugate if the targetNode is truncated
                dist <- model$getDistribution(firstNodeName)

                conjugacyObj <- conjugacys[[dist]]
                if(is.null(conjugacyObj)) next
                
                # NO: insert logic here to check a single dependency and do next if can't be conjugate
                #model$getDependencies('mu[1]',self=F,stochOnly=T)
                #for( loop through deps )
                #    conjugacyObj$checkConjugacyOneDep(model, targetNode, depNode)
                #    if(not conj) next

                # now try to guess if finding paths will be more intensive than simply looking at target-dependent pairs, to avoid path finding when there is nested structure such as stick-breaking
                numPaths <- sapply(nodeIDsFromOneDecl, model$getDependencyPathCountOneNode)
                deps <- lapply(nodeIDsFromOneDecl, function(x) model$getDependencies(x, stochOnly = TRUE, self = FALSE))
                numDeps <- sapply(deps, length)

                if(max(numPaths) > sum(numDeps)) {
                    # max(numPaths) is reasonable guess at number of unique (by node) paths (though it overestimates number of unique (by declaration ID) paths; if we have to evaluate conjugacy for more paths than we would by simply looking at all pairs of target-dependent nodes, then just use node pairs
                    # note that it's not clear what criterion to use here since computational time is combination of time for finding all paths and then for evaluating conjugacy for unique (by declaration ID) paths, but the hope is to make a crude cut here that avoids path calculations when there would be a lot of them
                    ansList[[length(ansList)+1]] <- lapply(seq_along(nodeIDsFromOneDecl),
                        function(index) {
                            targetNode <- maps$graphID_2_nodeName[nodeIDsFromOneDecl[index]]
                            depEnds <- deps[[index]]
                            depTypes <- sapply(depEnds, function(x) conjugacyObj$checkConjugacyOneDep(model, targetNode, x))
                            if(!length(depTypes)) return(NULL)
                            if(!any(sapply(depTypes, is.null))) {
                                uniqueDepTypes <- unique(depTypes)
                                control <- lapply(uniqueDepTypes,
                                                  function(oneType) {
                                                      boolMatch <- depTypes == oneType
                                                      depEnds[boolMatch]
                                                  })
                                names(control) <- uniqueDepTypes
                                return(list(prior = conjugacyObj$prior, type = conjugacyObj$samplerType, target = targetNode, control = control))
                            } else return(NULL)
                        })
                    names(ansList[[length(ansList)]]) <- maps$graphID_2_nodeName[nodeIDsFromOneDecl]
                    
                } else {
                # determine conjugacy based on unique (by declaration ID) paths
                    depPathsByNode <- lapply(nodeIDsFromOneDecl, getDependencyPaths, maps = maps)  ## make list (by nodeID) of lists of paths through graph
                    depPathsByNode <- depPathsByNode[!unlist(lapply(depPathsByNode, function(x) is.null(x) || (length(x)==0)))]
                    depPathsByNodeLabels <- lapply(depPathsByNode, function(z)                     ## make character labels that match for same path through graph
                        unlist(lapply(z,
                                      function(x)
                                          paste(maps$graphID_2_declID[x[,1]], x[,2], collapse = '\r', sep='\r'))))
                    
                    depPathsByNodeUnlisted <- unlist(depPathsByNode, recursive = FALSE)
                    depPathsByNodeLabelsUnlisted <- unlist(depPathsByNodeLabels)
                    ##  uniquePaths <- unique(depPathsByNodeLabelsUnlisted)
                    uniquePathsUnlistedIndices <- split(seq_along(depPathsByNodeLabelsUnlisted), depPathsByNodeLabelsUnlisted)
                    
                    conjDepTypes <- character(length(uniquePathsUnlistedIndices))
                    for(j in seq_along(uniquePathsUnlistedIndices)) {
                        firstDepPath <- depPathsByNodeUnlisted[[ uniquePathsUnlistedIndices[[j]][1] ]]
                        targetNode <- maps$graphID_2_nodeName[firstDepPath[1,1]]
                        depNode <- maps$graphID_2_nodeName[firstDepPath[nrow(firstDepPath), 1]]
                        oneDepType <- conjugacyObj$checkConjugacyOneDep(model, targetNode, depNode)
                        conjDepTypes[j] <- if(is.null(oneDepType)) "" else oneDepType
                    }
                    
                    conjBool <- conjDepTypes != ""
                    names(conjDepTypes) <- names(conjBool) <- names(uniquePathsUnlistedIndices)
                    if(any(conjBool)) {
                        targetNodes <- unlist(lapply(depPathsByNode, function(x) if(is.null(x)) '_NO_DEPS_' else maps$graphID_2_nodeName[x[[1]][1,1]]))
                        ansList[[length(ansList)+1]] <- mapply(
                            function(targetNode, depPathsOneNode, depPathsLabelsOneNode) {
                                if(targetNode == '_NO_DEPS_') return(NULL) ## these should have already been weeded out
                                if(all(conjBool[depPathsLabelsOneNode])) {
                                    depTypes <- conjDepTypes[depPathsLabelsOneNode]
                                    depEnds <- maps$graphID_2_nodeName[ unlist(lapply(depPathsOneNode, function(x) x[nrow(x)])) ]
                                    uniqueDepTypes <- unique(depTypes)
                                    control <- lapply(uniqueDepTypes,
                                                      function(oneType) {
                                                          boolMatch <- depTypes == oneType
                                                          ## prevent multiple instances of same dependent node name
                                                          ## (via different graph dependency paths   -DT Oct 2016
                                                          ##depEnds[boolMatch]
                                                          unique(depEnds[boolMatch])
                                                      })
                                    names(control) <- uniqueDepTypes
                                    list(prior = conjugacyObj$prior, type = conjugacyObj$samplerType, target = targetNode, control = control)
                                }
                            },
                            targetNodes, depPathsByNode, depPathsByNodeLabels, USE.NAMES = TRUE, SIMPLIFY = FALSE)
                    }
                }
            }
            if(length(ansList) > 0) ansList <- do.call('c', ansList)
            ansList <- ansList[!sapply(ansList, is.null)]  # strips out any NULL values
            if(!length(ansList)) return(list()) else return(ansList)  # replaces empty named list with empty list
        },
        ## older version of checkConjugacy(), which checks nodes one at a time (much slower)
        ## deprecated
        ## -DT Nov. 2016
        ##checkConjugacy = function(model, nodes) {
        ##    ## checks conjugacy of multiple nodes at once.
        ##    ## the return object is a named list, containing the conjugacyResult lists
        ##    ## *only* for nodes which are conjugate
        ##    conjugacyResultsAll <- list()
        ##    declarationIDs <- model$getDeclID(nodes)
        ##    nodesSplitByDeclaration <- split(nodes, declarationIDs)
        ##    for(theseNodes in nodesSplitByDeclaration)
        ##        conjugacyResultsAll <- c(conjugacyResultsAll, checkConjugacy_singleDeclaration(model, theseNodes))
        ##    return(conjugacyResultsAll)
        ##},
        ##checkConjugacy_singleDeclaration = function(model, nodes) {
        ##    ##browser()   ## removed by DT, July 2015, not sure why this was here
        ##    if(model$isTruncated(nodes[1])) return(list())   ## we say non-conjugate if the targetNode is truncated
        ##    dist <- model$getDistribution(nodes[1])
        ##    if(!dist %in% names(conjugacys)) return(list())
        ##    conjugacyObj <- conjugacys[[dist]]
        ##    ## temporary -- but works fine!
        ##    retList <- list()
        ##    for(node in nodes) {
        ##        result <- conjugacyObj$checkConjugacy(model, node)
        ##        if(!is.null(result))   retList[[node]] <- result
        ##    }
        ##    return(retList)
        ##    ## END temporary -- but works fine!
        ##    ## next line: this would be the new, more efficient approach -- not yet implemented
        ##    ##conjugacyObj$checkConjugacyAll(model, nodes) -- not yet implemented
        ##},
        ## update May 2016: old (non-dynamic) system is no longer supported -DT
        ##generateConjugateSamplerDefinitions = function() {
        ##    conjugateSamplerDefinitions <- list()
        ##    for(conjugacyObj in conjugacys) {  # conjugacyObj is a conjugacyClass object
        ##        samplerName <- cc_makeConjugateSamplerName(conjugacyObj$samplerType)
        ##        conjugateSamplerDefinitions[[samplerName]] <- conjugacyObj$generateConjugateSamplerDef()    ## workhorse for creating conjugate sampler nimble functions
        ##    }
        ##    return(conjugateSamplerDefinitions)
        ##},
        generateDynamicConjugateSamplerDefinition = function(prior, dependentCounts) {
            ## conjugateSamplerDefinitions[[paste0('sampler_conjugate_', conjugacyResult$prior)]]  ## using original (non-dynamic) conjugate sampler functions
            conjugacys[[prior]]$generateConjugateSamplerDef(dynamic = TRUE, dependentCounts = dependentCounts)
        }
    )
)

setMethod(
    '[[',
    'conjugacyRelationshipsClass',
    function(x, i)   return(x$conjugacys[[i]])
)

conjugacyClass <- setRefClass(
    Class = 'conjugacyClass',
    fields = list(
        samplerType =         'ANY',   ## name of the sampler for this conjugacy class, e.g. 'conjugate_dnorm'
        prior =               'ANY',   ## name of the prior distribution, e.g. 'dnorm'
        link =                'ANY',   ## the link ('linear', 'multiplicative', or 'identity')
        dependents =          'ANY',   ## (named) list of dependentClass objects, each contains conjugacy information specific to a particular sampling distribution (name is sampling distribution name)
        dependentDistNames =  'ANY',   ## character vector of the names of all allowable dependent sampling distributions.  same as: names(dependents)
        posteriorObject =     'ANY',   ## an object of posteriorClass
        needsLinearityCheck = 'ANY',   ## logical specifying whether we need to do the linearity check; if the link is 'multiplicative' or 'linear'
        model =               'ANY'    ## ONLY EXISTS TO PREVENT A WARNING for '<<-', in the code for generating the conjugate sampler function
    ),
    methods = list(
        initialize = function(cr) {
            dependents <<- list()
            samplerType <<- cc_makeSamplerTypeName(cr$prior)
            prior <<- cr$prior
            link <<- cr$link
            initialize_addDependents(cr$dependents)
            needsLinearityCheck <<- link %in% c('multiplicative', 'linear')
            posteriorObject <<- posteriorClass(cr$posterior, prior)
            model <<- NA
            },

        initialize_addDependents = function(depList) {
            for(i in seq_along(depList)) {
                dependents[[i]] <<- dependentClass(depList[[i]], names(depList)[i])
            }
            names(dependents) <<- names(depList)
            dependentDistNames <<- names(dependents)
        },

        ## used by new checkConjugacy() system
        ## see checkConjugacy for more explanation of each step
        checkConjugacyOneDep = function(model, targetNode, depNode) {
            if(model$getDistribution(targetNode) != prior)     return(NULL)    # check prior distribution of targetNode
            if(model$isTruncated(depNode)) return(NULL)   # if depNode is truncated, then not conjugate
            depNodeDist <- model$getDistribution(depNode)
            if(!(depNodeDist %in% dependentDistNames))     return(NULL)    # check sampling distribution of depNode
            dependentObj <- dependents[[depNodeDist]]
            linearityCheckExpr <- model$getParamExpr(depNode, dependentObj$param)   # extracts the expression for 'param' from 'depNode'
            linearityCheckExpr <- cc_expandDetermNodesInExpr(model, linearityCheckExpr)
            if(!cc_nodeInExpr(targetNode, linearityCheckExpr))                return(NULL)
            if(cc_vectorizedComponentCheck(targetNode, linearityCheckExpr))   return(NULL)   # if targetNode is vectorized, make sure none of its components appear in expr
            linearityCheck <- cc_checkLinearity(linearityCheckExpr, targetNode)   # determines whether paramExpr is linear in targetNode
            if(!cc_linkCheck(linearityCheck, link))                           return(NULL)
            if(!cc_otherParamsCheck(model, depNode, targetNode))              return(NULL)   # ensure targetNode appears in only *one* depNode parameter expression
            return(paste0('dep_', depNodeDist))
        },
        ## workhorse for checking conjugacy
        ## is this still used?????? (since transition to the "new" checkConjugacy, which
        ## checks multiple nodes from BUGS declarations at a time)
        ## should investigate if still used.  I'm honestly not certain.
        ## -DT Nov. 2016
        checkConjugacy = function(model, targetNode) {
            if(model$getDistribution(targetNode) != prior)     return(NULL)    # check prior distribution of targetNode
            control <- list()

            depNodes <- model$getDependencies(targetNode, stochOnly = TRUE, self = FALSE)
            if(length(depNodes) == 0)  return(NULL)   # no dependent stochastic nodes: not conjugate, return NULL

            for(depNode in depNodes) {
                if(model$isTruncated(depNode)) return(NULL)   # if depNode is truncated, then not conjugate
                depNodeDist <- model$getDistribution(depNode)
                if(!(depNodeDist %in% dependentDistNames))     return(NULL)    # check sampling distribution of depNode
                dependentObj <- dependents[[depNodeDist]]
                linearityCheckExpr <- model$getNodeExpr(depNode, dependentObj$param)   # extracts the expression for 'param' from 'depNode'
                linearityCheckExpr <- cc_expandDetermNodesInExpr(model, linearityCheckExpr)
                ## next line prevents the following potential error:
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
                control <- addDependentNodeToControl(control, depNodeDist, depNode)
            }
            return(list(type=samplerType, target=targetNode, control=control))   # all dependent nodes passed the conjugacy check
        },

        addDependentNodeToControl = function(control, depNodeDist, depNode) {
            listName <- paste0('dep_', depNodeDist)
            control[[listName]] <- c(control[[listName]], depNode)
            control
        },

        ## workhorse for creating conjugate sampler nimble functions
        generateConjugateSamplerDef = function(dynamic = FALSE, dependentCounts) {
            if(!dynamic) stop('something went wrong, should never have dynamic = FALSE here')
            substitute(
                nimbleFunction(contains = sampler_BASE,
                               setup    = SETUPFUNCTION,
                               run      = RUNFUNTION,
                               methods  = list(getPosteriorLogDensity = GETPOSTERIORLOGDENSITYFUNCTION,
                                               reset                  = function() {}),
                               where    = getLoadingNamespace()
                ),
                list(SETUPFUNCTION                  = genSetupFunction(dependentCounts = dependentCounts),
                     RUNFUNTION                     = genRunFunction(dependentCounts = dependentCounts),
                     GETPOSTERIORLOGDENSITYFUNCTION = genGetPosteriorLogDensityFunction(dependentCounts = dependentCounts)
                )
            )
        },

        genSetupFunction = function(dependentCounts) {
            functionBody <- codeBlockClass()

            functionBody$addCode({
                calcNodes       <- model$getDependencies(target)
                calcNodesDeterm <- model$getDependencies(target, determOnly = TRUE)
            })

            ## if this conjugate sampler is for a multivariate node (i.e., nDim > 0), then we need to determine the size (d)
            if(getDimension(prior) > 0) {
                functionBody$addCode(d <- max(determineNodeIndexSizes(target)))
            }

            ## make nodeFunction lists (target and dependents)
            ## create lists dependent nodeFunctions and their lengths
            ## NEWNODEFXN changes
            for(iDepCount in seq_along(dependentCounts)) {
                distName <- names(dependentCounts)[iDepCount]
                ##depNodeValueNdim <- getDimension(distName)
                functionBody$addCode({
                    DEP_NODENAMES <- control$DEP_CONTROL_NAME
                    N_DEP <- length(control$DEP_CONTROL_NAME)
                }, list(DEP_NODENAMES    = as.name(paste0(  'dep_', distName, '_nodeNames')),
                        N_DEP            = as.name(paste0('N_dep_', distName)),
                        DEP_CONTROL_NAME = as.name(paste0(  'dep_', distName))))
                if(getDimension(distName) > 0) {
                    functionBody$addCode({
                        DEP_NODESIZES <- sapply(DEP_NODENAMES, function(node) max(determineNodeIndexSizes(node)), USE.NAMES = FALSE)
                        if(length(DEP_NODESIZES) == 1) DEP_NODESIZES <- c(DEP_NODESIZES, -1)    ## guarantee to be a vector, for indexing and size processing
                        DEP_NODESIZEMAX <- max(DEP_NODESIZES)
                    }, list(DEP_NODESIZES   = as.name(paste0('dep_', distName, '_nodeSizes')),
                            DEP_NODENAMES   = as.name(paste0('dep_', distName, '_nodeNames')),
                            DEP_NODESIZEMAX = as.name(paste0('dep_', distName, '_nodeSizeMax'))))
                }

                ## uncomment this block to move from declare() to setup outputs for some variables
                ## functionBody$addCode({
                ##     DEP_VALUES_VAR <- array(0, dim = DECLARE_SIZE)
                ## },
                ##                      list(DEP_VALUES_VAR         = as.name(paste0('dep_', distName, '_values')),
                ##                           DECLARE_SIZE           = makeDeclareSizeField(as.name(paste0('N_dep_',distName)), depNodeValueNdim)   ## won't run since adding another argument to makeDeclareSizeField
                ##                           ))
                ## neededParams <- dependents[[distName]]$neededParamsForPosterior
                ## for(param in neededParams) {
                ##     depNodeParamNdim <- getDimension(distName, param) 
                ##     ## NEWNODEFXN
                ##     functionBody$addCode(DEP_PARAM_VAR <- array(0, dim = DECLARE_SIZE),
                ##                          list(DEP_PARAM_VAR      = as.name(paste0('dep_', distName, '_', param)),              ## DECLARE() statement
                ##                               DECLARE_SIZE       = makeDeclareSizeField(as.name(paste0('N_dep_',distName)), depNodeParamNdim)))   ## won't run since adding another argument to makeDeclareSizeField
                ## }

                ## functionBody$addCode({
                ##     N_DEP <- length(control$DEP_CONTROL_NAME)
                ##     DEP_NODEFUNCTIONS <- nimbleFunctionList(NF_VIRTUAL)
                ##     for(iDep in 1:N_DEP)
                ##         DEP_NODEFUNCTIONS[[iDep]] <- model$nodes[[control$DEP_CONTROL_NAME[iDep]]]
                ## }, list(N_DEP             = as.name(paste0('N_dep_', distName)),
                ##         DEP_CONTROL_NAME  = as.name(paste0(  'dep_', distName)),
                ##         DEP_NODEFUNCTIONS = as.name(paste0(  'dep_', distName, '_nfs')),
                ##         NF_VIRTUAL        = as.name(paste0('node_stoch_', distName))))
            }

            functionDef <- quote(function(model, mvSaved, target, control) {})
            functionDef[[3]] <- functionBody$getCode()
            functionDef[[4]] <- NULL   ## removes the 'scrref' attribute
            return(functionDef)
        },

        genRunFunction = function(dependentCounts) {
            functionBody <- codeBlockClass()

            ## only if we're verifying conjugate posterior distributions: get initial targetValue, and modelLogProb -- getLogProb(model, calcNodes)

            if(getNimbleOption('verifyConjugatePosteriors')) {
                functionBody$addCode({
                    modelLogProb0 <- getLogProb(model, calcNodes)
                    origValue <- model[[target]]
                })
            }

            addPosteriorQuantitiesGenerationCode(functionBody = functionBody, dependentCounts = dependentCounts)    ## adds code to generate the quantities prior_xxx, and contribution_xxx

            ## generate new value, store, calculate, copy, etc...
            functionBody$addCode(posteriorObject$prePosteriorCodeBlock, quote = FALSE)
            functionBody$addCode({
                newValue <- RPOSTERIORCALL
                model[[target]] <<- newValue
                calculate(model, calcNodes)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
            }, list(RPOSTERIORCALL = posteriorObject$rCallExpr))
            ## only if we're verifying conjugate posterior distributions: figure out if conjugate posterior distribution is correct
            if(nimbleOptions()$verifyConjugatePosteriors) {
                functionBody$addCode({
                    modelLogProb1 <- getLogProb(model, calcNodes)
                    posteriorLogDensity0 <- DPOSTERIORCALL_ORIG
                    posteriorLogDensity1 <- DPOSTERIORCALL_NEW
                    posteriorVerification <- modelLogProb0 - posteriorLogDensity0 - modelLogProb1 + posteriorLogDensity1
                    if(abs(posteriorVerification) > 1e-8)     {
                        nimPrint('conjugate posterior density appears to be wrong, off by ', posteriorVerification)
                    }
                }, list(DPOSTERIORCALL_ORIG = eval(substitute(substitute(expr, list(VALUE=quote(origValue))), list(expr=posteriorObject$dCallExpr))),
                        DPOSTERIORCALL_NEW  = eval(substitute(substitute(expr, list(VALUE=quote(newValue))),  list(expr=posteriorObject$dCallExpr)))))
            }

            functionDef <- quote(function() {})
            functionDef[[3]] <- functionBody$getCode()
            functionDef[[4]] <- NULL   ## removes the 'scrref' attribute
            return(functionDef)
        },

        genGetPosteriorLogDensityFunction = function(dependentCounts) {
            functionBody <- codeBlockClass()

            addPosteriorQuantitiesGenerationCode(functionBody = functionBody, dependentCounts = dependentCounts)    ## adds code to generate the quantities prior_xxx, and contribution_xxx

            ## calculate and return the (log)density for the current value of target
            functionBody$addCode(posteriorObject$prePosteriorCodeBlock, quote = FALSE)
            functionBody$addCode({
                targetValue <- model[[target]]
                posteriorLogDensity <- DPOSTERIORCALL
                returnType(double())
                return(posteriorLogDensity)
            }, list(DPOSTERIORCALL = eval(substitute(substitute(expr, list(VALUE=quote(targetValue))), list(expr=posteriorObject$dCallExpr)))))

            functionDef <- quote(function() {})
            functionDef[[3]] <- functionBody$getCode()
            functionDef[[4]] <- NULL   ## removes the 'scrref' attribute
            return(functionDef)
        },

        addPosteriorQuantitiesGenerationCode = function(functionBody = functionBody, dependentCounts = dependentCounts) {

            ## get current value of prior parameters which appear in the posterior expression
            for(priorParam in posteriorObject$neededPriorParams) {
                functionBody$addCode(PRIOR_PARAM_VAR <- model$getParam(target[1], PARAM_NAME),
                                     list(PRIOR_PARAM_VAR = as.name(paste0('prior_', priorParam)),
                                          PARAM_NAME      =                          priorParam))
            }

            for(iDepCount in seq_along(dependentCounts)) {
                distName <- names(dependentCounts)[iDepCount]
                neededParams <- dependents[[distName]]$neededParamsForPosterior
                depNodeValueNdim <- getDimension(distName)

                forLoopBody <- codeBlockClass()

                ## DECLARE() statement for dependent node values
                ## NEWNODEFXN: no change needed in this clause (can be moved to setup)
                functionBody$addCode(declare(DEP_VALUES_VAR, double(DEP_VALUES_VAR_NDIM, DECLARE_SIZE)),                       ## DECLARE() statement
                                     list(DEP_VALUES_VAR         = as.name(paste0('dep_', distName, '_values')),               ## DECLARE() statement
                                          DEP_VALUES_VAR_NDIM    = 1 + depNodeValueNdim,                                       ## DECLARE() statement
                                          DECLARE_SIZE           = makeDeclareSizeField(as.name(paste0('N_dep_', distName)), as.name(paste0('dep_', distName, '_nodeSizeMax')), as.name(paste0('dep_', distName, '_nodeSizeMax')), depNodeValueNdim)))
                ##functionBody$addCode(thisNodeSize <- 0) ## annoyingly this is to avoid a windows compiler warning about possible uninitialized use of thisNodeSize
                ## get *value* of each dependent node
                ## NEWNODEFXN
                if(getDimension(distName) > 0) {
                    forLoopBody$addCode(thisNodeSize <- DEP_NODESIZES[iDep],
                                        list(DEP_NODESIZES = as.name(paste0('dep_', distName, '_nodeSizes'))))
                }
                forLoopBody$addCode(DEP_VALUES_VAR_INDEXED <- model$getParam(DEP_NODENAMES[iDep], 'value'),
                                    list(DEP_VALUES_VAR_INDEXED = makeIndexedVariable(as.name(paste0('dep_', distName, '_values')), depNodeValueNdim, indexExpr = quote(iDep), secondSize = quote(thisNodeSize), thirdSize = quote(thisNodeSize)),
                                         DEP_NODENAMES = as.name(paste0('dep_', distName,'_nodeNames'))))

                for(param in neededParams) {
                    depNodeParamNdim <- getDimension(distName, param)
                    ## DECLARE() statement for each dependent node *parameter* value
                    ## NEWNODEFXN - no change needed here (can be moved to setup)
                    functionBody$addCode(declare(DEP_PARAM_VAR, double(DEP_PARAM_VAR_NDIM, DECLARE_SIZE)),                     ## DECLARE() statement
                                         list(DEP_PARAM_VAR      = as.name(paste0('dep_', distName, '_', param)),              ## DECLARE() statement
                                              DEP_PARAM_VAR_NDIM = 1 + depNodeParamNdim,                                       ## DECLARE() statement
                                              DECLARE_SIZE       = makeDeclareSizeField(as.name(paste0('N_dep_', distName)), as.name(paste0('dep_', distName, '_nodeSizeMax')), as.name(paste0('dep_', distName, '_nodeSizeMax')), depNodeParamNdim)))

                    ## get *parameter values* for each dependent node
                    ## NEWNODEFXN
                    forLoopBody$addCode(DEP_PARAM_VAR_INDEXED <- model$getParam(DEP_NODENAMES[iDep], PARAM_NAME),
                                        list(DEP_PARAM_VAR_INDEXED = makeIndexedVariable(as.name(paste0('dep_', distName, '_', param)), depNodeParamNdim, indexExpr = quote(iDep), secondSize = quote(thisNodeSize), thirdSize = quote(thisNodeSize)),
                                             DEP_NODENAMES = as.name(paste0('dep_', distName,'_nodeNames')),
                                             PARAM_NAME    = param))
                }

                functionBody$addCode(for(iDep in 1:N_DEP) FORLOOPBODY,
                                     list(N_DEP       = as.name(paste0('N_dep_', distName)),
                                          FORLOOPBODY = forLoopBody$getCode()))
            }

            ## if we need to determine 'coeff' and/or 'offset'
            if(needsLinearityCheck) {
                targetNdim <- getDimension(prior)
                targetCoeffNdim <- switch(as.character(targetNdim), `0`=0, `1`=2, `2`=2, stop())

                for(iDepCount in seq_along(dependentCounts)) {
                    distName <- names(dependentCounts)[iDepCount]
                    ## NEWNODEFXN - no change needed here
                    functionBody$addCode({                                                               ## DECLARE() statement
                        declare(DEP_OFFSET_VAR, double(DEP_OFFSET_VAR_NDIM, DECLARE_SIZE_OFFSET))        ## DECLARE() statement
                        declare(DEP_COEFF_VAR,  double(DEP_COEFF_VAR_NDIM,  DECLARE_SIZE_COEFF))         ## DECLARE() statement
                    }, list(DEP_OFFSET_VAR      = as.name(paste0('dep_', distName, '_offset')),          ## DECLARE() statement
                            DEP_COEFF_VAR       = as.name(paste0('dep_', distName, '_coeff')),           ## DECLARE() statement
                            DEP_OFFSET_VAR_NDIM = 1 + targetNdim,                                        ## DECLARE() statement
                            DEP_COEFF_VAR_NDIM  = 1 + targetCoeffNdim,                                   ## DECLARE() statement
                            DECLARE_SIZE_OFFSET = makeDeclareSizeField(as.name(paste0('N_dep_', distName)), as.name(paste0('dep_', distName, '_nodeSizeMax')), as.name(paste0('dep_', distName, '_nodeSizeMax')), targetNdim),
                            DECLARE_SIZE_COEFF  = makeDeclareSizeField(as.name(paste0('N_dep_', distName)), as.name(paste0('dep_', distName, '_nodeSizeMax')), quote(d),                                          targetCoeffNdim)))
                }

                switch(as.character(targetNdim),
                       `0` = {
                           functionBody$addCode({
                               model[[target]] <<- 0
                               model$calculate(calcNodesDeterm)
                           })

                           for(iDepCount in seq_along(dependentCounts)) {
                               distName <- names(dependentCounts)[iDepCount]
                               ## NEWNODEFXN
                               functionBody$addCode(
                                   for(iDep in 1:N_DEP)
                                       DEP_OFFSET_VAR[iDep] <- model$getParam(DEP_NODENAMES[iDep], PARAM_NAME),
                                   list(N_DEP          = as.name(paste0('N_dep_', distName)),
                                        DEP_OFFSET_VAR = as.name(paste0('dep_', distName, '_offset')),
                                        DEP_NODENAMES  = as.name(paste0('dep_', distName,'_nodeNames')),
                                        PARAM_NAME     = dependents[[distName]]$param))
                           }

                           functionBody$addCode({
                               model[[target]] <<- 1
                               model$calculate(calcNodesDeterm)
                           })

                           for(iDepCount in seq_along(dependentCounts)) {
                               distName <- names(dependentCounts)[iDepCount]
                               ## NEWNODEFXN
                               functionBody$addCode(
                                   for(iDep in 1:N_DEP)
                                       DEP_COEFF_VAR[iDep] <- model$getParam(DEP_NODENAMES[iDep], PARAM_NAME) - DEP_OFFSET_VAR[iDep],
                                   list(N_DEP             = as.name(paste0('N_dep_', distName)),
                                        DEP_COEFF_VAR     = as.name(paste0('dep_', distName, '_coeff')),
                                        DEP_NODENAMES     = as.name(paste0('dep_', distName, '_nodeNames')),
                                        PARAM_NAME        = dependents[[distName]]$param,
                                        DEP_OFFSET_VAR    = as.name(paste0('dep_', distName, '_offset'))))
                           }
                       },
                       `1` = {
                           functionBody$addCode({
                               model[[target]] <<- rep(0, d)
                               model$calculate(calcNodesDeterm)
                           })

                           for(iDepCount in seq_along(dependentCounts)) {
                               distName <- names(dependentCounts)[iDepCount]
                               ## NEWNODEFXN - forgot to copy and comment old code
                               functionBody$addCode({
                                   for(iDep in 1:N_DEP) {
                                       thisNodeSize <- DEP_NODESIZES[iDep]
                                       DEP_OFFSET_VAR[iDep, 1:thisNodeSize] <- model$getParam(DEP_NODENAMES[iDep], PARAM_NAME)
                                   }
                               },
                                                    list(N_DEP          = as.name(paste0('N_dep_', distName)),
                                                         DEP_NODESIZES  = as.name(paste0('dep_', distName, '_nodeSizes')),
                                                         DEP_OFFSET_VAR = as.name(paste0('dep_', distName, '_offset')),
                                                         DEP_NODENAMES  = as.name(paste0('dep_', distName,'_nodeNames')),
                                                         PARAM_NAME     = dependents[[distName]]$param))
                           }

                           functionBody$addCode(unitVector <- rep(0, d))

                           forLoopBody <- codeBlockClass()
                           forLoopBody$addCode({
                               unitVector[sizeIndex] <- 1
                               model[[target]] <<- unitVector
                               unitVector[sizeIndex] <- 0
                               calculate(model, calcNodesDeterm)
                           })

                           for(iDepCount in seq_along(dependentCounts)) {
                               distName <- names(dependentCounts)[iDepCount]
                               ## NEWNODEFXN
                               forLoopBody$addCode({
                                   for(iDep in 1:N_DEP) {
                                       thisNodeSize <- DEP_NODESIZES[iDep]
                                       DEP_COEFF_VAR[iDep, 1:thisNodeSize, sizeIndex] <- model$getParam(DEP_NODENAMES[iDep], PARAM_NAME) - DEP_OFFSET_VAR[iDep, 1:thisNodeSize]
                                   }
                               },
                                                   list(N_DEP          = as.name(paste0('N_dep_', distName)),
                                                        DEP_NODESIZES  = as.name(paste0('dep_', distName, '_nodeSizes')),
                                                        DEP_COEFF_VAR  = as.name(paste0('dep_', distName, '_coeff')),
                                                        DEP_NODENAMES  = as.name(paste0('dep_', distName, '_nodeNames')),
                                                        PARAM_NAME     = dependents[[distName]]$param,
                                                        DEP_OFFSET_VAR = as.name(paste0('dep_', distName, '_offset'))))
                           }

                           functionBody$addCode(for(sizeIndex in 1:d) FORLOOPBODY,
                                                list(FORLOOPBODY = forLoopBody$getCode()))
                       },
                       `2` = {
                           functionBody$addCode({
                               ## I <- model[[target]] * 0
                               ## for(sizeIndex in 1:d)   { I[sizeIndex, sizeIndex] <- 1 }
                               ## model[[target]] <<- I   ## initially, propogate through X = I
                               I <- identityMatrix(d)
                               model[[target]] <<- I   ## initially, propogate through X = I
                               calculate(model, calcNodesDeterm)
                           })

                           for(iDepCount in seq_along(dependentCounts)) {
                               distName <- names(dependentCounts)[iDepCount]
                               functionBody$addCode({
                                   for(iDep in 1:N_DEP) {
                                       ##thisNodeSize <- DEP_NODESIZES[iDep]  ## not needed for targetDim=2 case ????
                                       DEP_OFFSET_VAR[iDep, 1:d, 1:d] <- model$getParam(DEP_NODENAMES[iDep], PARAM_NAME)  ## DEP_OFFSET_VAR = A+B(I) = A+B
                                   }
                               },
                                                    list(N_DEP          = as.name(paste0('N_dep_', distName)),
                                                         ##DEP_NODESIZES  = as.name(paste0('dep_', distName, '_nodeSizes')),  ## not needed for targetDim=2 case ????
                                                         DEP_OFFSET_VAR = as.name(paste0('dep_', distName, '_offset')),
                                                         DEP_NODENAMES  = as.name(paste0('dep_', distName,'_nodeNames')),
                                                         PARAM_NAME     = dependents[[distName]]$param))
                           }

                           functionBody$addCode({
                               model[[target]] <<- I * 2   ## now, propogate through X = 2I
                               calculate(model, calcNodesDeterm)
                           })

                           for(iDepCount in seq_along(dependentCounts)) {
                               distName <- names(dependentCounts)[iDepCount]
                               functionBody$addCode({
                                   for(iDep in 1:N_DEP) {
                                       ##thisNodeSize <- DEP_NODESIZES[iDep]  ## not needed for targetDim=2 case ????
                                       DEP_COEFF_VAR[iDep, 1:d, 1:d] <- model$getParam(DEP_NODENAMES[iDep], PARAM_NAME) - DEP_OFFSET_VAR[iDep, 1:d, 1:d]   ## DEP_COEFF_VAR = (A+2B)-(A+B) = B
                                   }
                               },
                                                    list(N_DEP          = as.name(paste0('N_dep_', distName)),
                                                         ##DEP_NODESIZES  = as.name(paste0('dep_', distName, '_nodeSizes')),  ## not needed for targetDim=2 case ????
                                                         DEP_COEFF_VAR  = as.name(paste0('dep_', distName, '_coeff')),
                                                         DEP_NODENAMES  = as.name(paste0('dep_', distName, '_nodeNames')),
                                                         PARAM_NAME     = dependents[[distName]]$param,
                                                         DEP_OFFSET_VAR = as.name(paste0('dep_', distName, '_offset'))))
                               functionBody$addCode({
                                   for(iDep in 1:N_DEP) {
                                       ##thisNodeSize <- DEP_NODESIZES[iDep]  ## not needed for targetDim=2 case ????
                                       DEP_OFFSET_VAR[iDep, 1:d, 1:d] <- DEP_OFFSET_VAR[iDep, 1:d, 1:d] - DEP_COEFF_VAR[iDep, 1:d, 1:d]   ## now, DEP_OFFSET_VAR = (A+B)-(B) = A
                                   }
                               },
                                                    list(N_DEP          = as.name(paste0('N_dep_', distName)),
                                                         ##DEP_NODESIZES  = as.name(paste0('dep_', distName, '_nodeSizes')),  ## not needed for targetDim=2 case ????
                                                         DEP_OFFSET_VAR = as.name(paste0('dep_', distName, '_offset')),
                                                         DEP_COEFF_VAR  = as.name(paste0('dep_', distName, '_coeff'))))
                           }

                       },
                       stop()
                       )

            } # end if(needsLinearityCheck)

            targetNdim <- getDimension(prior)
            targetCoeffNdim <- switch(as.character(targetNdim), `0`=0, `1`=2, `2`=2, stop())
            ## adding declarations for the contribution terms, to remove Windows compiler warnings, DT August 2015
            for(contributionName in posteriorObject$neededContributionNames) {
                contribNdim <- posteriorObject$neededContributionDims[[contributionName]]
                ## still also need declare() statements?? new addition August 2015, for multivarate case
                ##if(contribNdim > 0)
                ##    functionBody$addCode(declare(CONTRIB_NAME, double(DIM, SIZES)),
                ##                         list(CONTRIB_NAME = as.name(contributionName),
                ##                              DIM          = contribNdim,
                ##                              SIZES        = if(contribNdim==1) quote(d) else if(contribNdim==2) quote(c(d,d)) else stop()))
                functionBody$addCode(CONTRIB_NAME <- CONTRIB_INITIAL_DECLARATION,
                                     list(CONTRIB_NAME                = as.name(contributionName),
                                          CONTRIB_INITIAL_DECLARATION = switch(as.character(contribNdim),
                                              `0` = 0, `1` = quote(numeric(length = d)), `2` = quote(array(dim = c(d, d))), stop())))
            }

            for(iDepCount in seq_along(dependentCounts)) {
                distName <- names(dependentCounts)[iDepCount]

                if(!any(posteriorObject$neededContributionNames %in% dependents[[distName]]$contributionNames))     next
                depParamsAvailable <- dependents[[distName]]$neededParamsForPosterior

                ## don't allow ragged dependencies for 2D conjugate case.
                ## no such cases exist, and it causes a runtime size check compiler warning.
                ## nonRaggedSizeExpr used to replace quote(thisNodeSize) below.
                ## August 2016
                nonRaggedSizeExpr <- if(targetNdim < 2) quote(thisNodeSize) else quote(d)
                subList <- lapply(depParamsAvailable, function(param)
                    makeIndexedVariable(as.name(paste0('dep_', distName, '_', param)), getDimension(distName, param), indexExpr = quote(iDep), secondSize = nonRaggedSizeExpr, thirdSize = nonRaggedSizeExpr))
                names(subList) <- depParamsAvailable
                
                subList$value  <- makeIndexedVariable(as.name(paste0('dep_', distName, '_values')), getDimension(distName), indexExpr = quote(iDep), secondSize = nonRaggedSizeExpr, thirdSize = nonRaggedSizeExpr)
                subList$offset <- makeIndexedVariable(as.name(paste0('dep_', distName, '_offset')), targetNdim, indexExpr = quote(iDep), secondSize = nonRaggedSizeExpr, thirdSize = nonRaggedSizeExpr)
                subList$coeff  <- makeIndexedVariable(as.name(paste0('dep_', distName, '_coeff')),  targetCoeffNdim, indexExpr = quote(iDep), secondSize = nonRaggedSizeExpr, thirdSize = quote(d))
                
                forLoopBody <- codeBlockClass()

                if(getDimension(distName) > 0) {
                    if(targetNdim == 1) ## 1D
                        forLoopBody$addCode(thisNodeSize <- DEP_NODESIZES[iDep],
                                            list(DEP_NODESIZES = as.name(paste0('dep_', distName, '_nodeSizes'))))
                    else ## 2D
                        forLoopBody$addCode(if(DEP_NODESIZES[iDep] != d) print('runtime error with sizes of 2D conjugate sampler'),
                                            list(DEP_NODESIZES = as.name(paste0('dep_', distName, '_nodeSizes'))))
                }

                for(contributionName in posteriorObject$neededContributionNames) {
                    if(!(contributionName %in% dependents[[distName]]$contributionNames))     next
                    contributionExpr <- eval(substitute(substitute(EXPR, subList), list(EXPR=dependents[[distName]]$contributionExprs[[contributionName]])))
                    forLoopBody$addCode(CONTRIB_NAME <- CONTRIB_NAME + CONTRIB_EXPR,
                                        list(CONTRIB_NAME = as.name(contributionName), CONTRIB_EXPR = contributionExpr))
                }
                functionBody$addCode(for(iDep in 1:N_DEP) FORLOOPBODY,
                                     list(N_DEP       = as.name(paste0('N_dep_', distName)),
                                          FORLOOPBODY = forLoopBody$getCode()))
            }
            ##}
        }
    )
)

dependentClass <- setRefClass(
    Class = 'dependentClass',
    fields = list(
        distribution =             'ANY',   ## the name of the (dependent) sampling distribution, e.g. 'dnorm'
        param =                    'ANY',   ## the name of the sampling distribution parameter in which target must appear
        contributionExprs =        'ANY',   ## a (named) list of expressions, giving the (additive) contribution to any parameters of the posterior. names correspond to variables in the posterior expressions
        contributionNames =        'ANY',   ## names of the contributions to the parameters of the posterior distribution.  same as names(posteriorExprs)
        neededParamsForPosterior = 'ANY'    ## names of all parameters appearing in the posteriorExprs
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
        prePosteriorCodeBlock =   'ANY',   ## a quoted {...} code block containing  DSL code to execute before making posterior call, possibly empty
        posteriorExpr =	          'ANY',   ## the full, parsed, posterior distribution expression, e.g. dnorm(mean = prior_mean + ..., sd = ...)
        rDistribution =           'ANY',   ## the *R* name of the posterior distribution, e.g. 'rnorm'
        dDistribution =           'ANY',   ## the *R* name of the posterior density distribution, e.g. 'dnorm'
        argumentExprs =           'ANY',   ## (named) list of expressions for each argument to the posterior distribution. names are the posterior distribution argument names
        argumentNames =           'ANY',   ## character vector of the argument names to the posterior distribution.  same as: names(argumentExprs)
        rCallExpr =               'ANY',   ## the actual 'rnorm(1, ...)' call, which will be substituted into the conjugate sampler function
        dCallExpr =               'ANY',   ## the 'dnorm(value, ...)' call, which can be used to get values of the posterior density
        neededPriorParams =       'ANY',   ## the names of any prior parameters (e.g., 'mean') which appear in the posterior expression as 'prior_mean'
        neededContributionNames = 'ANY',   ## the names of contributions from dependent nodes, such as 'contribution_scale'
        neededContributionDims =  'ANY'    ## a named list of contribution dimensions (0, 1, or 2). List element names are, e.g., 'contribution_scale'
    ),
    methods = list(
        initialize = function(posteriorText, prior) {
            parsedTotalPosterior <- parse(text = posteriorText)[[1]]
            if(parsedTotalPosterior[[1]] != '{') parsedTotalPosterior <- substitute({POST}, list(POST = parsedTotalPosterior))
            prePosteriorCodeBlock <<- parsedTotalPosterior[-length(parsedTotalPosterior)]
            posteriorExpr <<- parsedTotalPosterior[[length(parsedTotalPosterior)]]
            rDistribution <<- cc_makeRDistributionName(as.character(posteriorExpr[[1]]))
            dDistribution <<- as.character(posteriorExpr[[1]])
            argumentExprs <<- as.list(posteriorExpr)[-1]
            argumentNames <<- names(argumentExprs)
            rCallExpr <<- as.call(c(as.name(rDistribution), 1, argumentExprs))
            dCallExpr <<- as.call(c(as.name(dDistribution), quote(VALUE), argumentExprs, log = 1))
            posteriorVars <- all.vars(parsedTotalPosterior)
            neededPriorParams <<- gsub('^prior_', '', posteriorVars[grepl('^prior_', posteriorVars)])
            neededContributionNames <<- posteriorVars[grepl('^contribution_', posteriorVars)]
            neededContributionDims <<- inferContributionTermDimensions(prior)
        },
        inferContributionTermDimensions = function(prior) {
            distToLookup <- if(dDistribution %in% distributions$namesVector) dDistribution else if(prior %in% distributions$namesVector) prior else stop('cannot locate prior or posterior distribution in conjugacy processing')
            targetNdim <- getDimension(distToLookup)
            ## if posterior distribution is univariate, assume all contributions are scalar
            if(targetNdim == 0) {
                theDims <- lapply(neededContributionNames, function(x) 0)
                names(theDims) <- neededContributionNames
                return(theDims)
            }
            ## if posterior distribution is multivariate, attempt to infer contribution dimensionality from the *name* of each contribution term
            theDims <- list()
            typeNamesAvailable <- getParamNames(distToLookup) 
            for(contribName in neededContributionNames) {
                contribNameBase <- gsub('^contribution_', '', contribName)
                if(contribNameBase %in% typeNamesAvailable) {
                    ## contribution base name matches a parameter name of the posterior
                    theDims[[contribName]] <- getDimension(distToLookup, contribNameBase)
                } else {
                    ## contribution base name doesn't match any parameter; can't easily infer the dimensionality
                    browser()
                    message('The NIMBLE conjugacy system is attempting to infer the dimensionality of the contribution term: ', contribName, '. However, since the posterior distribution is multivariate, and the contribution name doesn\'t match any parameter names of the posterior distribution, NIMBLE can\'t infer this one. This means the conjugacy system might need to be extended, to allow users to provide the dimensionality of contribution terms. Or perhaps something more clever. -DT August 2015')
                }
            }
            return(theDims)
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



## expands all deterministic nodes in expr, to create a single expression with only stochastic nodes
cc_expandDetermNodesInExpr <- function(model, expr) {
    if(is.numeric(expr)) return(expr)     # return numeric
    if(is.name(expr) || (is.call(expr) && (expr[[1]] == '['))) { # expr is a name, or an indexed name
        exprText <- deparse(expr)
        expandedNodeNamesRaw <- try(model$expandNodeNames(exprText), silent=TRUE)  # causes error when expr is the name of an array memberData object, which isn't a node name
        if(inherits(expandedNodeNamesRaw, 'try-error')) {
            ## this case should no longer ever occur, I believe, under the newNodeFxns system -DT May 2016
            stop('something wrong with Daniel\'s understanding of newNodeFxns system')
            ##### at this point, should only be a 'name', representing an array memberData object
            ##### if it's an indexed name, we'll throw an error.
            ###if(is.call(expr)) stop('something went wrong with Daniel\'s understanding of newNimbleModel')
            ###return(expr) # expr is the name of an array memberData object; rather than throw an error, return expr
        }
        ## if exprText is a node itself (and also part of a larger node), then we only want the expansion to be the exprText node:
        expandedNodeNames <- if(exprText %in% expandedNodeNamesRaw) exprText else expandedNodeNamesRaw
        if(length(expandedNodeNames) == 1 && (expandedNodeNames == exprText)) {
            ## expr is a single node in the model
            type <- model$getNodeType(exprText)
            if(length(type) > 1) {
                ## if exprText is a node itself (and also part of a larger node), then we only want the expansion to be the exprText node:
                if(exprText %in% expandedNodeNamesRaw) type <- type[which(exprText == expandedNodeNamesRaw)]
                else stop('something went wrong with Daniel\'s understanding of newNimbleModel')
            }
            if(type == 'stoch') return(expr)
            if(type == 'determ') {
                newExpr <- model$getValueExpr(exprText)
                return(cc_expandDetermNodesInExpr(model, newExpr))
            }
            if(type == 'RHSonly') return(expr)
            stop('something went wrong with Daniel\'s understanding of newNimbleModel')
        }
        ## next line no longer necessary? (DT, May 2015)
        ## if(is.name(expr)) return(expr) # rather than throw an error, return expr; for the case where expr is the name of an array memberData object
        newExpr <- cc_createStructureExpr(model, exprText)
        for(i in seq_along(newExpr)[-1])
            newExpr[[i]] <- cc_expandDetermNodesInExpr(model, newExpr[[i]])
        return(newExpr)
    }
    if(is.call(expr)) {
        for(i in seq_along(expr)[-1])
            expr[[i]] <- cc_expandDetermNodesInExpr(model, expr[[i]])
        return(expr)
    }
    stop(paste0('something went wrong processing: ', deparse(expr)))
}

## special name used to represent vectors / arrays defined in terms of other stoch/determ nodes
cc_structureExprName <- quote(structureExpr)

## creates an expression of the form [cc_structureExprName](element11, element12, etc...) to represent vectors / arrays defined in terms of other stoch/determ nodes,
cc_createStructureExpr <- function(model, exprText) {
  expandedNodeNamesVector <- model$expandNodeNames(exprText)
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
    paramsList <- as.list(model$getValueExpr(depNode)[-1])       # extracts the list of all parameters, for the distribution of depNode
    timesFound <- 0   ## for success, we'll find targetNode in only *one* parameter expression
    for(i in seq_along(paramsList)) {
        expr <- cc_expandDetermNodesInExpr(model, paramsList[[i]])
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

compareConjugacyLists <- function(C1, C2) {
    if(identical(C1, C2)) return(TRUE)
    if(!identical(names(C1), names(C2))) {cat('Names do not match\n'); return(FALSE)}
    for(i in seq_along(C1)) {
        if(!identical(C1[[i]]$type, C2[[i]]$type)) cat(paste0('type mismatch for i =',i))
        if(!identical(C1[[i]]$target, C2[[i]]$target)) cat(paste0('target mismatch for i =',i))
        if(!identical(C1[[i]]$target, C2[[i]]$target)) cat(paste0('target mismatch for i =',i))
        if(!identical(names(C1[[i]]$control), names(C2[[i]]$control))) cat(paste0('control names mismatch for i =',i,'. Skipping node comparison'))
        else {
            for(j in seq_along(C1[[i]]$control)) {
                if(!identical(sort(C1[[i]]$control[[j]]), sort(C2[[i]]$control[[j]]))) cat(paste0('target mismatch for i =',i, 'j =', j))
            }
        }
    }
}

createDynamicConjugateSamplerName <- function(prior, dependentCounts) {
    ##depString <- paste0(dependentCounts, names(dependentCounts), collapse='_')  ## including the numbers of dependents
    depString <- paste0(names(dependentCounts), collapse='_')                     ## without the numbers of each type of dependent node
    paste0('sampler_conjugate_', prior, '_', depString)
}

makeDeclareSizeField <- function(firstSize, secondSize, thirdSize, nDim) {
    eval(substitute(switch(as.character(nDim),
                           `0` = quote(FIRSTSIZE),
                           `1` = quote(c(FIRSTSIZE, SECONDSIZE)),
                           `2` = quote(c(FIRSTSIZE, SECONDSIZE, THIRDSIZE)),
                           stop()),
                    list(FIRSTSIZE  = firstSize,
                         SECONDSIZE = secondSize,
                         THIRDSIZE  = thirdSize)))
}

makeIndexedVariable <- function(varName, nDim, indexExpr, secondSize, thirdSize) {
    eval(substitute(switch(as.character(nDim),
                           `0` = quote(VARNAME[INDEXEXPR]),
                           `1` = quote(VARNAME[INDEXEXPR, 1:SECONDSIZE]),
                           `2` = quote(VARNAME[INDEXEXPR, 1:SECONDSIZE, 1:THIRDSIZE]),
                           stop()),
                    list(VARNAME    = varName,
                         INDEXEXPR  = indexExpr,
                         SECONDSIZE = secondSize,
                         THIRDSIZE  = thirdSize)))
}


##############################################################################################
##############################################################################################
## create object: conjugacyRelationshipsObject
## also, generate all conjugate sampler nimbleFunctions
## and a function to rebuild conjugate sampler functions
##############################################################################################
##############################################################################################


## this is still *necessary* (and exported):
conjugacyRelationshipsObject <- conjugacyRelationshipsClass(conjugacyRelationshipsInputList)


## update May 2016: old (non-dynamic) system is no longer supported -DT
## this is still created (and exported) because it's handy:
##conjugateSamplerDefinitions <- conjugacyRelationshipsObject$generateConjugateSamplerDefinitions()
# Rebuild conjugate sampler functions
##buildConjugateSamplerFunctions <- function(writeToFile = NULL) {
##    conjugacyRelationshipsObject <- conjugacyRelationshipsClass(conjugacyRelationshipsInputList)
##    conjugateSamplerDefinitions <- conjugacyRelationshipsObject$generateConjugateSamplerDefinitions()
##    createNamedObjectsFromList(conjugateSamplerDefinitions, writeToFile = writeToFile, envir = parent.frame())
##}
##buildConjugateSamplerFunctions(writeToFile = 'TEMP_conjugateSamplerDefinitions.R')


## here after is for handling of dynamic conjugate sampler function

dynamicConjugateSamplerDefinitionsEnv <- new.env()
dynamicConjugateSamplerFunctionsEnv <- new.env()

dynamicConjugateSamplerExists <- function(name) {
    return(name %in% ls(dynamicConjugateSamplerDefinitionsEnv))
}

dynamicConjugateSamplerAdd <- function(name, def) {
    dynamicConjugateSamplerDefinitionsEnv[[name]] <- def
    dynamicConjugateSamplerFunctionsEnv[[name]] <- eval(def)
}

dynamicConjugateSamplerGet <- function(name) {
    return(dynamicConjugateSamplerFunctionsEnv[[name]])
}

dynamicConjugateSamplerWrite <- function(file = 'TEMP_dynamicConjugateSamplerDefinitions.R') {
    ## environment to create functions in is throw-away, since they've all already been created in dynamicConjugateSamplerFunctionsEnv
    createNamedObjectsFromList(as.list(dynamicConjugateSamplerDefinitionsEnv), writeToFile = file, envir = new.env())
}
















