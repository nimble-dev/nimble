
conjugacyRelationshipsInputList <- list(

    ## beta
    list(prior = 'dbeta',
         link = 'identity',
         dependents = list(
             dbern   = list(param = 'prob', contribution_shape1 = 'value', contribution_shape2 = '1 - value'   ),
             dbin    = list(param = 'prob', contribution_shape1 = 'value', contribution_shape2 = 'size - value'),
             dnegbin = list(param = 'prob', contribution_shape1 = 'size',  contribution_shape2 = 'value'       ),
             ## using 'stickbreaking' not 'stick_breaking' as link is now appended to dist with a '_' in realized conjugacy system
             ## Note that 'stick_breaking' is still the user-facing **function** name.
             dcat    = list(param = 'prob', link = 'stickbreaking', contribution_shape1 = 'value == offset',
                                                                    contribution_shape2 = 'value > offset')),  ## offset is set to be where in the broken stick the target is
         posterior = 'dbeta(shape1 = prior_shape1 + contribution_shape1,
                            shape2 = prior_shape2 + contribution_shape2)'),

    ## dirichlet
    list(prior = 'ddirch',
         link = 'identity',
         dependents = list(
             dmulti = list(param = 'prob', contribution_alpha = 'value'),
             dcat   = list(param = 'prob', contribution_alpha = 'calc_dcatConjugacyContributions(k, value)')),
         posterior = 'ddirch(alpha = prior_alpha + contribution_alpha)'),

    ## flat
    list(prior = 'dflat',
         link = 'linear',
         dependents = list(
             dnorm  = list(param = 'mean',    contribution_mean = 'coeff * (value-offset) * tau',         contribution_tau = 'coeff^2 * tau'),
             dlnorm = list(param = 'meanlog', contribution_mean = 'coeff * (log(value)-offset) * taulog', contribution_tau = 'coeff^2 * taulog')),
         posterior = 'dnorm(mean = contribution_mean / contribution_tau,
                            sd   = contribution_tau^(-0.5))'),
                                        
    ## halfflat - first possible conjugacy; user can achieve this with gamma(1, 0) as we handle that
    ## list(prior = 'dhalfflat',
    ##      link = 'multiplicative',
    ##      dependents = list(
    ##          dpois     = list(param = 'lambda', contribution_shape = 'value', contribution_rate = 'coeff'                           ),
    ##          dnorm     = list(param = 'tau',    contribution_shape = '1/2',   contribution_rate = 'coeff/2 * (value-mean)^2'        ),
    ##          dlnorm    = list(param = 'taulog', contribution_shape = '1/2',   contribution_rate = 'coeff/2 * (log(value)-meanlog)^2'),
    ##          dgamma    = list(param = 'rate',   contribution_shape = 'shape', contribution_rate = 'coeff   * value'                 ),
    ##          dinvgamma = list(param = 'scale',  contribution_shape = 'shape', contribution_rate = 'coeff   / value'                 ),
    ##          dexp      = list(param = 'rate',   contribution_shape = '1',     contribution_rate = 'coeff   * value'                 ),
    ##          dweib     = list(param = 'lambda', contribution_shape = '1',     contribution_rate = 'coeff   * value^shape'           )),
    ##          ## ddexp  = list(param = 'rate',   contribution_shape = '1',     contribution_rate = 'coeff   * abs(value-location)'   )
    ##          ## dpar = list(...)    ## contribution_shape=1; contribution_rate=coeff*log(value/c) 'c is 2nd param of pareto'
    ##      posterior = 'dgamma(shape = 1 + contribution_shape,
    ##                          scale = 1 / contribution_rate)'),

    ## halfflat - second possible conjugacy; note that user will be able to do invgamma(-1,0) and
    ## get conjugacy once we fix issue #314
    ## list(prior = 'dinvgamma',
    ##      link = 'multiplicative',
    ##      dependents = list(
    ##          dnorm     = list(param = 'var',    contribution_shape = '1/2',   contribution_scale = '(value-mean)^2 / (coeff * 2)'      ),
    ##          dlnorm    = list(param = 'varlog', contribution_shape = '1/2',   contribution_scale = '(log(value)-meanlog)^2 / (coeff*2)'),
    ##          dgamma    = list(param = 'scale',  contribution_shape = 'shape', contribution_scale = 'value / coeff'                     ),
    ##          dinvgamma = list(param = 'rate',   contribution_shape = 'shape', contribution_scale = '1 / (coeff * value)'               ),
    ##          dexp      = list(param = 'scale',  contribution_shape = '1',     contribution_scale = 'value / coeff'                     )),
    ##          ## add ddexp
    ##      posterior = 'dinvgamma(shape = -1 + contribution_shape,
    ##                          rate = 1 / contribution_scale)'),

    ## halfflat - third possible conjugacy - we use this because it corresponds to the
    ## Gelman (2006) recommended uniform on sd scale prior for variance components
    ## and current NIMBLE conjugacy system only allows one possible form of conjugacy.
    ## Note that sd ~ U(0,Inf) equivalent to var ~ IG(-1/2, 0).
    ## Also note that if conj system could detect 'squared' dependency, then
    ## we could allow dnorm with param = 'var'.                                     
    list(prior = 'dhalfflat',
         link = 'multiplicative',
         dependents = list(
             dnorm  = list(param = 'sd',    contribution_shape = '1/2',   contribution_scale = '(value-mean)^2 / (coeff^2 * 2)'        ),
             dlnorm = list(param = 'sdlog', contribution_shape = '1/2',   contribution_scale = '(log(value)-meanlog)^2 / (coeff^2 * 2)')),
         posterior = 'dsqrtinvgamma(shape = -1/2 + contribution_shape,
                             rate = 1 / contribution_scale)'),
         
    ## gamma
    list(prior = 'dgamma',
         link = 'multiplicative',
         dependents = list(
             dpois     = list(param = 'lambda', contribution_shape = 'value', contribution_rate = 'coeff'                           ),
             dnorm     = list(param = 'tau',    contribution_shape = '1/2',   contribution_rate = 'coeff/2 * (value-mean)^2'        ),
             dlnorm    = list(param = 'taulog', contribution_shape = '1/2',   contribution_rate = 'coeff/2 * (log(value)-meanlog)^2'),
             dgamma    = list(param = 'rate',   contribution_shape = 'shape', contribution_rate = 'coeff   * value'                 ),
             dinvgamma = list(param = 'scale',  contribution_shape = 'shape', contribution_rate = 'coeff   / value'                 ),
             dexp      = list(param = 'rate',   contribution_shape = '1',     contribution_rate = 'coeff   * value'                 ),
             dweib     = list(param = 'lambda', contribution_shape = '1',     contribution_rate = 'coeff   * value^shape'           ),
             ddexp     = list(param = 'rate',   contribution_shape = '1',     contribution_rate = 'coeff   * abs(value-location)'   )),
             ## dpar = list(...)    ## contribution_shape=1; contribution_rate=coeff*log(value/c) 'c is 2nd param of pareto'
         posterior = 'dgamma(shape = prior_shape + contribution_shape,
                             scale = 1 / (prior_rate + contribution_rate))'),

    ## invgamma
    list(prior = 'dinvgamma',
         link = 'multiplicative',
         dependents = list(
             dnorm     = list(param = 'var',    contribution_shape = '1/2',   contribution_scale = '(value-mean)^2 / (coeff * 2)'      ),
             dlnorm    = list(param = 'varlog', contribution_shape = '1/2',   contribution_scale = '(log(value)-meanlog)^2 / (coeff*2)'),
             dgamma    = list(param = 'scale',  contribution_shape = 'shape', contribution_scale = 'value / coeff'                     ),
             dinvgamma = list(param = 'rate',   contribution_shape = 'shape', contribution_scale = '1 / (coeff * value)'               ),
             dexp      = list(param = 'scale',  contribution_shape = '1',     contribution_scale = 'value / coeff'                     ),
             ddexp     = list(param = 'scale',  contribution_shape = '1',     contribution_scale = 'abs(value-location) / coeff'       )),
         posterior = 'dinvgamma(shape = prior_shape + contribution_shape,
                             rate = 1 / (prior_scale + contribution_scale))'),
    
    ## normal
    list(prior = 'dnorm',
         link = 'linear',
         dependents = list(
             dnorm  = list(param = 'mean',    contribution_mean = 'coeff * (value-offset) * tau',         contribution_tau = 'coeff^2 * tau'   ),
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

    ## multivariate-normal; note that to avoid unnecessary computations, there is special processing of the contribution_mean and contribution_prec
    ## in the generation of the conjugate sampler run method, so the specification here should not be changed without looking at that processing.
    list(prior = 'dmnorm',
         link = 'linear',
         dependents = list(
           ##dmnorm = list(param = 'mean', contribution_mean = '(t(coeff) %*% prec %*% asCol(value-offset))[,1]', contribution_prec = 't(coeff) %*% prec %*% coeff')),
             dmnorm = list(param = 'mean', contribution_mean = '(calc_dmnormConjugacyContributions(coeff, prec, value-offset, 1, 0))[,1]', contribution_prec = 'calc_dmnormConjugacyContributions(coeff, prec, value-offset, 2, 0)')),
         ## LINK will be replaced with appropriate link via code processing
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
         ## changing to only use link='identity' case, since the link='linear' case was not correct.
         ## -DT March 2017
         ## link = 'linear',
         link = 'multiplicativeScalar',  # we only handle scalar 'coeff'; this naming is slightly awkward since for univar dists, link is of course scalar
         dependents = list(
             ## parentheses added to the contribution_R calculation:
             ## colVec * (rowVec * matrix)
             ## Chris is checking to see whether this makes a difference for Eigen
             ## -DT April 2016
             ## changing to only use link='identity' case, since the link='linear' case was not correct
             ## -DT March 2017
             ## dmnorm = list(param = 'prec', contribution_R = 'asCol(value-mean) %*% (asRow(value-mean) %*% coeff)', contribution_df = '1')),
             dmnorm = list(param = 'prec', contribution_R = 'coeff * asCol(value-mean) %*% asRow(value-mean)', contribution_df = '1')),
         posterior = 'dwish_chol(cholesky    = chol(prior_R + contribution_R),
                                 df          = prior_df + contribution_df,
                                 scale_param = 0)'),

    ## inverse wishart
    list(prior = 'dinvwish',
         link = 'multiplicativeScalar',  # we only handle scalar 'coeff'; this naming is slightly awkward since for univar dists, link is of course scalar
         dependents = list(
             dmnorm = list(param = 'cov', contribution_S = 'asCol(value-mean) %*% asRow(value-mean) / coeff', contribution_df = '1')),
         posterior = 'dinvwish_chol(cholesky    = chol(prior_S + contribution_S),
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
        addOneConjugacy = function(XXXXXXXXXXXXXXXX) ??????
        checkConjugacy = function(model, nodeIDs, restrictLink = NULL) {
            if(isTRUE(getNimbleOption("oldConjugacyChecking")))
                checkConjugacy_original(model, nodeIDs, restrictLink)
            else
                checkConjugacy_new(model, nodeIDs, restrictLink)
        },
        checkConjugacy_new = function(model, nodeIDs, restrictLink = NULL) {
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

                # now try to guess if finding paths will be more intensive than simply looking at target-dependent pairs, to avoid path finding when there is nested structure such as stickbreaking
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
                            depTypes <- sapply(depEnds, function(x) conjugacyObj$checkConjugacyOneDep(model, targetNode, x, restrictLink))
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
                    depPathsByNode <- lapply(nodeIDsFromOneDecl, model$getDependencyPaths)  ## make list (by nodeID) of lists of paths through graph
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
                        oneDepType <- conjugacyObj$checkConjugacyOneDep(model, targetNode, depNode, restrictLink)
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

        checkConjugacy_original = function(model, nodeIDs, restrictLink = NULL) {
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

                # now try to guess if finding paths will be more intensive than simply looking at target-dependent pairs, to avoid path finding when there is nested structure such as stickbreaking
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
                            depTypes <- sapply(depEnds, function(x) conjugacyObj$checkConjugacyOneDep(model, targetNode, x, restrictLink))
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
                        oneDepType <- conjugacyObj$checkConjugacyOneDep(model, targetNode, depNode, restrictLink)
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
        generateDynamicConjugateSamplerDefinition = function(prior, dependentCounts, doDependentScreen = FALSE) {
            ## conjugateSamplerDefinitions[[paste0('sampler_conjugate_', conjugacyResult$prior)]]  ## using original (non-dynamic) conjugate sampler functions
            conjugacys[[prior]]$generateConjugateSamplerDef(dynamic = TRUE, dependentCounts = dependentCounts,
                                                            doDependentScreen = doDependentScreen)
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
        link =                'ANY',   ## the link ('linear', 'additive', 'multiplicative', 'multiplicativeScalar', or 'identity')
        dependents =          'ANY',   ## (named) list of dependentClass objects, each contains conjugacy information specific to a particular sampling distribution (name is sampling distribution name)
        dependentDistNames =  'ANY',   ## character vector of the names of all allowable dependent sampling distributions.  same as: names(dependents)
        posteriorObject =     'ANY',   ## an object of posteriorClass
        ## needsLinearityCheck = 'ANY',   ## logical specifying whether we need to do the linearity check; if the link is 'multiplicative' or 'linear'
        ## needsStickbreakingCheck = 'ANY',   ## logical specifying whether we need to do the stickbreaking check
        model                  = 'ANY',   ## these fields ONLY EXIST TO PREVENT A WARNING for '<<-',
        DEP_VALUES_VAR_INDEXED = 'ANY',   ## in the code for generating the conjugate sampler functions
        DEP_PARAM_VAR_INDEXED  = 'ANY',   ##
        DEP_OFFSET_VAR         = 'ANY',   ##
        DEP_COEFF_VAR          = 'ANY',   ##
        CONTRIB_NAME           = 'ANY'    ##
    ),
    methods = list(
        initialize = function(cr) {
            dependents <<- list()
            samplerType <<- cc_makeSamplerTypeName(cr$prior)
            prior <<- cr$prior
            link <<- cr$link
            initialize_addDependents(cr$dependents)
            ## removed because link info in specific conjugacy can now override default link
            ## needsLinearityCheck <<- link %in% c('multiplicative', 'linear')
            ## needsStickbreakingCheck <<- link %in% c('stickbreaking')
            posteriorObject <<- posteriorClass(cr$posterior, prior)
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
        checkConjugacyOneDep = function(model, targetNode, depNode, restrictLink = NULL) {
            if(model$getDistribution(targetNode) != prior)     return(NULL)    # check prior distribution of targetNode
            if(model$isTruncated(depNode)) return(NULL)   # if depNode is truncated, then not conjugate
            depNodeDist <- model$getDistribution(depNode)
            if(!(depNodeDist %in% dependentDistNames))     return(NULL)    # check sampling distribution of depNode
            dependentObj <- dependents[[depNodeDist]]
            if(is.null(restrictLink)) {
                if(!is.null(dependentObj$link)) currentLink <- dependentObj$link else currentLink <- link # handle multiple link case introduced for beta stickbreaking
            } else currentLink = restrictLink
            if(currentLink != 'stickbreaking') {
                depNodeParamName <- dependentObj$param
                linearityCheckExprRaw <- model$getParamExpr(depNode, depNodeParamName)   # extracts the expression for 'param' from 'depNode'
                linearityCheckExpr <- cc_expandDetermNodesInExpr(model, linearityCheckExprRaw, targetNode = targetNode)
                if(!cc_nodeInExpr(targetNode, linearityCheckExpr))                return(NULL)
                if(cc_vectorizedComponentCheck(targetNode, linearityCheckExpr))   return(NULL)   # if targetNode is vectorized, make sure none of its components appear in expr
                linearityCheck <- cc_checkLinearity(linearityCheckExpr, targetNode)   # determines whether paramExpr is linear in targetNode
                realizedLink <- cc_linkCheck(linearityCheck, currentLink)
                if(is.null(realizedLink)) return(NULL)
                ## ensure targetNode appears in only *one* depNode parameter expression
                if(!cc_otherParamsCheck(model, depNode, targetNode, depNodeExprExpanded = linearityCheckExpr, depParamNodeName = depNodeParamName))   return(NULL)
            } else {
                stickbreakingCheckExpr <- model$getParamExpr(depNode, dependentObj$param)   # extracts the expression for 'param' from 'depNode'
                stickbreakingCheckExpr <- cc_expandDetermNodesInExpr(model, stickbreakingCheckExpr, targetNode = targetNode)
                if(is.null(cc_checkStickbreaking(stickbreakingCheckExpr, targetNode))) return(NULL)
                realizedLink <- 'stickbreaking'
            }
            return(paste0('dep_', depNodeDist, '_', realizedLink))
        },

        addDependentNodeToControl = function(control, depNodeDist, depNode) {
            listName <- paste0('dep_', depNodeDist)
            control[[listName]] <- c(control[[listName]], depNode)
            control
        },

        ## workhorse for creating conjugate sampler nimble functions
        generateConjugateSamplerDef = function(dynamic = FALSE, dependentCounts, doDependentScreen = FALSE) {
            if(!dynamic) stop('something went wrong, should never have dynamic = FALSE here')
            substitute(
                nimbleFunction(contains = sampler_BASE,
                               setup    = SETUPFUNCTION,
                               run      = RUNFUNCTION,
                               methods  = list(getPosteriorLogDensity = GETPOSTERIORLOGDENSITYFUNCTION,
                                               reset                  = function() {})
                ),
                list(SETUPFUNCTION                  = genSetupFunction(dependentCounts = dependentCounts, doDependentScreen = doDependentScreen),
                     RUNFUNCTION                    = genRunFunction(dependentCounts = dependentCounts, doDependentScreen = doDependentScreen),
                     GETPOSTERIORLOGDENSITYFUNCTION = genGetPosteriorLogDensityFunction(dependentCounts = dependentCounts, doDependentScreen = doDependentScreen)
                )
            )
        },

        genSetupFunction = function(dependentCounts, doDependentScreen = FALSE) {
            functionBody <- codeBlockClass()
            functionBody$addCode({
                calcNodes       <- model$getDependencies(target)
                calcNodesDeterm <- model$getDependencies(target, determOnly = TRUE)
            })

            ## if this conjugate sampler is for a multivariate node (i.e., nDim > 0), then we need to determine the size (d)
            if(getDimension(prior) > 0) {
                functionBody$addCode(d <- max(determineNodeIndexSizes(target)))
            }
            for(iDepCount in seq_along(dependentCounts)) {
                distLinkName <- names(dependentCounts)[iDepCount]
                tmp <- strsplit(names(dependentCounts)[iDepCount], "_")[[1]]
                distName <- tmp[[1]]
                currentLink <- tmp[[2]]
                functionBody$addCode({
                    DEP_NODENAMES <- control$DEP_CONTROL_NAME
                    N_DEP <- length(control$DEP_CONTROL_NAME)
                }, list(DEP_NODENAMES    = as.name(paste0(  'dep_', distLinkName, '_nodeNames')),
                        N_DEP            = as.name(paste0('N_dep_', distLinkName)),
                        DEP_CONTROL_NAME = as.name(paste0(  'dep_', distLinkName))))
                if(any(getDimension(distName, includeParams = TRUE) > 0)) {
                    if(getDimension(distName) > 0) {  ## usual case: dependent node is multivariate -> get sizes from dependent node names
                        functionBody$addCode({
                            DEP_NODESIZES <- sapply(DEP_NODENAMES, function(node) max(determineNodeIndexSizes(node)), USE.NAMES = FALSE)
                        }, list(DEP_NODESIZES   = as.name(paste0('dep_', distLinkName, '_nodeSizes')),
                                DEP_NODENAMES   = as.name(paste0('dep_', distLinkName, '_nodeNames'))))
                    } else {  ## unusual case: dependent is univariate -> get sizes from dependent node's param (e.g., ddirch-dcat conjugacy)
                        functionBody$addCode({
                            DEP_NODESIZES <- sapply(DEP_NODENAMES, function(node) max(nimDim(model$getParam(node, MULTIVARIATE_PARAM_NAME))), USE.NAMES = FALSE)
                        }, list(DEP_NODESIZES           = as.name(paste0('dep_', distLinkName, '_nodeSizes')),
                                DEP_NODENAMES           = as.name(paste0('dep_', distLinkName, '_nodeNames')),
                                MULTIVARIATE_PARAM_NAME = names(which.max(getDimension(distName, includeParams = TRUE)))))
                    }
                    functionBody$addCode({
                        if(length(DEP_NODESIZES) == 1) DEP_NODESIZES <- c(DEP_NODESIZES, -1)    ## guarantee to be a vector, for indexing and size processing
                        DEP_NODESIZEMAX <- max(DEP_NODESIZES)
                    }, list(DEP_NODESIZES   = as.name(paste0('dep_', distLinkName, '_nodeSizes')),
                            DEP_NODESIZEMAX = as.name(paste0('dep_', distLinkName, '_nodeSizeMax'))))
                }
                
                ## declare() statements are removed from run() code,
                ## and were replaced with setup output array() calls below.
                ## July 2017
                depNodeValueNdim <- getDimension(distName)
                functionBody$addCode(DEP_VALUES_VAR <- array(0, dim = DECLARE_SIZE),
                                     list(DEP_VALUES_VAR = as.name(paste0('dep_', distLinkName, '_values')),
                                          DECLARE_SIZE   = makeDeclareSizeField(as.name(paste0('N_dep_', distLinkName)), as.name(paste0('dep_', distLinkName, '_nodeSizeMax')), as.name(paste0('dep_', distLinkName, '_nodeSizeMax')), depNodeValueNdim)))
                neededParams <- dependents[[distName]]$neededParamsForPosterior
                for(param in neededParams) {
                    depNodeParamNdim <- getDimension(distName, param)
                    functionBody$addCode(DEP_PARAM_VAR <- array(0, dim = DECLARE_SIZE),
                                         list(DEP_PARAM_VAR = as.name(paste0('dep_', distLinkName, '_', param)),
                                              DECLARE_SIZE  = makeDeclareSizeField(as.name(paste0('N_dep_', distLinkName)), as.name(paste0('dep_', distLinkName, '_nodeSizeMax')), as.name(paste0('dep_', distLinkName, '_nodeSizeMax')), depNodeParamNdim)))
                }
            }
            
            ## more new array() setup outputs, instead of declare() statements, for offset and coeff variables
            ## July 2017
            targetNdim <- getDimension(prior)
            targetCoeffNdim <- switch(as.character(targetNdim), `0`=0, `1`=2, `2`=2, stop())
            if(targetCoeffNdim == 2 && link == 'multiplicativeScalar')   ## Handles wish/invwish. There are no cases where we allow non-scalar 'coeff'.
                targetCoeffNdim <- 0
            for(iDepCount in seq_along(dependentCounts)) {
                distLinkName <- names(dependentCounts)[iDepCount]
                tmp <- strsplit(names(dependentCounts)[iDepCount], "_")[[1]]
                distName <- tmp[[1]]
                currentLink <- tmp[[2]]
                if(currentLink %in% c('additive', 'multiplicative', 'multiplicativeScalar', 'linear') || (nimbleOptions()$allowDynamicIndexing && doDependentScreen)) {
                    ## the 2's here are *only* to prevent warnings about assigning into member variable names using
                    inputList <-  list(DEP_OFFSET_VAR2     = as.name(paste0('dep_', distLinkName, '_offset')),  ## local assignment '<-', so changed the names to "...2"
                                       DEP_COEFF_VAR2      = as.name(paste0('dep_', distLinkName, '_coeff')),   ## so it doesn't recognize the ref class field name
                                       DECLARE_SIZE_OFFSET = makeDeclareSizeField(as.name(paste0('N_dep_', distLinkName)), as.name(paste0('dep_', distLinkName, '_nodeSizeMax')), as.name(paste0('dep_', distLinkName, '_nodeSizeMax')), targetNdim),
                                       DECLARE_SIZE_COEFF  = makeDeclareSizeField(as.name(paste0('N_dep_', distLinkName)), as.name(paste0('dep_', distLinkName, '_nodeSizeMax')), quote(d),                                          targetCoeffNdim))
                    if(currentLink == 'additive') 
                        functionBody$addCode(
                            DEP_OFFSET_VAR2  <- array(0, dim = DECLARE_SIZE_OFFSET),
                            inputList) 
                    if(currentLink %in% c('multiplicative', 'multiplicativeScalar')) 
                        functionBody$addCode(
                            DEP_COEFF_VAR2  <- array(0, dim = DECLARE_SIZE_COEFF),
                            inputList) 
                    if(currentLink == 'linear' || (nimbleOptions()$allowDynamicIndexing && doDependentScreen))
                        functionBody$addCode({
                            DEP_OFFSET_VAR2  <- array(0, dim = DECLARE_SIZE_OFFSET)
                            DEP_COEFF_VAR2  <- array(0, dim = DECLARE_SIZE_COEFF)
                        }, inputList) 
                }
            }

            for(iDepCount in seq_along(dependentCounts)) {
                distLinkName <- names(dependentCounts)[iDepCount]
                tmp <- strsplit(names(dependentCounts)[iDepCount], "_")[[1]]
                distName <- tmp[[1]]
                currentLink <- tmp[[2]]
                if(!is.null(dependents[[distName]]$link)) currentLink <- dependents[[distName]]$link
                if(currentLink == 'stickbreaking') {       
                    functionBody$addCode({
                        DEP_OFFSET_VAR2 <- array(0, dim = DECLARE_SIZE_OFFSET)                   ## the 2's here are *only* to prevent warnings about
                    }, list(DEP_OFFSET_VAR2     = as.name(paste0('dep_', distLinkName, '_offset')),  ## local assignment '<-', so changed the names to "...2"
                            DECLARE_SIZE_OFFSET = makeDeclareSizeField(as.name(paste0('N_dep_', distLinkName)), as.name(paste0('dep_', distLinkName, '_nodeSizeMax')), as.name(paste0('dep_', distLinkName, '_nodeSizeMax')), 0))
                    )
                    functionBody$addCode({
                        stickbreakingCheckExpr <- model$getValueExpr(calcNodesDeterm)
                        stickbreakingCheckExpr <- cc_expandDetermNodesInExpr(model, stickbreakingCheckExpr, targetNode = target)
                    })
                    functionBody$addCode(DEP_OFFSET_VAR2 <- rep(cc_checkStickbreaking(stickbreakingCheckExpr, target)$offset, DEP_OFFSET_SIZE), 
                                         list(DEP_OFFSET_VAR2 = as.name(paste0('dep_', distLinkName, '_offset')),
                                              DEP_OFFSET_SIZE = as.name(paste0('N_dep_', distLinkName))))
                }
            }
            ## adding declarations for the contribution terms, to remove Windows compiler warnings, DT August 2015
            ## moved these numeric() and array() declarations for contribution terms to setup outputs, July 2017
            for(contributionName in posteriorObject$neededContributionNames) {
                contribNdim <- posteriorObject$neededContributionDims[[contributionName]]
                functionBody$addCode(CONTRIB_NAME2 <- CONTRIB_INITIAL_DECLARATION,                   ## the 2's here are *only* to prevent warnings about
                                     list(CONTRIB_NAME2               = as.name(contributionName),   ## local assignment '<-' versus '<<-'
                                          CONTRIB_INITIAL_DECLARATION = switch(as.character(contribNdim),
                                              `0` = 0, `1` = quote(rep(0, length = d)), `2` = quote(array(0, dim = c(d, d))), stop())))
            }
            
            functionDef <- quote(function(model, mvSaved, target, control) {})
            functionDef[[3]] <- functionBody$getCode()
            functionDef[[4]] <- NULL   ## removes the 'scrref' attribute
            return(functionDef)
        },

        genRunFunction = function(dependentCounts, doDependentScreen = FALSE) {
            functionBody <- codeBlockClass()

            ## only if we're verifying conjugate posterior distributions: get initial targetValue, and modelLogProb -- model$getLogProb(calcNodes)

            if(getNimbleOption('verifyConjugatePosteriors')) {
                functionBody$addCode({
                    modelLogProb0 <- model$getLogProb(calcNodes)
                    origTargetValue <- model[[target]]
                })
            }

            ## adds code to generate the quantities prior_xxx, and contribution_xxx
            addPosteriorQuantitiesGenerationCode(functionBody = functionBody, dependentCounts = dependentCounts, doDependentScreen = doDependentScreen)

            ## generate new value, store, calculate, copy, etc...
            functionBody$addCode(posteriorObject$prePosteriorCodeBlock, quote = FALSE)
            functionBody$addCode({
                newTargetValue <- RPOSTERIORCALL
                model[[target]] <<- newTargetValue
                model$calculate(calcNodes)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
            }, list(RPOSTERIORCALL = posteriorObject$rCallExpr))
            ## only if we're verifying conjugate posterior distributions: figure out if conjugate posterior distribution is correct
            if(nimbleOptions()$verifyConjugatePosteriors) {
                functionBody$addCode({
                    modelLogProb1 <- model$getLogProb(calcNodes)
                    posteriorLogDensity0 <- DPOSTERIORCALL_ORIG
                    posteriorLogDensity1 <- DPOSTERIORCALL_NEW
                    posteriorVerification <- modelLogProb0 - posteriorLogDensity0 - modelLogProb1 + posteriorLogDensity1
                    if(abs(posteriorVerification) > 1e-8)     {
                        nimPrint('conjugate posterior density appears to be wrong, off by ', posteriorVerification)
                    }
                }, list(DPOSTERIORCALL_ORIG = eval(substitute(substitute(expr, list(VALUE=quote(origTargetValue))), list(expr=posteriorObject$dCallExpr))),
                        DPOSTERIORCALL_NEW  = eval(substitute(substitute(expr, list(VALUE=quote(newTargetValue))),  list(expr=posteriorObject$dCallExpr)))))
            }

            functionDef <- quote(function() {})
            functionDef[[3]] <- functionBody$getCode()
            functionDef[[4]] <- NULL   ## removes the 'scrref' attribute
            return(functionDef)
        },

        genGetPosteriorLogDensityFunction = function(dependentCounts, doDependentScreen = FALSE) {
            functionBody <- codeBlockClass()
            
            functionBody$addCode(origTargetValue <- model[[target]])
            
            ## adds code to generate the quantities prior_xxx, and contribution_xxx
            addPosteriorQuantitiesGenerationCode(functionBody = functionBody, dependentCounts = dependentCounts, doDependentScreen = doDependentScreen)

            ## calculate and return the (log)density for the current value of target
            functionBody$addCode(posteriorObject$prePosteriorCodeBlock, quote = FALSE)
            if(link != 'identity') {
                functionBody$addCode({
                    model[[target]] <<- origTargetValue
                    model$calculate(calcNodesDeterm)})
            }
            functionBody$addCode({
                posteriorLogDensity <- DPOSTERIORCALL
                returnType(double())
                return(posteriorLogDensity)
            }, list(DPOSTERIORCALL = eval(substitute(substitute(expr, list(VALUE=quote(origTargetValue))), list(expr=posteriorObject$dCallExpr)))))

            functionDef <- quote(function() {})
            functionDef[[3]] <- functionBody$getCode()
            functionDef[[4]] <- NULL   ## removes the 'scrref' attribute
            return(functionDef)
        },

        addPosteriorQuantitiesGenerationCode = function(functionBody = functionBody, dependentCounts = dependentCounts, doDependentScreen = FALSE) {
            ## get current value of prior parameters which appear in the posterior expression
            for(priorParam in posteriorObject$neededPriorParams) {
                functionBody$addCode(PRIOR_PARAM_VAR <- model$getParam(target[1], PARAM_NAME),
                                     list(PRIOR_PARAM_VAR = as.name(paste0('prior_', priorParam)),
                                          PARAM_NAME      =                          priorParam))
            }

            for(iDepCount in seq_along(dependentCounts)) {
                distLinkName <- names(dependentCounts)[iDepCount]
                tmp <- strsplit(names(dependentCounts)[iDepCount], '_')[[1]]
                distName <- tmp[[1]]
                currentLink <- tmp[[2]]
                neededParams <- dependents[[distName]]$neededParamsForPosterior
                depNodeValueNdim <- getDimension(distName)

                forLoopBody <- codeBlockClass()

                ## get *value* of each dependent node
                if(any(getDimension(distName, includeParams = TRUE) > 0)) {
                    forLoopBody$addCode(thisNodeSize <- DEP_NODESIZES[iDep],
                                        list(DEP_NODESIZES = as.name(paste0('dep_', distLinkName, '_nodeSizes'))))
                }
                forLoopBody$addCode(DEP_VALUES_VAR_INDEXED <<- model$getParam(DEP_NODENAMES[iDep], 'value'),
                                    list(DEP_VALUES_VAR_INDEXED = makeIndexedVariable(as.name(paste0('dep_', distLinkName, '_values')), depNodeValueNdim, indexExpr = quote(iDep), secondSize = quote(thisNodeSize), thirdSize = quote(thisNodeSize)),
                                         DEP_NODENAMES = as.name(paste0('dep_', distLinkName,'_nodeNames'))))

                for(param in neededParams) {
                    depNodeParamNdim <- getDimension(distName, param)
                    ## get *parameter values* for each dependent node
                    forLoopBody$addCode(DEP_PARAM_VAR_INDEXED <<- model$getParam(DEP_NODENAMES[iDep], PARAM_NAME),
                                        list(DEP_PARAM_VAR_INDEXED = makeIndexedVariable(as.name(paste0('dep_', distLinkName, '_', param)), depNodeParamNdim, indexExpr = quote(iDep), secondSize = quote(thisNodeSize), thirdSize = quote(thisNodeSize)),
                                             DEP_NODENAMES = as.name(paste0('dep_', distLinkName,'_nodeNames')),
                                             PARAM_NAME    = param))
                }

                functionBody$addCode(for(iDep in 1:N_DEP) FORLOOPBODY,
                                     list(N_DEP       = as.name(paste0('N_dep_', distLinkName)),
                                          FORLOOPBODY = forLoopBody$getCode()))
            }

            targetNdim <- getDimension(prior)
            allCurrentLinks <- sapply(names(dependentCounts), function(x) strsplit(x, '_')[[1]][[2]])
            switch(as.character(targetNdim),
                   `0` = {
                       if(any(allCurrentLinks %in% c('additive', 'linear')) || (nimbleOptions()$allowDynamicIndexing && doDependentScreen)) {
                           functionBody$addCode({
                               model[[target]] <<- 0
                               model$calculate(calcNodesDeterm)
                           })
                           for(iDepCount in seq_along(dependentCounts)) {
                               distLinkName <- names(dependentCounts)[iDepCount]
                               tmp <- strsplit(names(dependentCounts)[iDepCount], '_')[[1]]
                               distName <- tmp[[1]]
                               currentLink <- tmp[[2]]
                               
                               if(currentLink  %in% c('additive', 'linear') || (nimbleOptions()$allowDynamicIndexing && doDependentScreen)) 
                                   functionBody$addCode(
                                       for(iDep in 1:N_DEP)
                                           DEP_OFFSET_VAR[iDep] <<- model$getParam(DEP_NODENAMES[iDep], PARAM_NAME),
                                       list(N_DEP          = as.name(paste0('N_dep_', distLinkName)),
                                            DEP_OFFSET_VAR = as.name(paste0('dep_', distLinkName, '_offset')),
                                            DEP_NODENAMES  = as.name(paste0('dep_', distLinkName,'_nodeNames')),
                                            PARAM_NAME     = dependents[[distName]]$param))
                           }
                       }

                       if(any(allCurrentLinks == 'multiplicativeScalar'))
                           stop("Found 'multiplicativeScalar' link for 0-dimensional case.")
                       if(any(allCurrentLinks %in% c('multiplicative', 'linear')) || (nimbleOptions()$allowDynamicIndexing && doDependentScreen)) {
                           functionBody$addCode({
                               model[[target]] <<- 1
                               model$calculate(calcNodesDeterm)
                           })
                           
                           for(iDepCount in seq_along(dependentCounts)) {
                               distLinkName <- names(dependentCounts)[iDepCount]
                               tmp <- strsplit(names(dependentCounts)[iDepCount], '_')[[1]]
                               distName <- tmp[[1]]
                               currentLink <- tmp[[2]]
                               if(currentLink %in% c('multiplicative', 'linear') || (nimbleOptions()$allowDynamicIndexing && doDependentScreen)) {
                                   inputList <- list(N_DEP             = as.name(paste0('N_dep_', distLinkName)),
                                                     DEP_COEFF_VAR     = as.name(paste0('dep_', distLinkName, '_coeff')),
                                                     DEP_NODENAMES     = as.name(paste0('dep_', distLinkName, '_nodeNames')),
                                                     PARAM_NAME        = dependents[[distName]]$param,
                                                     DEP_OFFSET_VAR    = as.name(paste0('dep_', distLinkName, '_offset')))
                                   if(currentLink == 'linear'  || (nimbleOptions()$allowDynamicIndexing && doDependentScreen)) {
                                       functionBody$addCode(
                                           for(iDep in 1:N_DEP)
                                                     DEP_COEFF_VAR[iDep] <<- model$getParam(DEP_NODENAMES[iDep], PARAM_NAME) - DEP_OFFSET_VAR[iDep],
                                           inputList)
                                   } else {
                                       functionBody$addCode(
                                           for(iDep in 1:N_DEP)
                                                     DEP_COEFF_VAR[iDep] <<- model$getParam(DEP_NODENAMES[iDep], PARAM_NAME),
                                           inputList)
                                   }
                               }
                           }
                       }
                   },
                   `1` = {
                       if(any(allCurrentLinks == 'multiplicativeScalar'))
                           stop("Found 'multiplicativeScalar' link for 1-dimensional case.")
                       if(any(allCurrentLinks %in% c('additive', 'linear')) || (nimbleOptions()$allowDynamicIndexing && doDependentScreen)) {
                           functionBody$addCode({
                               model[[target]] <<- rep(0, d)
                               model$calculate(calcNodesDeterm)
                           })
                           
                           for(iDepCount in seq_along(dependentCounts)) {
                               distLinkName <- names(dependentCounts)[iDepCount]
                               tmp <- strsplit(names(dependentCounts)[iDepCount], '_')[[1]]
                               distName <- tmp[[1]]
                               currentLink <- tmp[[2]]

                               if(currentLink  %in% c('additive', 'linear') || (nimbleOptions()$allowDynamicIndexing && doDependentScreen)) 
                                   functionBody$addCode({
                                       for(iDep in 1:N_DEP) {
                                           thisNodeSize <- DEP_NODESIZES[iDep]
                                           DEP_OFFSET_VAR[iDep, 1:thisNodeSize] <<- model$getParam(DEP_NODENAMES[iDep], PARAM_NAME)
                                       }
                                   },
                                                        list(N_DEP          = as.name(paste0('N_dep_', distLinkName)),
                                                             DEP_NODESIZES  = as.name(paste0('dep_', distLinkName, '_nodeSizes')),
                                                             DEP_OFFSET_VAR = as.name(paste0('dep_', distLinkName, '_offset')),
                                                             DEP_NODENAMES  = as.name(paste0('dep_', distLinkName,'_nodeNames')),
                                                             PARAM_NAME     = dependents[[distName]]$param))
                           }
                       }
                       if(any(allCurrentLinks %in% c('multiplicative', 'linear')) || (nimbleOptions()$allowDynamicIndexing && doDependentScreen)) { 
                           functionBody$addCode(unitVector <- rep(0, d))
                           
                           forLoopBody <- codeBlockClass()
                           forLoopBody$addCode({
                               unitVector[sizeIndex] <- 1
                               model[[target]] <<- unitVector
                               unitVector[sizeIndex] <- 0
                               model$calculate(calcNodesDeterm)
                           })
                           
                           for(iDepCount in seq_along(dependentCounts)) {
                               distLinkName <- names(dependentCounts)[iDepCount]
                               tmp <- strsplit(names(dependentCounts)[iDepCount], '_')[[1]]
                               distName <- tmp[[1]]
                               currentLink <- tmp[[2]]
                               if(currentLink %in% c('multiplicative', 'linear') || (nimbleOptions()$allowDynamicIndexing && doDependentScreen)) {
                                   inputList <- list(N_DEP          = as.name(paste0('N_dep_', distLinkName)),
                                                     DEP_NODESIZES  = as.name(paste0('dep_', distLinkName, '_nodeSizes')),
                                                     DEP_COEFF_VAR  = as.name(paste0('dep_', distLinkName, '_coeff')),
                                                     DEP_NODENAMES  = as.name(paste0('dep_', distLinkName, '_nodeNames')),
                                                     PARAM_NAME     = dependents[[distName]]$param,
                                                     DEP_OFFSET_VAR = as.name(paste0('dep_', distLinkName, '_offset')))
                                   if(currentLink == 'linear' || (nimbleOptions()$allowDynamicIndexing && doDependentScreen)) {
                                       forLoopBody$addCode(
                                           for(iDep in 1:N_DEP) {
                                               thisNodeSize <- DEP_NODESIZES[iDep]
                                               DEP_COEFF_VAR[iDep, 1:thisNodeSize, sizeIndex] <<- model$getParam(DEP_NODENAMES[iDep], PARAM_NAME) - DEP_OFFSET_VAR[iDep, 1:thisNodeSize]
                                           }, inputList)
                                   } else {
                                       forLoopBody$addCode(
                                              for(iDep in 1:N_DEP) {
                                                  thisNodeSize <- DEP_NODESIZES[iDep]
                                                  DEP_COEFF_VAR[iDep, 1:thisNodeSize, sizeIndex] <<- model$getParam(DEP_NODENAMES[iDep], PARAM_NAME)
                                              }, inputList)
                                   }
                               }
                           }
                           functionBody$addCode(for(sizeIndex in 1:d) FORLOOPBODY,
                                                list(FORLOOPBODY = forLoopBody$getCode()))
                       }
                   },
                   `2` = {
                       if(!all(allCurrentLinks %in% c('identity', 'multiplicativeScalar')))
                           stop("Found non-multiplicative link for 2-d variable.")

                       if(any(allCurrentLinks == 'multiplicativeScalar') || (nimbleOptions()$allowDynamicIndexing && doDependentScreen)) {
                           functionBody$addCode({
                               model[[target]] <<- diag(d)
                               model$calculate(calcNodesDeterm)
                           })
                           
                           ## Use _COEFF_VAR to store value; we need this for stoch indexing case
                           ## where we determine that coeff = 0 because the potential dependency is not a dependency
                           ## given current index values.
                           for(iDepCount in seq_along(dependentCounts)) {
                               distLinkName <- names(dependentCounts)[iDepCount]
                               tmp <- strsplit(names(dependentCounts)[iDepCount], '_')[[1]]
                               distName <- tmp[[1]]
                               currentLink <- tmp[[2]]
                               if(currentLink == 'multiplicativeScalar' || (nimbleOptions()$allowDynamicIndexing && doDependentScreen)) 
                                   functionBody$addCode({
                                       for(iDep in 1:N_DEP) {
                                           DEP_COEFF_VAR[iDep] <<- model$getParam(DEP_NODENAMES[iDep], PARAM_NAME)[1, 1]   ## DEP_COEFF_VAR = (A+2B)-(A+B) = B
                                       }
                                   },
                                                        list(N_DEP          = as.name(paste0('N_dep_', distLinkName)),
                                                             DEP_COEFF_VAR  = as.name(paste0('dep_', distLinkName, '_coeff')),
                                                             DEP_NODENAMES  = as.name(paste0('dep_', distLinkName, '_nodeNames')),
                                                             PARAM_NAME     = dependents[[distName]]$param))
                           }
                           
                           functionBody$addCode({
                               ## Can't use zeros matrix as Cholesky fails; this is solely to determine if
                               ## potential dependency is not a dependency of target (due to stochastic indexing).
                               model[[target]] <<- 2*diag(d)  
                               model$calculate(calcNodesDeterm)
                           })
                           
                           for(iDepCount in seq_along(dependentCounts)) {
                               distLinkName <- names(dependentCounts)[iDepCount]
                               tmp <- strsplit(names(dependentCounts)[iDepCount], '_')[[1]]
                               distName <- tmp[[1]]
                               currentLink <- tmp[[2]]
                               if(currentLink == 'multiplicativeScalar'|| (nimbleOptions()$allowDynamicIndexing && doDependentScreen))
                                   functionBody$addCode({
                                       for(iDep in 1:N_DEP) {
                                           DEP_COEFF_VAR[iDep] <<- model$getParam(DEP_NODENAMES[iDep], PARAM_NAME)[1, 1] - DEP_COEFF_VAR[iDep]
                                       }
                                   },
                                                        list(N_DEP          = as.name(paste0('N_dep_', distLinkName)),
                                                             DEP_COEFF_VAR  = as.name(paste0('dep_', distLinkName, '_coeff')),
                                                             DEP_NODENAMES  = as.name(paste0('dep_', distLinkName, '_nodeNames')),
                                                             PARAM_NAME     = dependents[[distName]]$param))
                           }
                       }
                   },
                   stop()
                   )

            targetNdim <- getDimension(prior)
            ## contribution terms have been moved to be setup function outputs,
            ## but we still need to *zero these variables out* before adding contribution terms into them
            for(contributionName in posteriorObject$neededContributionNames) {
                contribNdim <- posteriorObject$neededContributionDims[[contributionName]]
                functionBody$addCode(CONTRIB_NAME <<- CONTRIB_ZERO_OUT,
                                     list(CONTRIB_NAME     = as.name(contributionName),
                                          CONTRIB_ZERO_OUT = switch(as.character(contribNdim),
                                              `0` = 0, `1` = quote(rep(0, length = d)), `2` = quote(array(0, dim = c(d, d))), stop())))
            }

            for(iDepCount in seq_along(dependentCounts)) {
                distLinkName <- names(dependentCounts)[iDepCount]
                tmp <- strsplit(names(dependentCounts)[iDepCount], '_')[[1]]
                distName <- tmp[[1]]
                currentLink <- tmp[[2]]
                targetCoeffNdim <- switch(as.character(targetNdim), `0`=0, `1`=2, `2`=2, stop())
                if(targetCoeffNdim == 2 && link == 'multiplicativeScalar')   ## There are no cases where we allow non-scalar 'coeff'.
                    targetCoeffNdim <- 0

                if(!any(posteriorObject$neededContributionNames %in% dependents[[distName]]$contributionNames))     next
                depParamsAvailable <- dependents[[distName]]$neededParamsForPosterior

                ## don't allow ragged dependencies for 2D conjugate case.
                ## no such cases exist, and it causes a runtime size check compiler warning.
                ## nonRaggedSizeExpr used to replace quote(thisNodeSize) below.
                ## August 2016
                nonRaggedSizeExpr <- if(targetNdim < 2) quote(thisNodeSize) else quote(d)
                subList <- lapply(depParamsAvailable, function(param)
                    makeIndexedVariable(as.name(paste0('dep_', distLinkName, '_', param)), getDimension(distName, param), indexExpr = quote(iDep), secondSize = nonRaggedSizeExpr, thirdSize = nonRaggedSizeExpr))
                names(subList) <- depParamsAvailable
                
                subList$value  <- makeIndexedVariable(as.name(paste0('dep_', distLinkName, '_values')), getDimension(distName), indexExpr = quote(iDep), secondSize = nonRaggedSizeExpr, thirdSize = nonRaggedSizeExpr)
                if(currentLink %in% c('additive', 'linear', 'stickbreaking'))
                    subList$offset <- makeIndexedVariable(as.name(paste0('dep_', distLinkName, '_offset')), targetNdim, indexExpr = quote(iDep), secondSize = nonRaggedSizeExpr, thirdSize = nonRaggedSizeExpr)
                if(currentLink %in% c('multiplicative', 'multiplicativeScalar', 'linear') || (nimbleOptions()$allowDynamicIndexing && doDependentScreen))
                    subList$coeff  <- makeIndexedVariable(as.name(paste0('dep_', distLinkName, '_coeff')),  targetCoeffNdim, indexExpr = quote(iDep), secondSize = nonRaggedSizeExpr, thirdSize = quote(d))
                
                forLoopBody <- codeBlockClass()
                
                if(any(getDimension(distName, includeParams = TRUE) > 0)) {
                    if(targetNdim == 1) ## 1D
                        forLoopBody$addCode(thisNodeSize <- DEP_NODESIZES[iDep],
                                            list(DEP_NODESIZES = as.name(paste0('dep_', distLinkName, '_nodeSizes'))))
                    if(targetNdim == 2) ## 2D  ## formerly this was 'else', but for 'dcat' we have targetNdim=0 while max(getDimension(distName, includeParams = TRUE)) is 1 so need explicit check for 2D
                        forLoopBody$addCode(if(DEP_NODESIZES[iDep] != d) print('runtime error with sizes of 2D conjugate sampler'),
                                            list(DEP_NODESIZES = as.name(paste0('dep_', distLinkName, '_nodeSizes'))))
                }
                for(contributionName in posteriorObject$neededContributionNames) {
                    if(!(contributionName %in% dependents[[distName]]$contributionNames))     next
                    contributionExpr <- dependents[[distName]]$contributionExprs[[contributionName]]
                    if(distName == 'dmnorm' && prior == 'dmnorm') {
                        ## need to deal with [,1] in contribution_mean
                        if(contributionName == 'contribution_mean') tmpExpr <- contributionExpr[[2]][[2]] else tmpExpr <- contributionExpr
                        if(nimbleOptions()$allowDynamicIndexing && doDependentScreen) {
                            ## Will always have 'coeff' but may not need to use it
                            if(currentLink %in% c('identity', 'additive')) 
                                tmpExpr[[6]] <- quote(0) else tmpExpr[[6]] <- quote(1)
                        } else {
                            if(currentLink %in%  c('identity', 'additive')) {
                                tmpExpr[[2]] <- quote(prec)  ## hack, since 'coeff' won't exist
                                tmpExpr[[6]] <- quote(0)
                            } else tmpExpr[[6]] <- quote(1)
                        }
                        tmpExpr[[4]] <- cc_stripExpr(tmpExpr[[4]], offset = currentLink %in% c('identity','multiplicative'), coeff = FALSE)  # strip 'offset'
                        if(contributionName == 'contribution_mean') contributionExpr[[2]][[2]] <- tmpExpr else contributionExpr <- tmpExpr
                    } else contributionExpr <- cc_stripExpr(contributionExpr, offset = currentLink %in% c('identity','multiplicative','multiplicativeScalar'),
                                                     coeff = currentLink %in% c('identity','additive'))
                    contributionExpr <- eval(substitute(substitute(EXPR, subList), list(EXPR=contributionExpr)))
                    if(nimbleOptions()$allowDynamicIndexing && doDependentScreen) { ## FIXME: would be nice to only have one if() here when we loop through multiple parameters
                        if(targetCoeffNdim == 0)
                            forLoopBody$addCode(if(COEFF_EXPR != 0) CONTRIB_NAME <<- CONTRIB_NAME + CONTRIB_EXPR,
                                                list(COEFF_EXPR = subList$coeff, CONTRIB_NAME = as.name(contributionName), CONTRIB_EXPR = contributionExpr))
                        else forLoopBody$addCode(if(min(COEFF_EXPR) != 0 | max(COEFF_EXPR) != 0) CONTRIB_NAME <<- CONTRIB_NAME + CONTRIB_EXPR,
                                                list(COEFF_EXPR = subList$coeff, CONTRIB_NAME = as.name(contributionName), CONTRIB_EXPR = contributionExpr))

                    } else forLoopBody$addCode(CONTRIB_NAME <<- CONTRIB_NAME + CONTRIB_EXPR,
                                        list(CONTRIB_NAME = as.name(contributionName), CONTRIB_EXPR = contributionExpr))
                }
                functionBody$addCode(for(iDep in 1:N_DEP) FORLOOPBODY,
                                     list(N_DEP       = as.name(paste0('N_dep_', distLinkName)),
                                          FORLOOPBODY = forLoopBody$getCode()))
            }
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
        neededParamsForPosterior = 'ANY',   ## names of all parameters appearing in the posteriorExprs
        link =                     'ANY',   ## optional specific link for a given conjugacy; currently only used for stickbreaking case
        userAdded =                'ANY',   ## boolean specifying whether this dependency is user-added via registerConjugacy
    ),
    methods = list(
        initialize = function(depInfoList, depDistName) {
        	contributionExprs <<- list()
                distribution <<- depDistName
                param <<- depInfoList$param
                link <<- depInfoList$link  ## will be NULL unless specific conjugacy overrides default link
                userAdded <<- XXXXXXXXXXXXXXX
                initialize_contributionExprs(depInfoList)
                initialize_neededParamsForPosterior()
        },
        initialize_contributionExprs = function(depInfoList) {
            depInfoList['param'] <- NULL
            depInfoList <- lapply(depInfoList, function(di) parse(text=di)[[1]])
            contributionExprs <<- depInfoList
            contributionExprs[['link']] <<- NULL   ## so link not included in neededParamsForPosterior
            contributionNames <<- names(depInfoList)
        },
        initialize_neededParamsForPosterior = function() {
            posteriorVars <- unique(unlist(lapply(contributionExprs, all.vars)))
            ## Formerly, we always included the target 'param' in 'neededParamsForPosterior'.
            ## Not certain why, it wasn't necessary in many cases, e.g., dbeta-dbinom conjugacy.  Excluding it now.
            ## -DT, August 2017
            ##neededParamsForPosterior <<- unique(c(param, posteriorVars[!(posteriorVars %in% c('value', 'coeff', 'offset'))]))
            neededParamsForPosterior <<- setdiff(posteriorVars, c('value', 'coeff', 'offset'))
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
                    stop(message('The NIMBLE conjugacy system is attempting to infer the dimensionality of the contribution term: ',
                                contribName, '. However, since the posterior distribution is multivariate, and the contribution name doesn\'t match any parameter names of the posterior distribution, NIMBLE can\'t infer this one. This means the conjugacy system might need to be extended, to allow users to provide the dimensionality of  contribution terms. Or perhaps something more clever. -DT August 2015'), call. = FALSE)
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
cc_expandDetermNodesInExpr <- function(model, expr, targetNode = NULL, skipExpansionsNode = NULL) {
    if(is.numeric(expr)) return(expr)     # return numeric
    if(is.name(expr) || (is.call(expr) && (expr[[1]] == '[') && is.name(expr[[2]]))) { # expr is a name, or an indexed name
        if(nimbleOptions()$allowDynamicIndexing) {
            ## this deals with having mu[k[1]] (which won't pass through expandNodeNames), replacing k[1] with the index from targetNode
            if(!is.name(expr)) {
                indexExprs <- expr[3:length(expr)]
                numericOrVectorIndices <- sapply(indexExprs,
                      function(x) is.numeric(x) || (length(x) == 3 && x[[1]] == ':'))
                if(!all(numericOrVectorIndices)) {
                    if(model$getVarNames(nodes = targetNode) == model$getVarNames(nodes = safeDeparse(expr, warn = TRUE))) {
                        ## expr var is same as target var, so plug in target indexes for
                        ## non-constant expr indexes to allow for possibility that
                        ## dynamic index will be that of target (in some cases based on allowed
                        ## values of dynamic index, this might actually not be possible, but
                        ## the only disdvantage should be failing to detect conjugacy where it is
                        ## present).
                        expr[which(!numericOrVectorIndices)+2] <- parse(text = targetNode)[[1]][3:length(expr)][!numericOrVectorIndices]
                        ## now check that target is not actually used in index expressions
                        newExpr <- as.call(c(cc_structureExprName, sapply(indexExprs[!numericOrVectorIndices], function(x) x)))
                        for(i in seq_along(newExpr)[-1])
                            newExpr[[i]] <- cc_expandDetermNodesInExpr(model, newExpr[[i]], targetNode, skipExpansionsNode)
                        if(cc_nodeInExpr(targetNode, newExpr))
                            return(as.call(c(cc_structureExprName, expr, sapply(indexExprs[!numericOrVectorIndices], function(x) x))))  # put 'expr' back in though shouldn't be needed downstream
                        ## otherwise continue with processing as in non-dynamic index case
                    } else {
                        ## in this case the expr var is a different var, so shouldn't really matter what
                        ## indexes get plugged in; plug in mins of indexes as a default
                        varInfo <- model$modelDef$varInfo[[safeDeparse(expr[[2]], warn = TRUE)]]
                        expr[which(!numericOrVectorIndices)+2] <- varInfo$mins[!numericOrVectorIndices]
                        ## sapply business gets rid of () at end of index expression
                        newExpr <- as.call(c(cc_structureExprName, expr, sapply(indexExprs[!numericOrVectorIndices], function(x) x)))
                        for(i in seq_along(newExpr)[-1])
                            newExpr[[i]] <- cc_expandDetermNodesInExpr(model, newExpr[[i]], targetNode, skipExpansionsNode)
                        return(newExpr)
                    }
                }  ## else continue with processing as in non-dynamic index case
            }
        }
        exprText <- safeDeparse(expr, warn = TRUE)
        expandedNodeNamesRaw <- model$expandNodeNames(exprText)
        if(!is.null(skipExpansionsNode) && (exprText %in% model$expandNodeNames(skipExpansionsNode, returnScalarComponents=TRUE))) return(expr)
        ## if exprText is a node itself (and also part of a larger node), then we only want the expansion to be the exprText node:
        expandedNodeNames <- if(exprText %in% expandedNodeNamesRaw) exprText else expandedNodeNamesRaw
        if(length(expandedNodeNames) == 1 && (expandedNodeNames == exprText)) {
            ## expr is a single node in the model
            type <- model$getNodeType(exprText)
            if(length(type) > 1) {
                ## if exprText is a node itself (and also part of a larger node), then we only want the expansion to be the exprText node:
                if(exprText %in% expandedNodeNamesRaw) type <- type[which(exprText == expandedNodeNamesRaw)]
                else stop('something went wrong with Daniel\'s understanding of newNimbleModel #1')
            }
            if(type == 'stoch') return(expr)
            if(type == 'determ') {
                newExpr <- model$getValueExpr(exprText)
                return(cc_expandDetermNodesInExpr(model, newExpr, targetNode, skipExpansionsNode))
            }
            if(type == 'RHSonly') return(expr)
            stop('something went wrong with Daniel\'s understanding of newNimbleModel #2')
        }
        newExpr <- cc_createStructureExpr(model, exprText)
        for(i in seq_along(newExpr)[-1])
            newExpr[[i]] <- cc_expandDetermNodesInExpr(model, newExpr[[i]], targetNode, skipExpansionsNode)
        return(newExpr)
    }
    if(is.call(expr)) {
        if(any(sapply(expr[-1], function(x) x=='')))
            stop('Found missing indexing in ', safeDeparse(expr), ' that prevents conjugacy processing for this particular model structure.')
        for(i in seq_along(expr)[-1])
            expr[[i]] <- cc_expandDetermNodesInExpr(model, expr[[i]], targetNode, skipExpansionsNode)
        return(expr)
    }
    stop(paste0('something went wrong processing: ', safeDeparse(expr)))
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
    linearityLinks <- c('identity', 'additive', 'multiplicative', 'multiplicativeScalar', 'linear') 
    if(!(link %in%  c(linearityLinks, 'stickbreaking')))
        stop(paste0('unknown link: \'', link, '\''))
    if(is.null(linearityCheck))    return(NULL)
    offset      <- linearityCheck$offset
    scale       <- linearityCheck$scale
    if(link %in% linearityLinks && offset == 0 && scale == 1)
        return('identity')
    if(link %in% c('additive', 'linear') && scale == 1)
        return('additive')
    ## For Wishart/inverse-Wishart, only handle scalar 'scale'.
    if(link == 'multiplicativeScalar' && offset == 0 && cc_checkScalar(scale))
        return('multiplicativeScalar')
    if(link %in% c('multiplicative', 'linear') && offset == 0)
        return('multiplicative')
    if(link == 'linear')
       return('linear')
    return(NULL)
}

cc_checkScalar <- function(expr) {
    if(length(expr) == 1) return(TRUE)
    if(expr[[1]] == ":") return(FALSE)
    ## We don't know the output dims of all functions,
    ## as there can be user-defined functions that produce non-scalar output from scalar inputs,
    ## so only say it could be scalar (based on input args) in specific known cases we enumerate.
    if(safeDeparse(expr[[1]], warn = TRUE) %in% c('(','[','*','+','-','/','exp','log','^','pow','sqrt')) {
        return(all(sapply(expr[2:length(expr)], cc_checkScalar)))
    } else return(FALSE)
}

## Removes 'coeff' and/or 'offset' from expression for use in simplifying conjugacy contribution
## expressions when realized link is 'simpler' than possible links.
cc_stripExpr <- function(expr, offset = TRUE, coeff = TRUE) {
    if(length(expr) == 1) {
        if(coeff && expr == 'coeff') return(1)
        if(offset && expr == 'offset') return(0)
        return(expr)
    }
    if(offset && (expr[[1]] == '+' || expr[[1]] == '-')) {
        if(expr[[2]] == 'offset' || expr[[2]] == 0)
            return(cc_stripExpr(expr[[3]], offset, coeff))
        if(expr[[3]] == 'offset' || expr[[3]] == 0)
            return(cc_stripExpr(expr[[2]], offset, coeff))
    }
    if(coeff && expr[[1]] == '*') {
        if(expr[[2]] == 'coeff' || expr[[2]] == 1)
            return(cc_stripExpr(expr[[3]], offset, coeff))
        if(expr[[3]] == 'coeff' || expr[[3]] == 1)
            return(cc_stripExpr(expr[[2]], offset, coeff))
    }
    if(coeff && expr[[1]] == '/' && (expr[[3]] == 'coeff' || expr[[3]] == 1)) 
        return(cc_stripExpr(expr[[2]], offset, coeff))
    if(coeff && expr[[1]] == '/' && expr[[2]] == 'coeff') {
        expr[[2]] <- 1
        expr[[3]] <- cc_stripExpr(expr[[3]], offset, coeff)
        return(expr)
    }
    if(coeff && (expr[[1]] == '^' || expr[[1]] == 'sqrt') && expr[[2]] == 'coeff') {
        return(1)
    }
    if(is.call(expr)) 
        for(i in seq_along(expr)[-1])
            expr[[i]] <- cc_stripExpr(expr[[i]], offset, coeff)
    return(expr)
}


## checks the parameter expressions in the stochastic distribution of depNode
## returns FALSE if we find 'targetNode' in ***more than one*** of these expressions
cc_otherParamsCheck <- function(model, depNode, targetNode, skipExpansionsNode = NULL, depNodeExprExpanded, depParamNodeName) {
    paramsList <- as.list(model$getValueExpr(depNode)[-1])       # extracts the list of all parameters, for the distribution of depNode
    timesFound <- 0   ## for success, we'll find targetNode in only *one* parameter expression
    for(i in seq_along(paramsList)) {
        if(!missing(depParamNodeName) && (names(paramsList)[i] == depParamNodeName)) {
            expr <- depNodeExprExpanded
        } else { expr <- cc_expandDetermNodesInExpr(model, paramsList[[i]], targetNode, skipExpansionsNode) }
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
## exporting as used in setup code of a nf called directly by user
#' @export
cc_getNodesInExpr <- function(expr) {
    if(is.numeric(expr)) return(character(0))   ## expr is numeric
    if(is.logical(expr)) return(character(0))   ## expr is logical
    if(is.name(expr) || (is.call(expr) && (expr[[1]] == '[') && is.name(expr[[2]]))) return(safeDeparse(expr, warn = TRUE))   ## expr is a node name
    if(is.call(expr)) return(unlist(lapply(expr[-1], cc_getNodesInExpr)))   ## expr is some general call
    stop(paste0('something went wrong processing: ', safeDeparse(expr)))
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

cc_replace01 <- function(expr) {
    ## replace {0,1} with {0,1}+1i so distinguished from offset,scale of 0,1 when match target
    if(is.numeric(expr) && (expr %in% c(0, 1))) 
        return(expr + 1i)
    if(is.call(expr))
        for(i in 2:length(expr))
            expr[[i]] <- cc_replace01(expr[[i]])
    return(expr)
}

cc_checkLinearity <- function(expr, targetNode) {

    ## targetNode doesn't appear in expr
    if(!cc_nodeInExpr(targetNode, expr)) {
        if(is.call(expr) && expr[[1]] == '(') return(cc_checkLinearity(expr[[2]], targetNode))
        # add +1i to tags 0s and 1s as not being from exact match to target
        return(list(offset = cc_replace01(expr), scale = 0))  
    }

    ## expr is exactly the targetNode
    if(identical(targetNode, safeDeparse(expr, warn = TRUE)))
        return(list(offset = 0, scale = 1))

    if(!is.call(expr))   stop('cc_checkLinearity: expression is not a call object')

    ## process the expression contents of the parentheses
    if(expr[[1]] == '(')
        return(cc_checkLinearity(expr[[2]], targetNode))

    if(expr[[1]] == '[') 
        return(cc_checkLinearity(expr[[2]], targetNode))

    ## Look for individual nodes in vectorized use or other strange cases.
    if(expr[[1]] == 'structureExpr') {
        ## Can't have target appear multiple times.
	if(sum(sapply(expr[2:length(expr)], function(x)
            cc_nodeInExpr(targetNode, x))) != 1)
            return(NULL)
        wh <- which(sapply(expr[2:length(expr)], function(x)
            cc_nodeInExpr(targetNode, x)))
        checkLinearityStrucExpr <- cc_checkLinearity(expr[[wh+1]], targetNode)
        if(is.null(checkLinearityStrucExpr)) return(NULL)
        return(list(offset = cc_combineExprsAddition(expr, checkLinearityStrucExpr$offset),  # was expr?,
                    scale = checkLinearityStrucExpr$scale))
    }

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
        stop('cc_checkLinearity: problem with negation expression')
    }

    if(expr[[1]] == '+') {
        checkLinearityLHS <- cc_checkLinearity(expr[[2]], targetNode)
        checkLinearityRHS <- cc_checkLinearity(expr[[3]], targetNode)
        if(is.null(checkLinearityLHS) || is.null(checkLinearityRHS)) return(NULL)
        return(list(offset = cc_combineExprsAddition(checkLinearityLHS$offset, checkLinearityRHS$offset),
                    scale  = cc_combineExprsAddition(checkLinearityLHS$scale,  checkLinearityRHS$scale)))
    }

    if(expr[[1]] == '*' || expr[[1]] == '%*%' || expr[[1]] == 'inprod' || expr[[1]] == 'sum') {
        ## X[,] %*% beta[] where beta[i] are nodes is not standard matrix multiplication, so there is offset and scale
        isMatrixMult <- ifelse(expr[[1]] == '%*%' &&
                               length(expr[[2]]) > 1 && expr[[2]][[1]] != 'structureExpr' &&
                               length(expr[[3]]) > 1 && expr[[3]][[1]] != 'structureExpr',
                               TRUE, FALSE)
        if(expr[[1]] == 'sum') 
            if(length(expr[[2]]) == 3 && expr[[2]][[1]] == '*') {
                tmpExpr <- quote(inprod(a, b))
                tmpExpr[[2]] <- expr[[2]][[2]]
                tmpExpr[[3]] <- expr[[2]][[3]]
                expr <- tmpExpr
            } else {
                return(NULL)  # cases such as sum(p[1:5]); there may be some unusual conjugacy cases here that we don't detect
            }
        if(length(expr) != 3) return(NULL)  # avoid error at next step for unexpected cases
        if(cc_nodeInExpr(targetNode, expr[[2]]) && cc_nodeInExpr(targetNode, expr[[3]])) return(NULL)
        checkLinearityLHS <- cc_checkLinearity(expr[[2]], targetNode)
        checkLinearityRHS <- cc_checkLinearity(expr[[3]], targetNode)
        if(is.null(checkLinearityLHS) || is.null(checkLinearityRHS)) return(NULL)
        if((checkLinearityLHS$scale != 0) && (checkLinearityRHS$scale != 0)) stop('cc_checkLinearity: incompatible scales in * operation')
        if((checkLinearityLHS$scale == 0) && (checkLinearityRHS$scale == 0)) {
            return(list(offset = cc_combineExprsMultiplication(checkLinearityLHS$offset, checkLinearityRHS$offset, isMatrixMult),
                        scale  = 0)) }
        if(checkLinearityLHS$scale != 0) {
            return(list(offset = cc_combineExprsMultiplication(checkLinearityLHS$offset, checkLinearityRHS$offset, isMatrixMult),
                        scale  = cc_combineExprsMultiplication(checkLinearityLHS$scale,  checkLinearityRHS$offset, isMatrixMult))) }
        if(checkLinearityRHS$scale != 0) {
            return(list(offset = cc_combineExprsMultiplication(checkLinearityLHS$offset, checkLinearityRHS$offset, isMatrixMult),
                        scale  = cc_combineExprsMultiplication(checkLinearityLHS$offset, checkLinearityRHS$scale , isMatrixMult))) }
        stop('cc_checkLinearity: something went wrong')
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

cc_checkStickbreaking <- function(expr, targetNode) {
    if(!is.call(expr) || expr[[1]] != 'stick_breaking' || !cc_nodeInExpr(targetNode, expr))
        return(NULL)
    expr <- expr[[2]]
    if(!is.call(expr) || expr[[1]] != 'structureExpr')
        return(NULL)
    offset <- which(targetNode == cc_getNodesInExpr(expr))
    if(length(offset) != 1) return(NULL) else return(list(offset = offset, scale = 1))
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

createDynamicConjugateSamplerName <- function(prior, dependentCounts, dynamicallyIndexed = FALSE) {
    ##depString <- paste0(dependentCounts, names(dependentCounts), collapse='_')  ## including the numbers of dependents
    depString <- paste0(sort(names(dependentCounts)), collapse='_')                     ## without the numbers of each type of dependent node
    paste0('sampler_conjugate_', prior, '_', depString, ifelse(dynamicallyIndexed, '_dynamicDeps', ''))
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





#' Register new conjugate relationships
#'
#' @author Daniel Turek
#' @export
registerConjugacy <- function(conjugacyDef, userEnv = parent.frame(), verbose = nimbleOptions('verbose')) {
    if(missing(conjugacyDef))   return(invisible(NULL))
    ## check syntax of the newConjugacyDef argument: 
    checkConjugacyDefinition(conjugacyDef)
    ## create nimbleUserNamespace$conjugacyRelationshipsObject, if it doesn't exist yet:
    if(!exists('conjugacyRelationshipsObject', nimbleUserNamespace, inherits = FALSE))
        nimbleUserNamespace$conjugacyRelationshipsObject <- conjugacyRelationshipsClass( list() )
    
    

    

    dupl <- nms[nms %in% getAllDistributionsInfo('namesVector', nimbleOnly = TRUE)]
    if(length(dupl)) {
        distributionsInput[dupl] <- NULL
        duplTogether <- paste0(dupl, collapse = ', ')
        messageIfVerbose("  [Warning] Ignoring the following user-supplied distributions as they have the same names as default NIMBLE distributions: ", duplTogether, ". Please rename to avoid the conflict.")
    }

    if(is.list(distributionsInput)) 
        sapply(distributionsInput, checkDistributionInput)
    sapply(distributionsInput, checkDistributionFunctions, userEnv = userEnv)
    if(is.character(distributionsInput)) {
        distributionsInput <- lapply(distributionsInput, prepareDistributionInput, userEnv = userEnv)
        names(distributionsInput) <- nms
    }
    
    if(exists('distributions', nimbleUserNamespace, inherits = FALSE)) {
        nimbleUserNamespace$distributions$add(distributionsInput)
    } else 
        nimbleUserNamespace$distributions <- distributionsClass(distributionsInput)
    virtualNodeFunctionDefinitions <- ndf_createVirtualNodeFunctionDefinitionsList(userAdded = TRUE)
    createNamedObjectsFromList(virtualNodeFunctionDefinitions, envir = .GlobalEnv)

    # note don't use rFunHandler as rUserDist nimbleFunction needs n as first arg so it works on R side, therefore we have n in the C version of the nimbleFunction and don't want to strip it out in Cpp generation

    return(invisible(NULL))
}


checkConjugacyDefinition <- function(newConjugacyDef) {
    message('XXX need to implement checking of newConjugacyDef')
}
