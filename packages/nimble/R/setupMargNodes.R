#' Organize model nodes for marginalization
#'
#' Process model to organize nodes for marginalization (integration over latent 
#' nodes or random effects) as by Laplace approximation.
#'
#' @param model A nimble model such as returned by \code{nimbleModel}.
#'
#' @param paramNodes A character vector of names of stochastic nodes that are
#'   parameters of nodes to be marginalized over (\code{randomEffectsNodes}).
#'   See details for default.
#'
#' @param randomEffectsNodes A character vector of nodes to be marginalized over
#'   (or "integrated out"). In the case of calculating the likelihood of a model
#'   with continuous random effects, the nodes to be marginalized over are the
#'   random effects, hence the name of this argument. However, one can
#'   marginalize over any nodes desired as long as they are continuous. 
#'   See details for default.
#'
#' @param calcNodes A character vector of nodes to be calculated as the
#'   integrand for marginalization. Typically this will include
#'   \code{randomEffectsNodes} and some data nodes. Se details for default.
#'
#' @param calcNodesOther A character vector of nodes to be calculated as part of
#'   the log likelihood that are not connected to the \code{randomEffectNodes}
#'   and so are not actually part of the marginalization. These are somewhat
#'   extraneous to the purpose of this function, but it is convenient to handle
#'   them here because often the purpose of marginalization is to calculate log
#'   likelihoods, including from "other" parts of the model.
#'
#' @param split A logical indicating whether to split \code{randomEffectsNodes}
#'   into conditionally independent sets that can be marginalized separately
#'   (\code{TRUE}) or to keep them all in one set for a single marginalization
#'   calculation.
#'
#' @param check A logical indicating whether to try to give reasonable warnings
#'   of badly formed inputs that might be missing important nodes or include
#'   unnecessary nodes.
#'
#' @param allowDiscreteLatent A logical indicating whether to
#'   allow discrete latent states. (default = \code{FALSE})
#'
#' @details
#'
#' This function is used by \code{buildLaplace} to organize model nodes into
#' roles needed for setting up the (approximate) marginalization done by Laplace
#' approximation. It is also possible to call this function directly and pass
#' the resulting list (possibly modified for your needs) to \code{buildLaplace}.
#'
#' Any of the input node vectors, when provided, will be processed using
#'   \code{nodes <- model$expandNodeNames(nodes)}, where \code{nodes} may be
#'   \code{paramNodes}, \code{randomEffectsNodes}, and so on. This step allows
#'   any of the inputs to include node-name-like syntax that might contain
#'   multiple nodes. For example, \code{paramNodes = 'beta[1:10]'} can be
#'   provided if there are actually 10 scalar parameters, 'beta[1]' through
#'   'beta[10]'. The actual node names in the model will be determined by the
#'   \code{exapndNodeNames} step.
#'
#' This function does not do any of the marginalization calculations. It only
#' organizes nodes into roles of parameters, random effects, integrand
#' calculations, and other log likelihood calculations.
#'
#' The checking done if `check=TRUE` tries to be reasonable, but it can't cover
#' all cases perfectly. If it gives an unnecessary warning, simply set `check=FALSE`.
#'
#' If \code{paramNodes} is not provided, its default depends on what other
#'   arguments were provided. If neither \code{randomEffectsNodes} nor
#'   \code{calcNodes} were provided, \code{paramNodes} defaults to all
#'   top-level, stochastic nodes, excluding any posterior predictive nodes
#'   (those with no data anywhere downstream). These are determined by
#'   \code{model$getNodeNames(topOnly = TRUE, stochOnly = TRUE,
#'   includePredictive = FALSE)}. If \code{randomEffectsNodes} was provided,
#'   \code{paramNodes} defaults to stochastic parents of
#'   \code{randomEffectsNodes}. In these cases, any provided \code{calcNodes} or
#'   \code{calcNodesOther} are excluded from default \code{paramNodes}. If
#'   \code{calcNodes} but not \code{randomEffectsNodes} was provided, then the
#'   default for \code{randomEffectsNodes} is determined first, and then
#'   \code{paramNodes} defaults to stochastic parents of
#'   \code{randomEffectsNodes}. Finally, any stochastic parents of
#'   \code{calcNodes} (whether provided or default) that are not in
#'   \code{calcNodes} are added to the default for \code{paramNodes}, but only
#'   after \code{paramNodes} has been used to determine the defaults for
#'   \code{randomEffectsNodes}, if necessary.
#'
#' Note that to obtain sensible defaults, some nodes must have been marked as
#'   data, either by the \code{data} argument in \code{nimbleModel} or by
#'   \code{model$setData}. Otherwise, all nodes will appear to be posterior
#'   predictive nodes, and the default \code{paramNodes} may be empty.
#'
#' For purposes of \code{buildLaplace}, \code{paramNodes} does not need to (but
#'   may) include deterministic nodes between the parameters and any
#'   \code{calcNodes}. Such deterministic nodes will be included in
#'   calculations automatically when needed.
#'
#' If \code{randomEffectsNodes} is missing, the default is a bit complicated: it
#'   includes all latent nodes that are descendants (or "downstream") of
#'   \code{paramNodes} (if provided) and are either (i) ancestors (or
#'   "upstream") of data nodes (if \code{calcNodes} is missing), or (ii)
#'   ancestors or elements of \code{calcNodes} (if \code{calcNodes} and
#'   \code{paramNodes} are provided), or (iii) elements of \code{calcNodes} (if
#'   \code{calcNodes} is provided but \code{paramNodes} is missing). In all
#'   cases, discrete nodes (with warning if \code{check=TRUE}), posterior
#'   predictive nodes and \code{paramNodes} are excluded.
#'
#' \code{randomEffectsNodes} should only include stochastic nodes.
#'
#' If \code{calcNodes} is missing, the default is \code{randomEffectsNodes} and
#'   their descendants to the next stochastic nodes, excluding posterior
#'   predictive nodes. These are determined by
#'   \code{model$getDependencies(randomEffectsNodes, includePredictive=FALSE)}.
#'
#' If \code{calcNodesOther} is missing, the default is all stochastic
#'   descendants of \code{paramNodes}, excluding posterior predictive nodes
#'   (from \code{model$getDependencies(paramNodes, stochOnly=TRUE, self=FALSE,
#'   includePosterior=FALSE)}) that are not part of \code{calcNodes}.
#'
#' For purposes of \code{buildLaplace}, neither \code{calcNodes} nor
#'   \code{calcNodesOther} needs to (but may) contain deterministic nodes
#'   between \code{paramNodes} and \code{calcNodes} or \code{calcNodesOther},
#'   respectively. These will be included in calculations automatically when
#'   needed.
#'
#' If \code{split} is \code{TRUE}, \code{model$getConditionallyIndependentSets}
#'   is used to determine sets of the \code{randomEffectsNodes} that can be
#'   independently marginalized. The \code{givenNodes} are the
#'   \code{paramNodes} and \code{calcNodes} excluding any
#'   \code{randomEffectsNodes} and their deterministic descendants. The
#'   \code{nodes} (to be split into sets) are the \code{randomEffectsNodes}.
#'
#' If \code{split} is a numeric vector, \code{randomEffectsNodes} will be split
#'   by \code{split}(\code{randomEffectsNodes}, \code{control$split}). The last
#'   option allows arbitrary control over how \code{randomEffectsNodes} are
#'   blocked.
#'
#' If \code{check=TRUE}, then defaults for each of the four categories of nodes
#'   are created even if the corresponding argument was provided. Then warnings
#'   are emitted if there are any extra (potentially unnecessary) nodes provided
#'   compared to the default or if there are any nodes in the default that were
#'   not provided (potentially necessary). These checks are not perfect and may
#'   be simply turned off if you are confident in your inputs.
#'
#' (If \code{randomEffectsNodes} was provided but \code{calcNodes} was not
#'   provided, the default (for purposes of \code{check=TRUE} only) for
#'   \code{randomEffectsNodes} differs from the above description. It uses
#'   stochastic descendants of \code{randomEffectsNodes} in place of the
#'   "data nodes" when determining ancestors of data nodes. And it uses item
#'   (ii) instead of (iii) in the list above.)
#'
#' @author Wei Zhang, Perry de Valpine, Paul van Dam-Bates
#' @return
#'
#' A list is returned with elements:
#'
#' \itemize{
#'
#' \item \code{paramNodes}: final processed version of \code{paramNodes}
#'
#' \item \code{randomEffectsNodes}: final processed version of \code{randomEffectsNodes}
#'
#' \item \code{calcNodes}: final processed version of \code{calcNodes}
#'
#' \item \code{calcNodesOther}: final processed version of \code{calcNodesOther}
#'
#' \item \code{givenNodes}: Input to \code{model$getConditionallyIndependentSets}, if \code{split=TRUE}.
#'
#' \item \code{randomEffectsSets}: Output from
#'   \code{model$getConditionallyIndependentSets}, if \code{split=TRUE}. This
#'   will be a list of vectors of node names. The node names in one list element
#'   can be marginalized independently from those in other list elements. The
#'   union of the list elements should be all of \code{randomEffectsNodes}. If
#'   \code{split=FALSE}, \code{randomEffectsSets} will be a list with one
#'   element, simply containing \code{randomEffectsNodes}. If \code{split} is a
#'   numeric vector,  \code{randomEffectsSets} will be the result of
#'   \code{split}(\code{randomEffectsNodes}, \code{control$split}).
#'
#' }
#'
#' @export
setupMargNodes <- function(model, paramNodes, randomEffectsNodes, calcNodes,
                           calcNodesOther,
                           split = TRUE,
                           check = TRUE,
                           allowDiscreteLatent = FALSE) {
  paramProvided     <- !missing(paramNodes)
  reProvided        <- !missing(randomEffectsNodes)
  calcProvided      <- !missing(calcNodes)
  calcOtherProvided <- !missing(calcNodesOther)

  normalizeNodes <- function(nodes, sort = FALSE) {
    if(is.null(nodes) || isFALSE(nodes)) character(0)
    else model$expandNodeNames(nodes, sort = sort)
  }
  if(paramProvided) paramNodes         <- normalizeNodes(paramNodes)
  if(reProvided)    randomEffectsNodes <- normalizeNodes(randomEffectsNodes)
  if(calcProvided)  calcNodes          <- normalizeNodes(calcNodes, sort = TRUE)
  if(calcOtherProvided) calcNodesOther <- normalizeNodes(calcNodesOther, sort = TRUE)

  if(reProvided) {
    if(check && !allowDiscreteLatent)
      if(any(model$isDiscrete(randomEffectsNodes)))
        messageIfVerbose("  [Warning] Some elements of `randomEffectsNodes` follow discrete distributions. That is likely to cause problems.")
  }

  # We considered a feature to allow params to be nodes without priors. This is a placeholder in case
  # we ever pursue that again.
  # allowNonPriors <- FALSE
  # We may need to use determ and stochastic dependencies of parameters multiple times below
  # Define these to avoid repeated computation
  # A note for future: determ nodes between parameters and calcNodes are needed inside buildOneAGHQuad
  # and buildOneAGHQuad1D. In the future, these could be all done here to be more efficient
  paramDetermDeps <- character(0)
  paramStochDeps  <- character(0)
  paramDetermDepsCalculated <- FALSE
  paramStochDepsCalculated  <- FALSE

  # 1. Default parameters are stochastic top-level nodes. (We previously
  #    considered an argument allowNonPriors, defaulting to FALSE. If TRUE, the
  #    default params would be all top-level stochastic nodes with no RHSonly
  #    nodes as parents and RHSonly nodes (handling of constants TBD, since
  #    non-scalars would be converted to data) that have stochastic dependencies
  #    (And then top-level stochastic nodes with RHSonly nodes as parents are
  #    essentially latent/data nodes, some of which would need to be added to
  #    randomEffectsNodes below.) However this got too complicated. It is
  #    simpler and clearer to require "priors" for parameters, even though prior
  #    probs may not be used.
  paramsHandled <- TRUE
  if(!paramProvided) {
    if(!reProvided) {
      if(!calcProvided) {
        paramNodes <- model$getNodeNames(topOnly = TRUE, stochOnly = TRUE, includePredictive = FALSE)
      } else {
        # calcNodes were provided, but RE nodes were not, so delay creating default params
        paramsHandled <- FALSE
      }
    } else {
      nodesToFindParentsFrom <- randomEffectsNodes
      paramNodes <- model$getParents(nodesToFindParentsFrom, self=FALSE, stochOnly=TRUE)
      # self=FALSE doesn't omit if one RE node is a parent of another, so we have to do the next step
      paramNodes <- setdiff(paramNodes, nodesToFindParentsFrom)
    }
    if(paramsHandled) {
      if(calcProvided) paramNodes <- setdiff(paramNodes, calcNodes)
      if(calcOtherProvided) paramNodes <- setdiff(paramNodes, calcNodesOther)
    }
  }

  # 2. Default random effects are latent nodes that are downstream stochastic dependencies of params.
  #    In step 3, default random effects are also limited to those that are upstream parents of calcNodes
  if((!reProvided) || check) {
    latentNodes <- model$getNodeNames(latentOnly = TRUE, stochOnly = TRUE,
                                      includeData = FALSE, includePredictive = FALSE)
    if(!allowDiscreteLatent) {
      latentDiscrete <- model$isDiscrete(latentNodes)
      if(any(latentDiscrete)) {
        if((!reProvided) && check) {
          messageIfVerbose("  [Note] In trying to determine default `randomEffectsNodes`, there are some nodes\n",
                  "         that follow discrete distributions. These will be omitted.")
        }
        latentNodes <- latentNodes[!latentDiscrete]
      }
    }
    if(paramsHandled) {
      paramDownstream <- model$getDependencies(paramNodes, stochOnly = TRUE, self = FALSE,
                                               downstream = TRUE, includePredictive = FALSE)
      #    paramStochDeps <- model$getDependencies(paramNodes, stochOnly = TRUE, self = FALSE)
      #    paramStochDepsCalculated <- TRUE
      reNodesDefault <- intersect(latentNodes, paramDownstream)
    } else {
      reNodesDefault <- latentNodes
    }
    # Next, if calcNodes were not provided, we create a temporary
    # dataNodesDefault for purposes of updating reNodesDefault if needed. The
    # idea is that reNodesDefault should be trimmed to include only nodes
    # upstream of "data" nodes, where "data" means nodes in the role of data for
    # purposes of marginalization.
    # The tempDataNodesDefault is either dependencies of RE nodes if provided, or
    # actual data nodes in the model if RE nodes not provided.
    # If calcNodes were provided, then they are used directly to trim reNodesDefault.
    if(!calcProvided) {
      if(reProvided)
        tempDataNodesDefault <- model$getDependencies(randomEffectsNodes, stochOnly = TRUE,
                                                      self = FALSE, includePredictive = FALSE)
      else
        tempDataNodesDefault <- model$getNodeNames(dataOnly = TRUE)
      if(paramsHandled)
        tempDataNodesDefault <- setdiff(tempDataNodesDefault, paramNodes)
      tempDataNodesDefaultParents <- model$getParents(tempDataNodesDefault, upstream = TRUE, stochOnly = TRUE)
      # See comment above about why this is necessary:
      tempDataNodesDefaultParents <- setdiff(tempDataNodesDefaultParents, tempDataNodesDefault)
      reNodesDefault <- intersect(reNodesDefault, tempDataNodesDefaultParents)
    } else {
      # Update reNodesDefault to exclude nodes that lack downstream connection to a calcNode
      if(paramsHandled) { # This means reProvided OR paramsProvided. Including parents allows checking
        # of potentially missing REs.
        reNodesDefault <- intersect(reNodesDefault,
                                    model$getParents(calcNodes, upstream=TRUE, stochOnly = TRUE))
      } else { # This means !paramsHandled and hence !reProvided AND !paramsProvided
        reNodesDefault <- intersect(reNodesDefault,
                                    calcNodes)
        reNodesDefault <- intersect(reNodesDefault,
                                    model$getParents(calcNodes, upstream=TRUE, stochOnly = TRUE))
      }
    }
  }

  # If only calcNodes were provided, we have now created reNodesDefault from calcNodes,
  # and are now ready to create default paramNodes
  if(!paramsHandled) {
    paramNodes <- model$getParents(reNodesDefault, self=FALSE, stochOnly=TRUE)
    # See comment above about why this is necessary:
    paramNodes <- setdiff(paramNodes, reNodesDefault)
    if(calcOtherProvided) paramNodes <- setdiff(paramNodes, calcNodesOther)
  }

  # 3. Optionally check random effects if they were provided (not default)
  if(reProvided && check) {
    # First check is for random effects that should have been included but weren't
    reCheck <- setdiff(reNodesDefault, randomEffectsNodes)
    if(length(reCheck)) {
      errorNodes <- paste0(head(reCheck, n = 4), sep = "", collapse = ", ")
      if(length(reCheck) > 4) errorNodes <- paste(errorNodes, "...")
      messageIfVerbose("  [Warning] There are some random effects (latent states) in the model that look like\n",
                       "            they should be included for the provided (or default) `paramNodes`,\n",
                       "            but are not included in `randomEffectsNodes`: ", errorNodes, ".\n",
                       "            To silence this warning, one can usually include `check = FALSE`\n",
                       "            (potentially in the control list) for the algorithm or as\n",
                       "            an argument to `setupMargNodes`.")
    }
    # Second check is for random effects that were included but look unnecessary
    reCheck <- setdiff(randomEffectsNodes, reNodesDefault)
    if(length(reCheck)) {
      # Top nodes should never trigger warning.
      # Descendants of top nodes that are in randomEffectsNodes should not trigger warning
      topNodes <- model$getNodeNames(topOnly=TRUE)
      reCheckTopNodes <- intersect(reCheck, topNodes)
      if(length(reCheckTopNodes)) {
        # Simple downstream=TRUE here is not a perfect check of connection among all nodes
        # but it will avoid false alarms
        reCheck <- setdiff(reCheck, model$getDependencies(reCheckTopNodes, downstream=TRUE, stochOnly=TRUE))
      }
      if(length(reCheck)) {
        errorNodes <- paste0(head(reCheck, n = 4), sep = "", collapse = ", ")
        if(length(reCheck) > 4) errorNodes <- paste(errorNodes, "...")
        extraMsg <- if(isTRUE(getNimbleOption('includeUnneededLatents'))) "" else "            They will be omitted, but one can force inclusion with\n            `nimbleOptions(includeUnneededLatents=TRUE)`.\n"
        messageIfVerbose("  [Warning] There are some `randomEffectsNodes` provided that look like\n",
                         "            they are not needed for the provided (or default) `paramNodes`:\n",
                         "            ", errorNodes, ".\n", extraMsg,
                         "            To silence this warning, one can usually include `check = FALSE`\n",
                         "            (potentially in the control list) for the algorithm or as\n",
                         "            an argument to `setupMargNodes`.")
        if(!isTRUE(getNimbleOption('includeUnneededLatents')))
            randomEffectsNodes <- setdiff(randomEffectsNodes, reCheck)
      }
    }
  }
  # Set final choice of randomEffectsNodes
  if(!reProvided) {
    randomEffectsNodes <- reNodesDefault
  }

  # Set actual default calcNodes. This time it has self=TRUE (default)
  if((!calcProvided) || check) {
    calcNodesDefault <- model$getDependencies(randomEffectsNodes, includePredictive = FALSE)
  }
  # 5. Optionally check calcNodes if they were provided (not default)
  if(calcProvided && check) {
    # First check is for calcNodes that look necessary but were omitted
    calcCheck <- setdiff(calcNodesDefault, calcNodes)
    if(length(calcCheck)) {
      errorNodes <- paste0(head(calcCheck, n = 4), sep = "", collapse = ", ")
      if(length(calcCheck) > 4) errorNodes <- paste(errorNodes, "...")
      messageIfVerbose("  [Warning] There are some model nodes that look like they should be\n",
                       "            included in the `calcNodes` because\n",
                       "            they are dependencies of some `randomEffectsNodes`: ", errorNodes, ".\n",
                       "            To silence this warning, one can usually include `check = FALSE`\n",
                       "            (potentially in the control list) for the algorithm or as\n",
                       "            an argument to `setupMargNodes`.")
    }
    # Second check is for calcNodes that look unnecessary
    # If some determ nodes between paramNodes and randomEffectsNodes are provided in calcNodes
    # then that's ok and we should not throw a warning message.
    calcCheck <- setdiff(calcNodes, calcNodesDefault)
    errorNodes <- calcCheck[model$getNodeType(calcCheck)=="stoch"]
    # N.B. I commented out this checking of deterministic nodes for now.
    #      Iterating through individual nodes for getDependencies can be slow
    #      and I'd like to think more about how to do this. -Perry
    ## determCalcCheck <- setdiff(calcCheck, errorNodes)
    ## lengthDetermCalcCheck <- length(determCalcCheck)
    ## # Check other determ nodes
    ## if(lengthDetermCalcCheck){
    ##   paramDetermDeps <- model$getDependencies(paramNodes, determOnly = TRUE, includePredictive = FALSE)
    ##   paramDetermDepsCalculated <- TRUE
    ##   for(i in 1:lengthDetermCalcCheck){
    ##     if(!(determCalcCheck[i] %in% paramDetermDeps) ||
    ##        !(any(model$getDependencies(determCalcCheck[i], self = FALSE) %in% calcNodesDefault))){
    ##       errorNodes <- c(errorNodes, determCalcCheck[i])
    ##     }
    ##   }
    ## }
    if(length(errorNodes)){
      outErrorNodes <- paste0(head(errorNodes, n = 4), sep = "", collapse = ", ")
      if(length(errorNodes) > 4) outErrorNodes <- paste(outErrorNodes, "...")
      messageIfVerbose("  [Warning] There are some `calcNodes` provided that look like\n",
                       "            they are not needed for the provided (or default) `randomEffectsNodes`:\n",
                       "            ", outErrorNodes, ".\n",
                       "            To silence this warning, one can usually include `check = FALSE`\n",
                       "            (potentially in the control list) for the algorithm or as\n",
                       "            an argument to `setupMargNodes`.")
    }
  }
  # Finish step 4
  if(!calcProvided){
    calcNodes <- calcNodesDefault
  }
  if(!paramProvided) {
    possibleNewParamNodes <- model$getParents(calcNodes, self=FALSE, stochOnly=TRUE, includeData=FALSE)
    # includeData=FALSE as data nodes cannot be parameters
    # self=FALSE doesn't omit if one node is a parent of another, so we have to do the next step
    possibleNewParamNodes <- setdiff(possibleNewParamNodes, calcNodesDefault)
    paramNodes <- unique(c(paramNodes, possibleNewParamNodes))
  }

  # 6. Default calcNodesOther: nodes needed for full model likelihood but
  #    that are not involved in the marginalization done by Laplace.
  #    Default is a bit complicated: All dependencies from paramNodes to
  #    stochastic nodes that are not part of calcNodes. Note that calcNodes
  #    does not necessarily contain deterministic nodes between paramNodes and
  #    randomEffectsNodes. We don't want to include those in calcNodesOther.
  #    (A deterministic that is needed for both calcNodes and calcNodesOther should be included.)
  #    So we have to first do a setdiff on stochastic nodes and then fill in the
  #    deterministics that are needed.
  if(!calcOtherProvided || check) {
    paramStochDeps <- model$getDependencies(paramNodes, stochOnly = TRUE, # Should this be dataOnly=TRUE?
                                            self = FALSE, includePredictive = FALSE)
    calcNodesOtherDefault <- setdiff(paramStochDeps, calcNodes)
  }
  if(calcOtherProvided) {
    if((length(calcNodesOther) > 0) && !any(model$getNodeType(calcNodesOther)=="stoch")){
      messageIfVerbose("  [Warning] There are no stochastic nodes in the `calcNodesOther` provided for Laplace or AGHQ approximation.")
    }
  }
  if(!calcOtherProvided){
    calcNodesOther <- calcNodesOtherDefault
  }
  if(calcOtherProvided && check) {
    calcOtherCheck <- setdiff(calcNodesOtherDefault, calcNodesOther)
    if(length(calcOtherCheck)) {
      # We only check missing stochastic nodes; determ nodes will be added below
      missingStochNodesInds <- which((model$getNodeType(calcOtherCheck)) == "stoch")
      lengthMissingStochNodes <- length(missingStochNodesInds)
      if(lengthMissingStochNodes){
        missingStochNodes <- calcOtherCheck[missingStochNodesInds]
        errorNodes <- paste0(head(missingStochNodes, n = 4), sep = "", collapse = ", ")
        if(lengthMissingStochNodes > 4) errorNodes <- paste(errorNodes, "...")
        messageIfVerbose("  [Warning] There are some model nodes (stochastic) that look like they should be\n",
                         "            included in the `calcNodesOther` for parts of the likelihood calculation\n",
                         "            outside of Laplace or AGHQ approximation: ", errorNodes, ".\n",
                         "            To silence this warning, include `check = FALSE` in the control list\n",
                         "            to `buildLaplace` or as an argument to `setupMargNodes`.")
      }
    }
    # Check redundant stochastic nodes
    calcOtherCheck <- setdiff(calcNodesOther, calcNodesOtherDefault)
    stochCalcOtherCheck <- calcOtherCheck[model$getNodeType(calcOtherCheck) == "stoch"]
    errorNodes <- stochCalcOtherCheck
    # Check redundant determ nodes
    # N.B. I commented-out this deterministic node checking for reasons similar to above. -Perry
    ## determCalcOtherCheck <- setdiff(calcOtherCheck, stochCalcOtherCheck)
    ## lengthDetermCalcOtherCheck <- length(determCalcOtherCheck)
    ## errorNodes <- character(0)
    ## if(lengthDetermCalcOtherCheck){
    ##   if(!paramDetermDepsCalculated) {
    ##     paramDetermDeps <- model$getDependencies(paramNodes, determOnly = TRUE, includePredictive = FALSE)
    ##     paramDetermDepsCalculated <- TRUE
    ##   }
    ##   for(i in 1:lengthDetermCalcOtherCheck){
    ##     if(!(determCalcOtherCheck[i] %in% paramDetermDeps) ||
    ##        !(any(model$getDependencies(determCalcOtherCheck[i], self = FALSE) %in% calcNodesOtherDefault))){
    ##       errorNodes <- c(errorNodes, determCalcOtherCheck[i])
    ##     }
    ##   }
    ## }
    ## errorNodes <- c(stochCalcOtherCheck, errorNodes)
    if(length(errorNodes)){
      outErrorNodes <- paste0(head(errorNodes, n = 4), sep = "", collapse = ", ")
      if(length(errorNodes) > 4) outErrorNodes <- paste(outErrorNodes, "...")
      messageIfVerbose("  [Warning] There are some nodes provided in `calcNodesOther` that look like\n",
                       "            they are not needed for parts of the likelihood calculation\n",
                       "            outside of Laplace or AGHQ approximation: ", outErrorNodes, ".\n",
                       "            To silence this warning, include `check = FALSE` in the control list\n",
                       "            to `buildLaplace` or as an argument to `setupMargNodes`.")
    }
  }
  # Check and add necessary (upstream) deterministic nodes into calcNodesOther
  # This ensures that deterministic nodes between paramNodes and calcNodesOther are used.
  num_calcNodesOther <- length(calcNodesOther)
  if(num_calcNodesOther > 0){
    if(!paramDetermDepsCalculated) {
      paramDetermDeps <- model$getDependencies(paramNodes, determOnly = TRUE, includePredictive = FALSE)
      paramDetermDepsCalculated <- TRUE
    }
    numParamDetermDeps <- length(paramDetermDeps)
    if(numParamDetermDeps > 0) {
      keep_paramDetermDeps <- logical(numParamDetermDeps)
      for(i in seq_along(paramDetermDeps)) {
        nextDeps <- model$getDependencies(paramDetermDeps[i])
        keep_paramDetermDeps[i] <- any(nextDeps %in% calcNodesOther)
      }
      paramDetermDeps <- paramDetermDeps[keep_paramDetermDeps]
    }
    calcNodesOther <- model$expandNodeNames(c(paramDetermDeps, calcNodesOther), sort = TRUE)
  }

  # 7. Do the splitting into sets (if given) or conditionally independent sets (if TRUE)
  givenNodes <- NULL
  reSets <- list()
  if(length(randomEffectsNodes)) {
    if(isFALSE(split)) {
      reSets <- list(randomEffectsNodes)
    } else {
      if(isTRUE(split)) {
        # givenNodes should only be stochastic
        givenNodes <- setdiff(c(paramNodes, calcNodes),
                              c(randomEffectsNodes,
                                model$getDependencies(randomEffectsNodes, determOnly=TRUE)))
        reSets <- model$getConditionallyIndependentSets(
          nodes = randomEffectsNodes, givenNodes = givenNodes,
          unknownAsGiven = TRUE)
      }
      else if(is.numeric(split)){
        reSets <- split(randomEffectsNodes, split)
      }
      else stop("setupMargNodes: Invalid value for `split`")
    }
  }
  list(paramNodes = paramNodes,
       randomEffectsNodes = randomEffectsNodes,
       calcNodes = calcNodes,
       calcNodesOther = calcNodesOther,
       givenNodes = givenNodes,
       randomEffectsSets = reSets
       )
}


splitLatents <- function(model, paramNodes, latentNodes, calcNodes, calcNodesOther,
    control = list()) {
    stochNodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
    discreteStochNodes <- model$isDiscrete(stochNodes)
    if (any(discreteStochNodes))
        stop("splitLatents: found discrete non-data stochastic nodes in processing nodes for quadrature-based posterior approximation: ",
            paste0(stochNodes[discreteStochNodes], collapse = ", "), ". Discrete non-data stochastic nodes cannot be handled by the posterior approximation algorithm.")
    split <- extractControlElement(control, "split", TRUE)
    check <- extractControlElement(control, "check", TRUE)
    margNodes <- setupMargNodes(model = model, paramNodes = paramNodes, randomEffectsNodes = latentNodes,
        calcNodes = calcNodes, calcNodesOther = calcNodesOther, split = split, check = check)
    if (missing(paramNodes) && missing(latentNodes)) {
        if (!missing(calcNodes) || !missing(calcNodesOther))
            messageIfVerbose("   [Note] Ignoring provide `calcNodes` and `calcNodesOther` because `paramNodes` and `latentNodes` not provided and are being determined automatically.")
        paramNodes <- margNodes$paramNodes
        latentNodes <- margNodes$randomEffectsNodes
        deps <- model$getDependencies(latentNodes, includeData = FALSE, self = FALSE)
        ## By default, we treat "siblings" of latent nodes as latents.
        ## This attempts to have fixed effects in latents,
        ## along with random effects.
        newLatents <- model$getParents(deps, stochOnly = TRUE, includeData = FALSE)
        paramNodes <- setdiff(paramNodes, newLatents)
        latentNodes <- unique(c(latentNodes, newLatents))
        margNodes <- setupMargNodes(model = model, paramNodes = paramNodes,
            randomEffectsNodes = latentNodes, split = split, check = check)
    }

    return(margNodes)
}

