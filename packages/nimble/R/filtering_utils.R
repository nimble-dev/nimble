particleFilter_splitModelSteps <- function(model,
                                           nodes,
                                           iNode,
                                           notFirst) {
  ## if !notFirst (it is the first node), prevNode will not be used anyway, so
  ## giving its value is a dummy.
  prevNode <- nodes[if(notFirst) iNode-1 else iNode]
  thisNode <- nodes[iNode]
  thisNodeExpanded <- model$expandNodeNames(thisNode, sort = TRUE)
  ## Set up steps for calculations internal to thisNode and for downstream simulation from thisNode
  thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
  if(length(thisDeterm) > 0) {
    thisDeterm_is_intermediate <- logical(length(thisDeterm))
    for(i in seq_along(thisDeterm)) {
      theseDeps <- model$getDependencies(thisDeterm[i], stochOnly = TRUE)
      thisDeterm_is_intermediate[i] <- any(theseDeps %in% thisNodeExpanded)
    }
    thisDeterm_self <- thisDeterm[ thisDeterm_is_intermediate ]
    thisDeterm <- thisDeterm[ !thisDeterm_is_intermediate ]
    calc_thisNode_self <-  model$expandNodeNames(c(thisNodeExpanded, thisDeterm_self), ## only for the sort
                                                 sort = TRUE)

  } else {
    calc_thisNode_self <- thisNodeExpanded
  }
  ## calc_thisNode_deps <- all_deps_from_thisNode[-(1:finalSelfIndex)] ## some of these could be unnecessary...
  
  ## Roughly we have prevNode -> prevDeterm -> thisNode -> thisDeterm -> thisData
  ## However, there is further sorting to do.
  ## prevDeterm may have (i) components not really needed for thisNode (leading to other parts of the graph).
  ##                     (ii) components that also depend on one part of thisNode, on which other parts of thisNode might depend.
  ## thisNode might have internal ordering and might comprise multiple nodes.
  ## thisDeterm might have (i) components not really needed for thisData (leading to other parts of the graph).
  ##                       (ii) components that also depend on one part of thisData, on which other parts of thisData might depend.
  ## thisData should be ordered and complete.
  ##
  ## The bootFstep only ever does thisDeterm and thisData in sequence, so they can be combined. Same for auxFstep
  ##
  ## In bootFstep, the prevDeterm and thisNode could be combined, but that is not the case in auxFstep,
  ## because it breaks apart thisNode for doing the lookahead
  prevDeterm <- model$getDependencies(prevNode, determOnly = TRUE)
  ## Weed out any prevDeterm elements that do not lead to part of thisNode
  if(length(prevDeterm) > 0) {
    keep_prevDeterm <- logical(length(prevDeterm))
    for(i in seq_along(prevDeterm)) {
      theseDeps <- model$getDependencies(prevDeterm[i])
      keep_prevDeterm[i] <- any(theseDeps %in% calc_thisNode_self)
    }
    prevDeterm <- prevDeterm[ keep_prevDeterm ]
  }
  ## Weed out any prevDeterm elements that are redundant with deterministic parts of calc_thisNode_self
  ##
  if(length(prevDeterm) > 0) {
    keep_prevDeterm <- logical(length(prevDeterm))
    for(i in seq_along(prevDeterm)) {
      theseDeps <- model$getDependencies(prevDeterm[i], determOnly = TRUE)
      keep_prevDeterm[i] <- !any(theseDeps %in% calc_thisNode_self)
    }
    prevDeterm <- prevDeterm[ keep_prevDeterm ]
  }
  
  ## So far we have prevNode -> prevDeterm -> calc_thisNode_self -> calc_thisNode_deps
  ## calc_thisNode_deps combines the old thisDeterm and thisData.
  ## Separating out determ and data calcs is not necessary, but we do need to skip calcs
  ## that lead to the next latent nodes.
  ## A question is whether we condition only on data -- I think yes.
  
  ## Do similar weeding of calc_thisNode_deps in relation to thisData.
  ## calc_thisNode_deps includes determ+data dependencies of thisNode but excludes intermediate deterministic
  ##    nodes that need to be calculated between multiple nodes of thisNode.
  ## We will separately get determ and data dependencies,
  ##    filter determ by which ones lead to data (vs. other states)
  ## The list of determ will be only deterministic nodes that are in calc_thisNode_deps
  ## What if there is a data node in calc_thisNode_self?
  ## thisDeterm <- model$getDependencies(thisNode, determOnly = TRUE)
  ## from above, thisDeterm are determ deps on thisNode that are not needed
  ##    internally, i.e. between elements of thisNode
  thisData   <- model$getDependencies(thisNode, dataOnly = TRUE)
  if(length(thisDeterm) > 0) {
    keep_thisDeterm <- logical(length(thisDeterm))
    for(i in seq_along(thisDeterm)) {
     theseDeps <- model$getDependencies(thisDeterm[i])
     keep_thisDeterm[i] <- any(theseDeps %in% thisData)
    }
    thisDeterm <- thisDeterm[keep_thisDeterm]
    calc_thisNode_deps <- model$expandNodeNames(c(thisDeterm, thisData), sort = TRUE) ## only for the sort
  } else {
    calc_thisNode_deps <- thisData
  }
  list(prevDeterm = prevDeterm,
       calc_thisNode_self = calc_thisNode_self,
       calc_thisNode_deps = calc_thisNode_deps)
}

fillIndices <- function(node, info, returnExpr = FALSE) {
    ## Fill missing indexes with full extent of that dimension.
    node <- parse(text = node)[[1]]
    if(info$nDim != length(node) - 2 || node[[1]] != '[')
        stop("findLatentNodes: invalid node expression: ", node, ".")
    for(i in seq_len(info$nDim)) { 
        if(node[[i+2]] == "") {
            node[[i+2]] <- substitute(A:B,
                                      list(A = info$mins[i], B = info$maxs[i]))
            if(info$mins[i] == info$maxs[i])  # avoid things like 3:3
                node[[i+2]] <- info$mins[i]
        }
    }
    if(!returnExpr)
        node <- deparse(node)
    return(node)    
}

findLatentNodes <- function(model, nodes, timeIndex = NULL) {
    ## Determine set of latent 'nodes', one per time point.
    ## Note that each time point might have one node or a set of nodes.
    varName <- sapply(nodes, function(x) {
        return(model$getVarNames(nodes = x))
    })
    if(length(unique(varName)) > 1){
        stop("findLatentNodes: all latent nodes must come from same variable.")
    }
    varName <- varName[1]
    info <- model$getVarInfo(varName)

    if(length(nodes) > 1) {
        ## Check for and fill in empty dimensions if more than 1 dimension.
        ## Otherwise, assume user-provided indexing is valid.
        nodes <- sapply(nodes, fillIndices, info, returnExpr = FALSE,
                        USE.NAMES = FALSE)
    } else {
        if(nodes == varName)
            ## 'nodes' is a variable, so setup indexing
            nodes <- paste0(varName, '[', paste0(rep(',', info$nDim - 1), collapse = ''), ']')
        
        nodeExpr <- fillIndices(nodes, info, returnExpr = TRUE)
        indexLengths <- sapply(nodeExpr[3:length(nodeExpr)],
                               function(x) length(eval(x)))
        ## Determine time index as longest dimension if not provided
        if(is.null(timeIndex)){
            maxLength <- max(indexLengths)
            if(sum(indexLengths == maxLength) > 1)         
                stop("findLatentNodes: unable to determine which dimension indexes time. Specify manually using the 'timeIndex' control list argument.")
            timeIndex <- which.max(indexLengths)
            timeLength <- maxLength
        } else{
            timeLength <- indexLengths[timeIndex]
        }
        
        timeIndices <- nodeExpr[[timeIndex+2]]  # timeIndices is unevaluated
        
        ## Expand so will have one 'node' per time point, overwriting time dimension.
        nodes <- rep('', timeLength)
        cnt <- 1
        for(i in eval(timeIndices)) {
            nodeExpr[[2+timeIndex]] <- as.numeric(i)
            nodes[cnt] <- deparse(nodeExpr)
            cnt <- cnt + 1
        }
    }
    return(nodes)
}


