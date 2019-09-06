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
