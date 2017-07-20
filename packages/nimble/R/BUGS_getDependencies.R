


## returns a list of all deterministic dependents up to the first stochastic dependent,
## omitting any nodes in 'omit', or down the path of 'omit' nodes
## works only in terms of vertex IDs, as in the igraph object.
gd_getDependencies_IDs <- function(graph, maps, nodes, omit, downstream) {
  # nonStochNodes <- which(maps$types != 'stoch') 		We should be able to speed things up by looking up by graphID instead of intersecting...

    nodes <- setdiff(nodes, omit)
    newNodes <- if(length(nodes) > 0)    unlist(maps$edgesFrom2To[nodes]) else integer(0)  ## first set of dependencies, including from LHSinferred
    newNodes <- setdiff(newNodes, omit)
    
    boolLHSinferred <- maps$types[nodes] == 'LHSinferred'
    LHSinferredNodes <- nodes[boolLHSinferred]
    nodes <- setdiff(nodes[!boolLHSinferred], omit) ## filter out LHSinferred
    
    if(length(LHSinferredNodes)>0) {
        ## something like x[1], an inferred piece of x[1:10] because x[1] appeared somewhere on its own
        ## include the node it is from (x[1:10]) and well as *non-inferred* dependencies of that node
        fullNodes <- unique(maps$vertexID_2_nodeID[ LHSinferredNodes ]) ## get the x[1:10]
        fullNodes <- setdiff(fullNodes, omit)  ## filter omits
        nodes <- c(nodes, fullNodes)           ## add to nodes

        fullNodesForRecursion <- fullNodes
      ##  fullNodesForRecursion <- if(downstream)   fullNodes   else  fullNodes[maps$types[fullNodes] != 'stoch'] ## find recursion nodes

        fullNodesDeps <- if(length(fullNodesForRecursion) > 0) unlist(maps$edgesFrom2To[fullNodesForRecursion]) else integer(0) ## get dependencies of x[1:10]
        
        fullNodesDepsLHSinferred <- maps$types[fullNodesDeps] == 'LHSinferred' ## filter out LHSinferred dependencies, e.g. x[2]
        fullNodesDeps <-  fullNodesDeps[!fullNodesDepsLHSinferred]
        fullNodesDeps <- setdiff(fullNodesDeps, omit) ## filter omits
        
        newNodes <- c(newNodes, fullNodesDeps)
    }
    
    while(length(newNodes) > 0) {
        nodes <- c(nodes, newNodes)
        newNodesForRecursion <- if(downstream)   newNodes   else  newNodes[maps$types[newNodes] != 'stoch'] 
        newNodes <- if(length(newNodesForRecursion) > 0)  unlist(maps$edgesFrom2To[newNodesForRecursion]) else integer(0) 
        newNodes <- setdiff(newNodes, omit)
    }
    nodes <- unique(nodes)
    nodes <- sort(nodes)    # topological sort
    return(nodes)

    ## old
    ## nodes <- setdiff(nodes, omit)
    ## newNodes <- if(length(nodes) > 0)    unlist(maps$edgesFrom2To[nodes]) else integer(0) 
    ## newNodes <- setdiff(newNodes, omit)
    ## while(length(newNodes) > 0) {
    ##     nodes <- c(nodes, newNodes)
    ##     newNodesForRecursion <- if(downstream)   newNodes   else  newNodes[maps$types[newNodes] != 'stoch'] 
    ##     newNodes <- if(length(newNodesForRecursion) > 0)  unlist(maps$edgesFrom2To[newNodesForRecursion]) else integer(0) 
    ##     newNodes <- setdiff(newNodes, omit)
    ## }
    ## nodes <- unique(nodes)
    ## nodes <- sort(nodes)    # topological sort
    ## return(nodes)

}


gd_allNeighbors <- function(graph, nodes) stop("shouldn't be calling gd_allNeighbors any more")



enhanceDepsForDerivs <- function(inputNodes, deps, model, nfv) {
  ## This function takes a set of dependencies and returns a list with the original dependencies and a
  ## set of enhanced information needed for chain-ruling derivatives
  ##
  ## deps is a vector of nodes returned by model$getDependencies(inputNodes)
  ## inputNodes should also be in the deps.  
  ## 
  ## convert inputNodes and deps from character to integer IDs
  if(!is.integer(inputNodes)) inputNodes <- model$modelDef$nodeName2GraphIDs(inputNodes)
  depIDs <- if(!is.integer(deps)) model$modelDef$nodeName2GraphIDs(deps) else deps
  
  maps <- model$modelDef$maps
  
  ## get the BUGS declaration ID for every node

  ## initialize the enhanced information
  ## Elements of depIndex_2_parentDepIndices correspond to elements of deps
  ## depIndex_2_parentDepIndices[[i]] will have one of two formats:
  ##    (1) a single negative integer.  This gives the (-) index of inputNodes corresponding to this node.
  ##     e.g. if inputNodes is c('x[1]', 'x[2]'), and these are elements 1 and 2 in deps, then
  ##      depIndex_2_parentDepIndices[[1]] will be -1
  ##      depIndex_2_parentDepIndices[[2]] will be -2
  ##    (2) a vector of integers giving the calculation index of deps corresponding to each input parameter
  ##     e.g. if deps[5] is y[3], whose first argument is beta and second argument is x[4], then
  ##       depIndex_2_parentDepIndices[[2]] will be c(0, 3)
  ##           The 0 means that beta is not part of deps
  ##           The 3 means that x[4] is deps[3]
  indexingInfo <- nfv$indexingInfo
  declIDs <- indexingInfo$declIDs
  numNodes <- length(declIDs)
  unrolledIndicesMatrixRows <- indexingInfo$unrolledIndicesMatrixRows
  
  ### A function that substitutes correct values of unrolledIndicesMatrix into symbolicParentNodesReplaced
  recurseReplaceIndices <- function(code, unrolledIndicesRow){
    replaceNames <- names(unrolledIndicesRow)
    if(length(code) > 1){
      for(i in seq_along(code)){
        if(length(code[[i]]) > 1){
          code[[i]] <- recurseReplaceIndices(code[[i]], unrolledIndicesRow)
        }
        else if(deparse(code[[i]]) %in% replaceNames){
          code[[i]] <- unrolledIndicesRow[deparse(code[[i]])]
        }
      }
    }
    else if(deparse(code) %in% replaceNames){
      code <- unrolledIndicesRow[deparse(code)]
    }
    return(code)
  }
  
  declIDlengths <- sapply(1:numNodes, function(x){
    length(model$expandNodeNames(
      lapply(model$modelDef$declInfo[[declIDs[x]]]$symbolicParentNodesReplaced, function(y){
        if(!unrolledIndicesMatrixRows[x] == 0){
          deparse(recurseReplaceIndices(y,
                                        model$modelDef$declInfo[[declIDs[x]]]$unrolledIndicesMatrix[unrolledIndicesMatrixRows[x],]))
        }
        else{
          deparse(y)
        }
        }))) + 1
  })
  

  depIndex_2_parentDepIndices <- lapply(declIDlengths, function(x){
    outList <- list()
    for(i in 1:x){
      outList[[i]] <- 0
    }
    return(outList)}
  )
  
  ## For each input depsID
  for(i in seq_along(depIDs)) {
    thisNode <- depIDs[i]
    ## Check if it is an input node
    if(thisNode %in% inputNodes) {
      depIndex_2_parentDepIndices[[i]][[1]] <- -which(inputNodes == thisNode) ## e.g. set -2 for 2nd input node
    } 
    else{
      depIndex_2_parentDepIndices[[i]][[1]] <- 0
    }
    ## Follow its descendents that are also in deps
    ## toNodes will be the children of thisNode
    toNodes <- maps$edgesFrom2To[[ thisNode ]]
    ## parentExprIDs will be the argument ID that thisNode represents to each of its child nodes
    parentExprIDs <- maps$edgesFrom2ParentExprID[[ thisNode ]]
    ## for each child node
    for(iTo in seq_along(toNodes)) {
      thisToNode <- toNodes[iTo]
      ## Check if this child is in depIDs
      if(thisToNode %in% depIDs) {
        ## Populate an entry in the resultss
        iThisNodeInDeps <- which(depIDs == thisToNode)
        thisParentExprID <- parentExprIDs[iTo]
        if(length(depIndex_2_parentDepIndices[[iThisNodeInDeps]][[ thisParentExprID + 1 ]]) == 1 && 
           depIndex_2_parentDepIndices[[iThisNodeInDeps]][[ thisParentExprID + 1 ]][1] == 0){
          depIndex_2_parentDepIndices[[iThisNodeInDeps]][[ thisParentExprID + 1 ]] <- i
        }
        else{
          depIndex_2_parentDepIndices[[iThisNodeInDeps]][[ thisParentExprID + 1 ]] <- c(depIndex_2_parentDepIndices[[iThisNodeInDeps]][[ thisParentExprID + 1 ]], i)
        }
      }
    }
  }
  list(deps, depIndex_2_parentDepIndices)
}

explainDerivContent <- function(enhancedDeps, model) {
    ## This function 
    writeLines('The following calculations would be done from this input:')
    deps <- enhancedDeps[[1]]
    depIDs <- if(!is.integer(deps)) model$modelDef$nodeName2GraphIDs(deps) else deps
    declIDs <- model$modelDef$maps$graphID_2_declID[depIDs]
    derivInfo <- enhancedDeps[[2]]
    for(i in seq_along(depIDs)) {
        thisDerivInfo <- derivInfo[[i]]
        if(length(thisDerivInfo) > 0) {
            if(all(thisDerivInfo < 0)) { ## initiating node
                if(length(thisDerivInfo) > 1) stop('Problem with thisDerivInfo')
                description <- paste0('Input parameter ', -thisDerivInfo)
            } else {
                used <- thisDerivInfo > 0
                if(sum(used) > 0) {
                    argumentNames <- lapply( model$modelDef$declInfo[[ declIDs[i] ]]$symbolicParentNodesReplaced, deparse)
                    description <- paste0('Argument ', seq_along(thisDerivInfo)[used], ' (', argumentNames[used],') comes from calculation ', thisDerivInfo[used])
                }
                else
                    description <- "(No arguments are from previous calculations)\n"
            }
        } else
            description <- ""
        BUGSline <- deparse(model$modelDef$declInfo[[ declIDs[i] ]]$codeReplaced)
        output <- paste0(i,': ', deps[i], ' (from ', BUGSline, ')\n', paste0('\t', description, collapse = '\n'))
        writeLines(output)
    }
}
