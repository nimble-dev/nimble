


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


