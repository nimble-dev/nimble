


## returns a list of all deterministic dependents up to the first stochastic dependent,
## omitting any nodes in 'omit', or down the path of 'omit' nodes
## works only in terms of vertex IDs, as in the igraph object.
gd_getDependencies_IDs <- function(graph, maps, nodes, omit, downstream) {
 #   nonStochNodes <- which(maps$types != 'stoch') 		We should be able to speed things up by looking up by graphID instead of intersecting...
    nodes <- setdiff(nodes, omit)
    newNodes <- if(length(nodes) > 0)    unlist(maps$edgesFrom2To[nodes]) else integer(0) 
    newNodes <- setdiff(newNodes, omit)
    while(length(newNodes) > 0) {
        nodes <- c(nodes, newNodes)
        newNodesForRecursion <- if(downstream)   newNodes   else  newNodes[maps$types[newNodes] != 'stoch'] 
        newNodes <- if(length(newNodesForRecursion) > 0)  unlist(maps$edgesFrom2To[newNodesForRecursion]) else integer(0) 
        newNodes <- setdiff(newNodes, omit)
    }
    nodes <- unique(nodes)
    nodes <- sort(nodes)    # topological sort
    return(nodes)
}


gd_allNeighbors <- function(graph, nodes) stop("shouldn't be calling gd_allNeighbors any more")

## gd_getDependencies_IDs <- function(graph, maps, nodes, omit, downstream) {
##  #   nonStochNodes <- which(maps$types != 'stoch') 		We should be able to speed things up by looking up by graphID instead of intersecting...
##     nodes <- setdiff(nodes, omit)
##     newNodes <- if(length(nodes) > 0)     gd_allNeighbors(graph, nodes)     else     numeric(0)
##     newNodes <- setdiff(newNodes, omit)
##     while(length(newNodes) > 0) {
##         nodes <- c(nodes, newNodes)
##         newNodesForRecursion <- if(downstream)   newNodes   else  newNodes[maps$types[newNodes] != 'stoch'] 			#intersect(newNodes, nonStochNodes)
##         newNodes <- if(length(newNodesForRecursion) > 0)     gd_allNeighbors(graph, newNodesForRecursion)     else     numeric(0)
##         newNodes <- setdiff(newNodes, omit)
##     }
##     nodes <- unique(nodes)
##     nodes <- sort(nodes)    # topological sort
##     return(nodes)
## }

## gd_allNeighbors <- function(graph, nodes) {
##     nei <- numeric(0)
##     for(n in nodes)  nei <- c(nei, neighbors(graph, n))
##     return(nei)
## }

