
## CASE 2: modelVariableAccessorVector
## copy(model1, xxx, nodes) becomes:
## model1_nodes_accessors <- modelVariableAccessorVector(model1, nodes)
## copy(model1_nodes_accessors, xxx)
## model1_nodes_accessors$getAccessors() returns a list of the modelVariableAccessor objects


modelVariableAccessorVector <- function(mv, nodeNames, logProb = FALSE, logProbOnly = FALSE) {
    ans <- list(mv, substitute(nodeNames), logProb, parent.frame(), logProbOnly)
    class(ans) <- c('modelVariableAccessorVector', 'valuesAccessorVector')
    ans
}

