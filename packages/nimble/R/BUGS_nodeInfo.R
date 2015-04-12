## The nodeInfoClass contains the fully-indexed information of one BUGS declaration (node)
## The indexedNodeInfo field of BUGSdeclClass is a list of nodeInfoClass objects

#' nodeInfoClass contains information for a node
nodeInfoClass <- setRefClass('nodeInfoClass',
                             fields = list(
                                 code = 'ANY',
                                 type = 'ANY',
                                 targetNodeExpr = 'ANY',
                                 targetNodeName = 'ANY',
                                 targetNodeIndexValuesList = 'ANY',
                                 targetNodeIndexSizes = 'ANY',
                                 logProbIndexValues = 'ANY',
                                 targetVarName = 'ANY',
                                 parentNodeExprs = 'ANY',
                                 parentNodeNames = 'ANY',
                                 nodeFunctionName = 'ANY',
                                 indexVariableExprs = 'ANY',
                                 indexVariableValues = 'ANY',
                                 codeReplaced = 'ANY',
                                 replacementValues = 'ANY',
                                 codeReplacedWithValues = 'ANY',
                                 logProbNodeReplacedWithValues = 'ANY',
                                 altParamExprs = 'ANY',
                                 altParamExprsWithValues = 'ANY'
                             ),
                             methods = list(
                                 getDistribution = function() {
                                     if(type != 'stoch')  stop('getting distribution of non-stochastic node')
                                     return(as.character(codeReplacedWithValues[[3]][[1]]))
                                 },
                                 getValueExpr = function() {
                                     return(codeReplacedWithValues[[3]])
                                 },
                                 getParamExpr = function(param) {
                                     if(type != 'stoch')  stop('getting parameter expression for a non-stochastic node')
                                     if(param %in% names(codeReplacedWithValues[[3]]))     return(codeReplacedWithValues[[3]][[param]])
                                     if(param %in% names(altParamExprsWithValues))         return(altParamExprsWithValues[[param]])
                                     stop('getting a parameter not present in stochastic node')
                                 },
                                 isDiscrete = function() {
                                     if(type != 'stoch')  stop('querying whether a non-stochastic node is \'discrete\'')
                                     dist <- getDistribution()
                                     return(getDistribution(dist)$discrete)
                                 }
                             ))
