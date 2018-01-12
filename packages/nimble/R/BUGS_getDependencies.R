


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



nimDerivsInfoClass <- setRefClass(
    'nimDerivsInfoClass',
    fields = list(
      allWrtAndCalcNodeNames = 'ANY',
      wrtNodeNames = 'ANY',
      calcNodeNames = 'ANY',
      parentIndicesList = 'ANY',
      stochNodeIndicators = 'ANY',
      calcNodeIndicators = 'ANY',
      cppWrtArgIndices = 'ANY',
      wrtNodeIndicators = 'ANY',
      wrtToIndices = 'ANY',
      wrtFromIndices = 'ANY',
      wrtLineIndices = 'ANY',
      wrtLineSize = 'ANY',
      wrtLineNums = 'ANY',
      lineWrtArgsAsCharacters = 'ANY',
      lineWrtArgSizeInfo = 'ANY',
      calcWithArgsCalls = 'ANY',
      nodeLengths = 'ANY',
      model =  'ANY',
      declIDs = 'ANY',
      rowIndices = 'ANY',
      topLevelWrtDeps = 'ANY',
      allNeededWRTCopyVars = 'ANY',
      isAddedScalarNode = 'ANY',
      thisAddedNodeJacobianList = 'ANY'
    ),
    methods = list(
      initialize = function(wrtNodes = NA, calcNodes = NA, thisModel = NA, cInfo = FALSE, ...){
        model <<- thisModel
        calcNodeNames <<- model$expandNodeNames(calcNodes)
        wrtNodeNames <<- model$expandNodeNames(wrtNodes, returnScalarComponents = TRUE)
        callSuper(...)
        ## This function takes a set of dependencies and returns a list with the original dependencies and a
        ## set of enhanced information needed for chain-ruling derivatives
        ##
        ## allWrtAndCalcNodeNames is a vector of nodes returned by model$getDependencies(wrtNodes)
        ## inputNodes should also be in allwrtAndCalcNodeNames
        ##
        ## convert inputNodes and deps from character to integer IDs
        allWrtAndCalcNodeNames <<- model$expandNodeNames(c(wrtNodeNames, calcNodeNames), sort = TRUE)
        nfv <- nodeFunctionVector(model, allWrtAndCalcNodeNames, sortUnique = FALSE)
        
        wrtNodes <- model$modelDef$nodeName2GraphIDs(model$expandNodeNames(wrtNodeNames))
        depIDs <- model$modelDef$nodeName2GraphIDs(allWrtAndCalcNodeNames)
        calcNodes <-  model$modelDef$nodeName2GraphIDs(model$expandNodeNames(calcNodeNames))
        maps <- model$modelDef$maps
 
        ## we want to check for children that are scalar representations of multivariate nodes
        depIdIndex <- 1
        isAddedScalarNode <<- c()
        thisAddedNodeJacobianList <<- list()
        addedScalarNodeParentInds <- c()
        for(i in seq_along(depIDs)) {
          thisNode <- depIDs[depIdIndex]
          ## Follow its descendents that are also in deps
          ## toNodes will be the children of thisNode
          toNodes <- maps$edgesFrom2To[[ thisNode ]]
          ## parentExprIDs will be the argument ID that thisNode represents to each of its child nodes
          # parentExprIDs <- maps$edgesFrom2ParentExprID[[ thisNode ]]
        addDepIds <- c()
        for(jNode in toNodes){
          if(maps$graphID_2_nodeName[jNode] %in% model$expandNodeNames(allWrtAndCalcNodeNames[depIdIndex], returnScalarComponents = TRUE)){
            if((length(maps$edgesFrom2To[[jNode]]) > 0) && !(maps$edgesFrom2To[[jNode]] %in% addDepIds)){
              addDepIds <- c(addDepIds, jNode)
              thisAddedNodeJacobianList[[length(thisAddedNodeJacobianList) + 1]] <<- matrix(0, nrow = 1, ncol = length(model[[maps$graphID_2_nodeName[thisNode]]]))
              thisAddedNodeJacobianList[[length(thisAddedNodeJacobianList)]][1, which(maps$graphID_2_nodeName[jNode] == model$expandNodeNames(allWrtAndCalcNodeNames[depIdIndex], returnScalarComponents = TRUE))] <<- 1
            }
          }
        }
        isAddedScalarNode <<- c(isAddedScalarNode, 0)
        if(length(addDepIds) > 0){
          isAddedScalarNode <<- c(isAddedScalarNode, rep(1, length(addDepIds)))
        }
        newDepIDs <- c(depIDs[1:depIdIndex], addDepIds)
        if(depIdIndex + 1 <= length(depIDs)){
          depIDs <- c(newDepIDs, depIDs[(depIdIndex + 1):length(depIDs)])
        }
        depIdIndex <- depIdIndex + 1 + length(addDepIds)
        
        }
        allWrtAndCalcNodeNames <<- maps$graphID_2_nodeName[depIDs]

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
        declIDs <<- indexingInfo$declIDs
        rowIndices <<-  indexingInfo$unrolledIndicesMatrixRows
        numNodes <- length(declIDs)

        declIDlengths <- sapply(1:numNodes, function(x){
          length(model$expandNodeNames(
            lapply(model$modelDef$declInfo[[declIDs[x]]]$symbolicParentNodesReplaced, function(y){
              if(!rowIndices[x] == 0){
                deparse(recurseReplaceIndices(y,
                                              model$modelDef$declInfo[[declIDs[x]]]$unrolledIndicesMatrix[rowIndices[x],]))
              }
              else{
                deparse(y)
              }
            }))) + 1
        })
        
        ## add length 1 addedScalarNodes to declIDlengths
        for(i in 1:length(depIDs)){
          if(isAddedScalarNode[i]){
            addDeclIDlengths <- c(declIDlengths[1:(i-1)], 2)
            if(i <= length(depIDs)) addDeclIDlengths <- c(addDeclIDlengths, declIDlengths[(i):length(declIDlengths)])
            declIDlengths <- addDeclIDlengths
          }
        }
        

        depIndex_2_parentDepIndices <- lapply(declIDlengths, function(x){
          outList <- list()
          for(i in 1:x){
            outList[[i]] <- 0
          }
          return(outList)}
        )

        stochNodeIndicators <<- model$getNodeType(allWrtAndCalcNodeNames) == 'stoch'
        wrtNodeIndicators <<- numeric(length(allWrtAndCalcNodeNames))
        calcNodeIndicators <<- numeric(length(allWrtAndCalcNodeNames))
        nodeLengths <<- numeric(length(allWrtAndCalcNodeNames))
        topLevelWrtDeps <<- depIndex_2_parentDepIndices
        allNeededWRTCopyVars <<- list()
        for(i in seq_along(topLevelWrtDeps))
          for(j in seq_along(topLevelWrtDeps[[i]]))
              topLevelWrtDeps[[i]][[j]] <<- list(0)
              
        ## For each input depsID
        for(i in seq_along(depIDs)) {
          if(isAddedScalarNode[i]){
            nodeLengths[i] <<- 1
            calcNodeIndicators[i] <<- 0
            depIndex_2_parentDepIndices[[i]][[1]] <- 0
          }
          else{
          nodeLengths[i] <<- length(values(model, allWrtAndCalcNodeNames[i]))
          thisNode <- depIDs[i]
          if(thisNode %in% calcNodes){
            calcNodeIndicators[i] <<- 1
          }
          else{
            calcNodeIndicators[i] <<- 0
          }
          if(thisNode %in% wrtNodes && (stochNodeIndicators[i] || (!calcNodeIndicators[i]))) {
            depIndex_2_parentDepIndices[[i]][[1]] <- -which(wrtNodes == thisNode) ## e.g. set -2 for 2nd wrt node
            topLevelWrtDeps[[i]][[1]][[1]] <<- which(wrtNodes == thisNode)
            wrtNodeIndicators[i] <<- 1
          }
          else{
            depIndex_2_parentDepIndices[[i]][[1]] <- 0
            # topLevelWrtDeps[[i]][[1]] <<- 0
          }
          }
        }
        
        ### next let's do wrt info
        scalarWrtNames <- model$expandNodeNames(wrtNodeNames, returnScalarComponents = TRUE)
        wrtToIndices <<- list()
        wrtFromIndices <<- list()
        wrtLineIndices <<- list()
        wrtLineSize <<- list()
        wrtLineNums <<- numeric(length(model$expandNodeNames(wrtNodeNames)))
        thisIndex <- 1
        
        for(i in seq_along(model$expandNodeNames(wrtNodeNames))){
          ## length of the node
          wrtLineSize[[i]] <<- length(model[[model$expandNodeNames(wrtNodeNames)[i]]])
          wrtLineNums[i] <<- which(model$expandNodeNames(wrtNodeNames)[i] == allWrtAndCalcNodeNames)
          ## function below, for each scalar element of the i'th wrt par, returns the index
          ## of that element in the vector of all wrt pars.
          thisWrtNodeInds <- sapply(model$expandNodeNames(model$expandNodeNames(wrtNodeNames)[i],
                                                          returnScalarComponents = TRUE),
                                    function(x){
                                      outInd <- which(x == scalarWrtNames)
                                      if(length(outInd) > 0){
                                        return(outInd)
                                      }
                                      else{return(0)}
                                    }
          )
          
          ## toIndices are the indices of the returned deriv element (e.g. gradient)
          ## that this wrt param will map to
          wrtToIndices[[i]] <<- thisWrtNodeInds[which(thisWrtNodeInds != 0)]
          ## lineIndices are the full indices of this wrt node (could be longer than toIndices if only one element of a multivar node is used for wrt)
          wrtLineIndices[[i]] <<- thisIndex:(thisIndex + wrtLineSize[[i]] - 1)
          ## fromIndices are the elements of the calculated derivative that will be placed into the toIndices of the output deriv.
          wrtFromIndices[[i]] <<- wrtLineIndices[[i]][which(thisWrtNodeInds != 0)]
          thisIndex <- thisIndex + wrtLineSize[[i]]
        }
        
        for(i in seq_along(depIDs)) {
          thisNode <- depIDs[i]
          ## Follow its descendents that are also in deps
          ## toNodes will be the children of thisNode
          toNodes <- maps$edgesFrom2To[[ thisNode ]]
          ## parentExprIDs will be the argument ID that thisNode represents to each of its child nodes
          parentExprIDs <- maps$edgesFrom2ParentExprID[[ thisNode ]]
          parentExprIDs[is.na(parentExprIDs)] <- 1
          # addToNodes <- c()
          # for(jNode in toNodes){
          #   if(maps$graphID_2_nodeName[jNode] %in% model$expandNodeNames(allWrtAndCalcNodeNames[i], returnScalarComponents = TRUE)){
          #     if((length(maps$edgesFrom2To[[jNode]]) > 0) && !(maps$edgesFrom2To[[jNode]] %in% addToNodes)){
          #       addToNodes <- c(addToNodes, maps$edgesFrom2To[[jNode]])
          #       parentExprIDs <- c(parentExprIDs, maps$edgesFrom2ParentExprID[[jNode]])
          #     }
          #   }
          # }
          # toNodes <- c(toNodes, addToNodes)
          ## for each child node
          for(iTo in seq_along(toNodes)) {
            thisToNode <- toNodes[iTo]
            ## Check if this child is in depIDs
            thisParentExprID <- parentExprIDs[iTo]
            
            if(thisToNode %in% depIDs) {
              ## Populate an entry in the results
              iThisToNode <- which(depIDs == thisToNode)
              if(wrtNodeIndicators[iThisToNode] == 1 ||
                 ((calcNodeIndicators[iThisToNode] == 1 || isAddedScalarNode[iThisToNode] == 1) &&
                  (wrtNodeIndicators[i] == 1 || stochNodeIndicators[i] == 0))){
                if(!is.na(thisParentExprID)){
                  if(length(depIndex_2_parentDepIndices[[iThisToNode]][[ thisParentExprID + 1 ]]) == 1 &&
                     depIndex_2_parentDepIndices[[iThisToNode]][[ thisParentExprID + 1 ]][1] == 0){
                    depIndex_2_parentDepIndices[[iThisToNode]][[ thisParentExprID + 1 ]] <- i
                  }
                  else{
                    depIndex_2_parentDepIndices[[iThisToNode]][[ thisParentExprID + 1 ]] <- c(depIndex_2_parentDepIndices[[iThisToNode]][[ thisParentExprID + 1 ]], i)
                  }
                }
              }
              if(wrtNodeIndicators[i] && stochNodeIndicators[i]){
                if(length(topLevelWrtDeps[[iThisToNode]][[ thisParentExprID + 1 ]]) == 1 && (topLevelWrtDeps[[iThisToNode]][[ thisParentExprID + 1 ]][[1]][1] == 0)){
                  topLevelWrtDeps[[iThisToNode]][[ thisParentExprID + 1 ]][[1]] <<- which(i == wrtLineNums)
                }
                else{
                  topLevelWrtDeps[[iThisToNode]][[ thisParentExprID + 1 ]][[length(topLevelWrtDeps[[iThisToNode]][[ thisParentExprID + 1 ]]) + 1]] <<- which(i == wrtLineNums)
                }
              } else if(stochNodeIndicators[i] == 0){
                if(length(topLevelWrtDeps[[iThisToNode]][[ thisParentExprID + 1 ]]) == 1 && (!is.na(topLevelWrtDeps[[iThisToNode]][[ thisParentExprID + 1 ]][[1]][1]) &&
                                                                                             topLevelWrtDeps[[iThisToNode]][[ thisParentExprID + 1 ]][[1]][1] == 0)){
                  if(!all(unlist(topLevelWrtDeps[[i]]) == 0)){
                    topLevelWrtDeps[[iThisToNode]][[ thisParentExprID + 1 ]][[1]] <<- setdiff(unique(unlist(topLevelWrtDeps[[i]])), c(0, NA))
                  }
                  else{
                    topLevelWrtDeps[[iThisToNode]][[ thisParentExprID + 1 ]][[1]] <<- NA
                  }
                }
                else{
                  if(!all(unlist(topLevelWrtDeps[[i]]) == 0)){
                    topLevelWrtDeps[[iThisToNode]][[ thisParentExprID + 1 ]][[length(topLevelWrtDeps[[iThisToNode]][[ thisParentExprID + 1 ]]) + 1]] <<-
                      setdiff(unique(unlist(topLevelWrtDeps[[i]])),c(0, NA))
                  }
                  else{
                    topLevelWrtDeps[[iThisToNode]][[ thisParentExprID + 1 ]][[length(topLevelWrtDeps[[iThisToNode]][[ thisParentExprID + 1 ]]) + 1]] <<-
                     NA
                  }
                }
              }
            }
          }
          allNeededWRTCopyVars[[i]] <<- sort(setdiff(unique(unlist(topLevelWrtDeps[[i]])), c(0, NA)))
        }
        for(i in seq_along(topLevelWrtDeps)){
          for(j in seq_along(topLevelWrtDeps[[i]])){
            for(k in seq_along(topLevelWrtDeps[[i]][[j]])){
              topLevelWrtDeps[[i]][[j]][[k]][is.na(topLevelWrtDeps[[i]][[j]][[k]])] <<- 0
            }
          }
        }
          

        ## Next lets get line wrt info.  Do both as characters (for R) and as indices (for C++)! with a flag to control which is returned.
        ## Then update the explainDerivContent function.  May also make more sense to make this function a ref class w/ fields for different return vals and such.
        lineWrtArgsAsCharacters <<- list()
        lineWrtArgSizeInfo <<- list()
        calcWithArgsCalls  <<- list()
        cppWrtArgIndices <<- list()
        tempDeclIDs <- declIDs
        for(i in seq_along(depIndex_2_parentDepIndices)){
          if(calcNodeIndicators[i] == 1){
            lineWrtArgSizeInfo[[i]]      <<- numeric(length(depIndex_2_parentDepIndices[[i]]))
            sizeAndDimInfo <- environment(model$nodeFunctions[[tempDeclIDs[1]]]$.generatorFunction)[['parentsSizeAndDims']]
            formalArgNames <- formals(model$nodeFunctions[[tempDeclIDs[1]]]$calculateWithArgs)
            unrolledIndicesMatrixRow <- model$modelDef$declInfo[[tempDeclIDs[1]]]$unrolledIndicesMatrix[ rowIndices[i], ]
            tempDeclIDs <- tempDeclIDs[-1]
            modelArgNames <- lapply(names(formalArgNames)[-1],
                                    function(x){parse(text = convertCalcArgNameToModelNodeName(x, sizeAndDimInfo, unrolledIndicesMatrixRow))[[1]]})
            calcWithArgsCalls[[i]] <<- as.call(c(list(as.name('calcWithArgs'), unrolledIndicesMatrixRow), modelArgNames))
            for(j in seq_along(depIndex_2_parentDepIndices[[i]])){
              if(j == 1 && wrtNodeIndicators[[i]] == 1){
                thisWrtLine <- which(wrtLineNums == i)
                lineWrtArgsAsCharacters[[i]]  <<- names(formalArgNames[2])
                lineWrtArgSizeInfo[[i]][1]  <<- wrtLineSize[[thisWrtLine]]
              }
              else if(depIndex_2_parentDepIndices[[i]][[j]][1] > 0){
                for(k in 1:length(depIndex_2_parentDepIndices[[i]][[j]])){
                  wrtInfoList <-  convertToWrtArg(allWrtAndCalcNodeNames[depIndex_2_parentDepIndices[[i]][[j]][k]],
                                                  modelArgNames[[j]],
                                                  names(formalArgNames)[j+1],
                                                  thisModel)
                  if(length(lineWrtArgsAsCharacters) < i){
                    lineWrtArgsAsCharacters[[i]] <<- wrtInfoList$wrtArg
                  }
                  else{
                    lineWrtArgsAsCharacters[[i]] <<- c(lineWrtArgsAsCharacters[[i]], wrtInfoList$wrtArg)
                  }
                  # lineWrtArgSizeInfo[[i]][j] <<- lineWrtArgSizeInfo[[i]][j] + wrtInfoList$argSize
                }
                lineWrtArgSizeInfo[[i]][j] <<- length(eval(modelArgNames[[j]]))
              }
            }
            if(length(lineWrtArgsAsCharacters) < i){
              lineWrtArgsAsCharacters[[i]] <<- NA
            }
            if(cInfo){
              functionArgsDimsList <- formalArgNames[-1]
              argCounter <- 1
              for(j in seq_along(sizeAndDimInfo)){
                for(k in seq_along(sizeAndDimInfo[[j]])){
                  functionArgsDimsList[[argCounter]] <- substitute(double(PARDIM, PARSIZES), 
                                                                   list(PARDIM = as.numeric(sizeAndDimInfo[[j]][[k]]$nDim), 
                                                                        PARSIZES = nndf_makeParentSizeExpr(sizeAndDimInfo[[j]][[k]]))) 
                  argCounter <- argCounter + 1
                }
              }
              functionArgsDimsList <- lapply(functionArgsDimsList, function(x){
                if(deparse(x) == ""){
                  return(parse(text = "double(0, 1)")[[1]])
                }
                else{
                  return(x)
                }
              })
              cppWrtArgIndices[[i]] <<- convertWrtArgToIndices(lineWrtArgsAsCharacters[[i]], functionArgsDimsList, fxnName = 'calculate')
            }
          }
          else if(isAddedScalarNode[i]){
            lineWrtArgsAsCharacters[[i]] <<- NA
            lineWrtArgSizeInfo[[i]]      <<- c(0, length(model[[allWrtAndCalcNodeNames[depIndex_2_parentDepIndices[[i]][[2]]]]]))
            calcWithArgsCalls[[i]] <<-   NA
            cppWrtArgIndices[[i]] <<- -1
          }
          else{
            tempDeclIDs <- tempDeclIDs[-1]
            lineWrtArgsAsCharacters[[i]] <<- NA
            lineWrtArgSizeInfo[[i]]      <<- -1
            calcWithArgsCalls[[i]] <<-   NA
            cppWrtArgIndices[[i]] <<- -1
          }
          if(length(lineWrtArgsAsCharacters) < i){
            lineWrtArgsAsCharacters[[i]] <<- NA
          }
        }
        if(cInfo){
          for(i in seq_along(depIndex_2_parentDepIndices)){
            for(j in seq_along(depIndex_2_parentDepIndices[[i]])){
              depIndex_2_parentDepIndices[[i]][[j]] <- depIndex_2_parentDepIndices[[i]][[j]] - 1
            }
          }
        }
        parentIndicesList <<- depIndex_2_parentDepIndices
      },

      # ### A function that substitutes correct values of unrolledIndicesMatrix
      # ### into symbolicParentNodesReplaced.
      recurseReplaceIndices = function(code, unrolledIndicesRow){
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
      },
      # 
      ## A function that converts the name of an argument to a calculateWithArgs
      ## function to a character string representing that argument in the model,
      ## a possible example:
      ## calcArgName = 'y_1'
      ## function output = model$y[1:2]
      convertCalcArgNameToModelNodeName = function(calcArgName, sizeAndDimInfo, unrolledIndicesMatrixRow){
        thisModelElementNum <- as.numeric(gsub(".*([0-9]+)$", "\\1", calcArgName)) ## Extract 1, 2, etc. from end of arg name.
        thisName <- sub("_[0-9]+$","",calcArgName) ## Extract node name from beginning of arg name.
        indexBracketInfo <- paste0('[',
                                   paste0(sapply(sizeAndDimInfo[[thisName]][[thisModelElementNum]]$indexExpr, function(x){
                                     if(length(x) == 1) return(deparse(x[[1]]))
                                     else if( deparse(x[[1]]) == 'getNodeFunctionIndexedInfo'){
                                       rowNum <- x[[3]]
                                       x <- eval(parse(text = paste0('unrolledIndicesMatrixRow[',rowNum,']'))[[1]])
                                     }
                                     else{
                                       return(paste0(deparse(x[[1]]), ':', deparse(x[[2]])))
                                     }}), collapse = ', '),
                                   ']')
        return(paste0('model$', thisName, indexBracketInfo))
      },
      # A function that takes a wrt name supplied by a user and returns a
      # wrt name that can be used in a call to nimDerivs(calcWithArgs()).
      # function output:  a list with elements:
      #    wrtArg:  a character string denoting the actual wrt argument
      #             that will be given to the call to nimDerivs(calcWithArgs())
      #    argSize: the flattened length of that argument
      convertToWrtArg = function(wrtName, modelArgName, fxnArgName, thisModel){
        modelName <- strsplit(deparse(modelArgName), "\\$")[[1]][2]
        overlapPars <- sapply(thisModel$expandNodeNames(wrtName, returnScalarComponents = TRUE), function(x){which(x == thisModel$expandNodeNames(modelName, returnScalarComponents = TRUE))})
        return(list(wrtArg = paste0(fxnArgName, '[c(', paste(overlapPars, collapse = ', ') ,')]'),
                    argSize = length(values(thisModel, wrtName))))
      },
      ## Explains some of the elements of the nimDerivsInfoClass
      explainDerivContent = function(enhancedDeps) {
        writeLines('The following calculations would be done from this input:')
        deps <- allWrtAndCalcNodeNames
        depIDs <- model$modelDef$nodeName2GraphIDs(deps)
        currentDeclIDs <- model$modelDef$maps$graphID_2_declID[depIDs]
        for(i in seq_along(depIDs)) {
          description <- c()
          derivInfo <- parentIndicesList
          thisDerivInfo <- derivInfo[[i]]
          argumentNames <- lapply( model$modelDef$declInfo[[ currentDeclIDs[i] ]]$symbolicParentNodesReplaced, deparse)
          for(j in seq_along(thisDerivInfo)){
            if(j == 1 && thisDerivInfo[[j]][1] < 0){
              description <- paste0('wrt parameter ', -thisDerivInfo[[1]][1])
            }
            else{
              if(thisDerivInfo[[j]][1] > 0){
                description <- c(description, 
                                 paste0('Argument ', j, ' (', argumentNames[j - 1],') comes from calculation(s) ', thisDerivInfo[[j]], collapse = ' '))
              }
            }
          }
          if(length(description) == 0){
            description <- "(No arguments are wrt or deterministic arguments from previous calculations)\n"
          }
          BUGSline <- deparse(model$modelDef$declInfo[[ currentDeclIDs[i] ]]$codeReplaced)
          output <- paste0(i,': ', deps[i], ' (from ', BUGSline, ')\n', paste0('\t', description, collapse = '\n'))
          writeLines(output)
        }
      }
    )
)
