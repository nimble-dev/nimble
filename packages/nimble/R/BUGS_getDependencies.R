


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
      model =  'ANY'
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

        stochNodeIndicators <<- model$getNodeType(allWrtAndCalcNodeNames) == 'stoch'
        wrtNodeIndicators <<- numeric(length(allWrtAndCalcNodeNames))
        calcNodeIndicators <<- numeric(length(allWrtAndCalcNodeNames))
        nodeLengths <<- numeric(length(allWrtAndCalcNodeNames))
        ## For each input depsID
        for(i in seq_along(depIDs)) {
          nodeLengths[i] <<- length(values(model, allWrtAndCalcNodeNames[i]))
          thisNode <- depIDs[i]
          if(thisNode %in% wrtNodes) {
            depIndex_2_parentDepIndices[[i]][[1]] <- -which(wrtNodes == thisNode) ## e.g. set -2 for 2nd wrt node
            wrtNodeIndicators[i] <<- 1
          }
          else{
            depIndex_2_parentDepIndices[[i]][[1]] <- 0
          }
          if(thisNode %in% calcNodes){
            calcNodeIndicators[i] <<- 1
          }
          else{
            calcNodeIndicators[i] <<- 0
          }
        }
        for(i in seq_along(depIDs)) {
          thisNode <- depIDs[i]
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
              ## Populate an entry in the results
              iThisToNode <- which(depIDs == thisToNode)
              if(wrtNodeIndicators[iThisToNode] == 1 ||
                 (calcNodeIndicators[iThisToNode] == 1 &&
                  (wrtNodeIndicators[i] == 1 || stochNodeIndicators[i] == 0))){
                thisParentExprID <- parentExprIDs[iTo]
                if(length(depIndex_2_parentDepIndices[[iThisToNode]][[ thisParentExprID + 1 ]]) == 1 &&
                   depIndex_2_parentDepIndices[[iThisToNode]][[ thisParentExprID + 1 ]][1] == 0){
                  depIndex_2_parentDepIndices[[iThisToNode]][[ thisParentExprID + 1 ]] <- i
                }
                else{
                  depIndex_2_parentDepIndices[[iThisToNode]][[ thisParentExprID + 1 ]] <- c(depIndex_2_parentDepIndices[[iThisToNode]][[ thisParentExprID + 1 ]], i)
                }
              }
            }
          }
        }

        parentIndicesList <<- depIndex_2_parentDepIndices

        ### next let's do wrt info
        scalarWrtNames <- model$expandNodeNames(wrtNodeNames, returnScalarComponents = TRUE)
        wrtToIndices <<- list()
        wrtFromIndices <<- list()
        wrtLineIndices <<- list()
        wrtLineSize <<- list()
        wrtLineNums <<- numeric(length(model$expandNodeNames(wrtNodeNames)))
        thisIndex <- 1

        for(i in seq_along(model$expandNodeNames(wrtNodeNames))){
          ## which node is it? let's say unnecessary for now
          ##wrtLineInfo[[i]]$lineNum <- which(model$expandNodeNames(wrtNodeNames)[i] == derivInfo[[1]])
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


        ## Next lets get line wrt info.  Do both as characters (for R) and as indices (for C++)! with a flag to control which is returned.
        ## Then update the explainDerivContent function.  May also make more sense to make this function a ref class w/ fields for different return vals and such.
        lineWrtArgsAsCharacters <<- list()
        lineWrtArgSizeInfo <<- list()
        calcWithArgsCalls  <<- list()
        cppWrtArgIndices <<- list()
        for(i in seq_along(depIndex_2_parentDepIndices)){
          if(calcNodeIndicators[i] == 1){
            lineWrtArgSizeInfo[[i]]      <<- numeric(length(depIndex_2_parentDepIndices[[i]]))
            sizeAndDimInfo <- environment(model$nodeFunctions[[declIDs[i]]]$.generatorFunction)[['parentsSizeAndDims']]
            formalArgNames <- formals(model$nodeFunctions[[declIDs[i]]]$calculateWithArgs)
            unrolledIndicesMatrixRow <- model$modelDef$declInfo[[declIDs[i]]]$unrolledIndicesMatrix[ unrolledIndicesMatrixRows[i], ]
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
                  lineWrtArgSizeInfo[[i]][j] <<- lineWrtArgSizeInfo[[i]][j] + wrtInfoList$argSize
                }
              }
            }
            if(cInfo){
              functionArgsDimsList <- formalArgNames[-1]
              argCounter <- 1
              for(j in seq_along(sizeAndDimInfo)){
                for(k in seq_along(sizeAndDimInfo[[j]])){
                  functionArgsDimsList[[argCounter]] <- substitute(double(PARDIM, PARSIZES), 
                                                                   list(PARDIM = as.numeric(sizeAndDimInfo[[j]][[k]]$nDim), 
                                                                        PARSIZES = nndf_makeParentSizeExpr(sizeAndDimInfo[[j]][[k]]))) 
                }
              }
              cppWrtArgIndices[[i]] <<- convertWrtArgToIndices(lineWrtArgsAsCharacters[[i]], functionArgsDimsList, fxnName = 'calculate')
            }
          }
          else{
            lineWrtArgsAsCharacters[[i]] <<- NA
            lineWrtArgSizeInfo[[i]]      <<- NA
            calcWithArgsCalls[[i]] <<-   NA
          }
          if(length(lineWrtArgsAsCharacters) < i){
            lineWrtArgsAsCharacters[[i]] <<- NA
          }
        }
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
        declIDs <- model$modelDef$maps$graphID_2_declID[depIDs]
        for(i in seq_along(depIDs)) {
          description <- c()
          derivInfo <- parentIndicesList
          thisDerivInfo <- derivInfo[[i]]
          argumentNames <- lapply( model$modelDef$declInfo[[ declIDs[i] ]]$symbolicParentNodesReplaced, deparse)
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
          BUGSline <- deparse(model$modelDef$declInfo[[ declIDs[i] ]]$codeReplaced)
          output <- paste0(i,': ', deps[i], ' (from ', BUGSline, ')\n', paste0('\t', description, collapse = '\n'))
          writeLines(output)
        }
      }
    )
)
