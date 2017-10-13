#include <nimble/nimbleCppAD.h>

nimSmartPtr<NIMBLE_ADCLASS>  NIM_DERIVS_CALCULATE(NodeVectorClassNew_derivs &nodes, NimArr<1, double> &derivOrders) {
  nimSmartPtr<NIMBLE_ADCLASS> ADlist = new NIMBLE_ADCLASS;
  const vector<NodeInstruction> &instructions = nodes.getInstructions();
  vector<NodeInstruction>::const_iterator iNode(instructions.begin());
  vector<NodeInstruction>::const_iterator iNodeEnd(instructions.end());
  bool hessianFlag = false;
  bool jacobianFlag = false;
  bool valueFlag = false;
  
  (*ADlist).value.setSize(1);

  for(int i = 0; i < derivOrders.dimSize(0); i++){
	  if(derivOrders[i] == 0){
		  valueFlag = true;
	  }
      if(derivOrders[i] == 1){
		  jacobianFlag = true;
	  }
	  if(derivOrders[i] == 2){
		  hessianFlag = true;
		  jacobianFlag = true;
	  }
  }
  if(valueFlag && (!jacobianFlag && !hessianFlag)){
	(*ADlist).value[0] += calculate(nodes);
    return(ADlist);
  }
  (*ADlist).gradient.setSize(1, nodes.totalOutWrtSize);
  if(hessianFlag){
	 (*ADlist).hessian.setSize(nodes.totalOutWrtSize, nodes.totalOutWrtSize, 1);
  } 
  
  bool isDeterminisitic;
  bool isWrtLine;
  bool isCalcNodeLine;
  int thisWrtLine;
  int thisNodeSize;
  int thisIndex;
  int wrtLineIndicator = 0;
  vector< NimArr<2, double> > parentJacobians;
  vector< NimArr<3, double> > parentHessians;
  vector< NimArr<2, double> > chainRuleJacobians;

  for(int i = 0; i < length(nodes.parentIndicesList); i++){
	  isDeterminisitic = nodes.stochNodeIndicators[i];
	  isWrtLine = (i == nodes.wrtLineNums[wrtLineIndicator]);
	  isCalcNodeLine = nodes.calcNodeIndicators[i];
	  thisWrtLine = nodes.cumulativeWrtLineNums[i];
	  thisNodeSize = nodes.nodeLengths[i];
	  if(isCalcNodeLine){
		parentJacobians.resize(length(nodes.parentIndicesList[i]));
        if(hessianFlag){
			parentHessians.resize(length(nodes.parentIndicesList[i]));
		}
        for(int j = 0; j < length(nodes.parentIndicesList[i]); j++){
			// we can pre-calculate sumParentDims earlier.  Do this.
			if(j == 0 && isWrtLine){
				parentJacobians[j] = NimArr<2, double>(nodes.wrtLineSize[thisWrtLine], nodes.totalWrtSize);
				for(int k = 0; k < nodes.wrtLineIndices[thisWrtLine].dimSize(0); k++){
					thisIndex = nodes.wrtLineIndices[thisWrtLine][k];
					parentJacobians[j](thisIndex, thisIndex) = 1;
				}
			}
			else if(nodes.parentIndicesList[i][j][0] > 0){
				int sumParentDims = 0;
				for(int k = 0; k < nodes.parentIndicesList[i][j].dimSize(0); k++){
					sumParentDims += chainRuleJacobians[nodes.parentIndicesList[i][j][k]].dimSize(0);
				}
				parentJacobians[j] = NimArr<2, double>(sumParentDims, nodes.totalWrtSize);
				int rowIndicator = 0;
				for(int k = 0; k < nodes.parentIndicesList[i][j].dimSize(0); k++){
					for(int l = 0; l <  chainRuleJacobians[nodes.parentIndicesList[i][j][k]].dimSize(0); l++){
						for(int m = 0; m <  chainRuleJacobians[nodes.parentIndicesList[i][j][k]].dimSize(1); m++){
							parentJacobians[j](rowIndicator, m) = chainRuleJacobians[nodes.parentIndicesList[i][j][k]](l, m);
						}
						rowIndicator++;
					}
				}
				if(hessianFlag){
					parentHessians[j] = NimArr<3, double>(nodes.totalWrtSize, nodes.totalWrtSize, sumParentDims);
					rowIndicator = 0;
					for(int k = 0; k < nodes.parentIndicesList[i][j].dimSize(0); k++){
						for(int l = 0; l <  chainRuleHessians[nodes.parentIndicesList[i][j][k]].dimSize(2); l++){
							for(int m = 0; m <  chainRuleHessians[nodes.parentIndicesList[i][j][k]].dimSize(0); m++){
								for(int n = 0; n <  chainRuleHessians[nodes.parentIndicesList[i][j][k]].dimSize(1); n++){
									parentHessians[j](m, n, rowIndicator) = chainRuleHessians[nodes.parentIndicesList[i][j][k]](m, n, l);
								}
							}
							rowIndicator++;
						}
					}
				}
			}
		}

        // if(!is.na(derivInfo$lineWrtArgsAsCharacters[[i]][1])){
          // derivList <- eval(substitute(nimDerivs(CALCCALL, DERIVORDERS, DROPARGS, WRT),
                                       // list(CALCCALL = derivInfo$calcWithArgsCalls[[i]],
                                            // DERIVORDERS = order,
                                            // DROPARGS = 'INDEXEDNODEINFO_',
                                            // WRT = derivInfo$lineWrtArgsAsCharacters[[i]])))

          // if(isDeterm){
            // derivList$value <- 0
            // model$nodeFunctions[[declID]]$calculate(unrolledIndicesMatrixRow)
          // }
          

		
	}
  }

  for(; iNode != iNodeEnd; iNode++){
    (*ADlist).value[0] += iNode->nodeFunPtr->calculateBlock(iNode->operand);
  }
  return(ADlist);  
}


// nimDerivs_calculate <- function(model, nodes = NA, order, wrtPars, silent = TRUE){
  // nfv <-  nodeFunctionVector(model, model$expandNodeNames(c(nodes, 
                                                            // model$expandNodeNames(wrtPars)), 
                                                          // sort = TRUE), 
                             // sortUnique = FALSE)
  // declIDs <- nfv$indexingInfo$declIDs
  // chainRuleDerivList <- list()
  // if(hessianFlag) chainRuleHessianList <- list()
  // ## totalOutWrtSize is the sum of the lengths of all ouput wrt parameters.
  // totalOutWrtSize <- sum(sapply(derivInfo$wrtToIndices, function(x){return(length(x))}))
  // ## outDerivList will be returned from this function.
  // outDerivList <- ADNimbleList$new()
  // if(valueFlag) outDerivList$value = 0
  // outDerivList$gradient = matrix(0, ncol = totalOutWrtSize, nrow = 1)
  // if(hessianFlag) outDerivList$hessian = array(0, dim = c(totalOutWrtSize, totalOutWrtSize, 1))
  // totalWrtSize <-  sum(sapply(derivInfo$wrtLineSize, function(x){return(x)}))
  // wrtLineNums <- derivInfo$wrtLineNums
  // calcNodesLineNums <- which(derivInfo$calcNodeIndicators == 1)
  // for(i in seq_along(derivInfo$allWrtAndCalcNodeNames)){
    // if(length(calcNodesLineNums) > 0){
      // declID <- declIDs[i]
      // isDeterm <- derivInfo$stochNodeIndicators[i] == 0
      // thisWrtLine <- which(wrtLineNums == i)
      // isWrtLine <-  derivInfo$wrtNodeIndicators[i] == 1
      // isCalcNodeLine <- i %in% calcNodesLineNums
      // ## Below shouldn't be necessary in c++ chain rule code.
      // ## If this node is a calulate node,
      // ## we need to take derivatives of its calculateWithArgs function.  The derivative function
      // ## call is evaluated below.
      // thisNodeSize <- length(values(model, derivInfo$allWrtAndCalcNodeNames[i]))  
      // unrolledIndicesMatrixRow <- model$modelDef$declInfo[[declID]]$unrolledIndicesMatrix[nfv$indexingInfo$unrolledIndicesMatrixRows[i], ]
      // if(isCalcNodeLine){
        // calcWithArgs <- model$nodeFunctions[[declID]]$calculateWithArgs
        
        // ## Below we construct two lists:
        // ## parentGradients, a list of all the gradients of the parent nodes of each argument of this node.
        // ## parentHessians, a list of all the hessians of the parent nodes of this node.  
        // parentGradients <- vector('list', length = length(derivInfo$parentIndicesList[[i]]))
        // if(hessianFlag) parentHessians <- vector('list', length = length(derivInfo$parentIndicesList[[i]]))
        // for(k in seq_along(derivInfo$parentIndicesList[[i]])){
          // if(k == 1 && isWrtLine){
            // ## The first argument (k = 1) of a node's calculateWithArgs function will always be the node itself.
            // ## If this node is a wrt node, we want to set the parent gradient of the first arg (the derivative of this node wrt itself)
            // ## to the identity.
            // parentGradients[[k]] <- matrix(0, nrow = derivInfo$wrtLineSize[[thisWrtLine]],
                                           // ncol = totalWrtSize)
            // parentGradients[[k]][, derivInfo$wrtLineIndices[[thisWrtLine]]] <- diag(derivInfo$wrtLineSize[[thisWrtLine]])
          // }
          // else if(derivInfo$parentIndicesList[[i]][[k]][1] > 0){
            // ## Otherwise, if argument k has parents that depend on a wrt node, grab the parent gradients (which will have already been calculated)
            // ## and combine them into a single matrix.
            // parentGradientsList <- chainRuleDerivList[derivInfo$parentIndicesList[[i]][[k]]]
            // parentGradients[[k]] <- parentGradientsList[[1]]
            // if(length(parentGradientsList) > 1){
              // for(i_listEntry in 2:length(parentGradientsList)){
                // parentGradients[[k]] <- rbind(parentGradients[[k]], parentGradientsList[[i_listEntry]])
              // }
            // }
            // ## Similarly, grab the parent Hessians (which will have already been calculated)
            // ## and combine them into a single array.
            // if(hessianFlag){
              // parentHessiansList <- chainRuleHessianList[derivInfo$parentIndicesList[[i]][[k]]]
              // parentHessiansDim <- sum(sapply(parentHessiansList, function(x){return(dim(x)[3])}))
              // parentHessians[[k]] <- array(NA, dim = c(dim(parentHessiansList[[1]])[1], dim(parentHessiansList[[1]])[1],
                                                       // parentHessiansDim))
              // thisDim <- 1
              // for(i_listEntry in 1:length(parentHessiansList)){
                // parentHessians[[k]][, , thisDim:(thisDim + dim(parentHessiansList[[i_listEntry]])[3] - 1)] <- parentHessiansList[[i_listEntry]]
                // thisDim <- thisDim +  dim(parentHessiansList[[i_listEntry]])[3]
              // }
            // }
          // }
        // }
        // if(!is.na(derivInfo$lineWrtArgsAsCharacters[[i]][1])){
          // derivList <- eval(substitute(nimDerivs(CALCCALL, DERIVORDERS, DROPARGS, WRT),
                                       // list(CALCCALL = derivInfo$calcWithArgsCalls[[i]],
                                            // DERIVORDERS = order,
                                            // DROPARGS = 'INDEXEDNODEINFO_',
                                            // WRT = derivInfo$lineWrtArgsAsCharacters[[i]])))

          // if(isDeterm){
            // derivList$value <- 0
            // model$nodeFunctions[[declID]]$calculate(unrolledIndicesMatrixRow)
          // }
          
          // ## The derivOutputFlag determines whether the derivatives of this node (node i): 
          // ## should be calculated for inclusion in the chain rule output (TRUE),
          // ## should be calculated for later use in the chain rule (FALSE), 
          // derivOutputFlag <- if(isDeterm) FALSE else TRUE
          
          // if(derivOutputFlag == TRUE){
            // ## If these derivatives will be included in output, we are taking derivative of a log prob. calculation,
            // ## so ouput will be length 1, and input args will be all wrt args.
            // chainRuleDerivList[[i]] <- matrix(0, nrow = 1, ncol = totalWrtSize)
            // if(hessianFlag) chainRuleHessianList[[i]] <- array(0, dim = c(totalWrtSize, totalWrtSize, 1))
          // }
          // else{
            // ## otherwise, we are taking derivative of a node value calculation (not log prob. calculation).
            // ## so ouput will be the length of this node, and input args will be all wrt args.
            // chainRuleDerivList[[i]] <- matrix(0, nrow = thisNodeSize, ncol = totalWrtSize)
            // if(hessianFlag) chainRuleHessianList[[i]] <- array(0, dim = c(totalWrtSize, totalWrtSize, thisNodeSize))
          // }
          // ## Iterate over all wrt params.
          // for(j in seq_along(wrtLineNums)){
            // thisArgIndex <- 0
            // ## Iterate over this line's parent nodes.
            // for(k in seq_along(derivInfo$parentIndicesList[[i]])){
              // if(!is.null(parentGradients[[k]])){
                // ## Calculate derivs of this node (i) wrt this parameter (j) for this parent node (k) via chain rule. 
                // chainRuleDerivList[[i]][,derivInfo$wrtLineIndices[[j]]] <- chainRuleDerivList[[i]][,derivInfo$wrtLineIndices[[j]]] +
                  // derivList$gradient[,(thisArgIndex + 1):(thisArgIndex + derivInfo$lineWrtArgSizeInfo[[i]][k]), drop = FALSE]%*%parentGradients[[k]][,derivInfo$wrtLineIndices[[j]], drop = FALSE]
                
              // }
              // thisArgIndex <- thisArgIndex + derivInfo$lineWrtArgSizeInfo[[i]][k]
            // }
            // if(derivOutputFlag == TRUE){
              // ## If this line is included in output, add the derivative of this line (i) wrt this param (j).
              // outDerivList$gradient[, derivInfo$wrtToIndices[[j]]] <- outDerivList$gradient[, derivInfo$wrtToIndices[[j]]]  +  
                // chainRuleDerivList[[i]][,derivInfo$wrtFromIndices[[j]]]
            // }
            // if(hessianFlag){
              // ## The Hessian is calculated below using Faà di Bruno's formula.
              // ## Second iteration over wrt parameters 
              // for(j_2 in j:length(wrtLineNums)){
                // thisArgIndex <- 0
                // ## Iterate over this line's parent nodes.
                // for(k in seq_along(derivInfo$parentIndicesList[[i]])){
                  // if(!is.null(parentHessians[[k]])){
                    // for(dim1 in derivInfo$WrtLineIndices[[j]]){
                      // for(dim2 in derivInfo$WrtLineIndices[[j_2]]){
                        // chainRuleHessianList[[i]][dim1, dim2, ] <- chainRuleHessianList[[i]][dim1, dim2, ] +
                          // c(derivList$gradient[ ,(thisArgIndex + 1):(thisArgIndex + derivInfo$lineWrtArgSizeInfo[[i]][k]), drop = FALSE]%*%parentHessians[[k]][dim1, dim2, , drop = FALSE])
                      // }
                    // }
                  // }
                  // thisArgIndex_2 <- 0
                  // for(k_2 in seq_along(derivInfo$parentIndicesList[[i]])){
                    // if(!is.null(parentGradients[[k]])){
                      // if(!is.null(parentGradients[[k_2]])){
                        // for(dim3 in 1:dim(derivList$hessian)[3]){
                          // chainRuleHessianList[[i]][derivInfo$wrtLineIndices[[j]], derivInfo$wrtLineIndices[[j_2]], dim3] <- chainRuleHessianList[[i]][derivInfo$wrtLineIndices[[j]], derivInfo$wrtLineIndices[[j_2]], dim3] +
                            // t(parentGradients[[k]][, derivInfo$wrtLineIndices[[j]], drop = FALSE])%*%derivList$hessian[(thisArgIndex + 1):(thisArgIndex + derivInfo$lineWrtArgSizeInfo[[i]][k]),(thisArgIndex_2 + 1):(thisArgIndex_2 + derivInfo$lineWrtArgSizeInfo[[i]][k_2]), dim3]%*%
                            // parentGradients[[k_2]][, derivInfo$wrtLineIndices[[j_2]], drop = FALSE]
                        // }
                      // }
                      
                    // }
                    // thisArgIndex_2 <- thisArgIndex_2 + derivInfo$lineWrtArgSizeInfo[[i]][k_2]
                  // }
                  // thisArgIndex <- thisArgIndex + derivInfo$lineWrtArgSizeInfo[[i]][k]
                // }
                // if(derivOutputFlag == TRUE){
                  // ## If this line is included in output, add the Hessian of this line (i) wrt this param #1 (j) and this param #2 (j_2).
                  // outDerivList$hessian[derivInfo$wrtToIndices[[j]], derivInfo$wrtToIndices[[j_2]], ] <-  outDerivList$hessian[derivInfo$wrtToIndices[[j]], derivInfo$wrtToIndices[[j_2]], ]   +
                    // chainRuleHessianList[[i]][derivInfo$wrtFromIndices[[j]], derivInfo$wrtFromIndices[[j_2]],]
                // }
              // }
            // }
          // }
        // }
        // else{ 
          // if(valueFlag) derivList <- eval(substitute(nimDerivs(CALCCALL, DERIVORDERS, DROPARGS),
                                                     // list(CALCCALL = derivInfo$calcWithArgsCalls[[i]],
                                                          // DERIVORDERS = c(0),
                                                          // DROPARGS = 'INDEXEDNODEINFO_')))
          // chainRuleDerivList[[i]] <- matrix(0, nrow = thisNodeSize, ncol = totalWrtSize)
          // if(hessianFlag) chainRuleHessianList[[i]] <- array(0, dim = c(totalWrtSize, totalWrtSize, thisNodeSize))
        // }
      // }
      
      // if(!isDeterm){
        // ## If this is a wrt node, we need to set the chainRule lists appropriately so that the chain rule
        // ## will work for dependent nodes of this node.  That means taking the first and second derivs of the
        // ## function f(x) = x, which will be the identity matrix and 0 respectively.
        // chainRuleDerivList[[i]] <- matrix(0,nrow = thisNodeSize, ncol = totalWrtSize)
        // if(isWrtLine)   chainRuleDerivList[[i]][,derivInfo$wrtLineIndices[[thisWrtLine]]] <- diag(thisNodeSize)
        // if(isWrtLine && hessianFlag) chainRuleHessianList[[i]] <- array(0, dim = c(totalWrtSize, totalWrtSize, thisNodeSize))
        // if(isCalcNodeLine){
          // if(valueFlag) outDerivList$value <- outDerivList$value + derivList$value
        // }
      // }
    // }
  // }
  
  // ## Reflect hessian across the diagonal
  // if(hessianFlag){
    // upperTriHess <- outDerivList$hessian[,,1]
    // upperTriHess[lower.tri(upperTriHess)] <-   t(upperTriHess)[lower.tri(upperTriHess)]
    // outDerivList$hessian[,,1] <-  upperTriHess
  // }
  // return(outDerivList)
// }



































nimSmartPtr<NIMBLE_ADCLASS>  NIM_DERIVS_CALCULATE(NodeVectorClassNew_derivs &nodes,  int iNodeFunction, NimArr<1, double> &derivOrders) {
  nimSmartPtr<NIMBLE_ADCLASS> ADlist = new NIMBLE_ADCLASS;
  // if(nodes.getInstructions().size() < static_cast<unsigned int>(iNodeFunction)) {
    // PRINTF("Warning in calculate: index of requested set of nodes is too large\n");
    // return(0);
  // }
  // const NodeInstruction &oneUseInfo = nodes.getInstructions()[iNodeFunction-1];
  // (*ADlist).value[1] += oneUseInfo.nodeFunPtr->calculateBlock(oneUseInfo.operand);
  return(ADlist); 
}
