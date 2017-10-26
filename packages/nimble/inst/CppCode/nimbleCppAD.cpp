#include <nimble/nimbleCppAD.h>


nimSmartPtr<NIMBLE_ADCLASS>  NIM_DERIVS_CALCULATE(NodeVectorClassNew_derivs &nodes, NimArr<1, double> &derivOrders) {
	nimSmartPtr<NIMBLE_ADCLASS> ansList = new NIMBLE_ADCLASS;
	nimSmartPtr<NIMBLE_ADCLASS> thisDerivList;
	const vector<NodeInstruction> &instructions = nodes.getInstructions();
	bool hessianFlag = false;
	bool jacobianFlag = false;
	bool valueFlag = false;

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
	if(!valueFlag && (!jacobianFlag && !hessianFlag)){
		printf("Error: must have at least one of 0, 1, or 2 in the 'order' argument.\n");
		return(ansList);
	}
	if(valueFlag){
		(*ansList).value.initialize(0, 1, 1);
	}
	if(valueFlag && (!jacobianFlag && !hessianFlag)){
		for(int i = 0; i < length(nodes.parentIndicesList); i++){
			if(nodes.calcNodeIndicators[i]){
				(*ansList).value[0] += instructions[i].nodeFunPtr->calculateBlock(instructions[i].operand);
			}
		}
		return(ansList);
	}
	derivOrders.setSize(valueFlag + jacobianFlag + hessianFlag);
	if(valueFlag){
		derivOrders[0] = 0;
	}
	if(jacobianFlag){
		derivOrders[valueFlag] = 1;
		if(hessianFlag){
			derivOrders[valueFlag + 1] = 2;
		}
	}
	(*ansList).gradient.initialize(0, 1, 1, nodes.totalOutWrtSize);
	if(hessianFlag){
		(*ansList).hessian.initialize(0, 1, nodes.totalOutWrtSize, nodes.totalOutWrtSize, 1);
	} 
	Map<MatrixXd> ansJacobian(0,0,0);
	new (&ansJacobian) Map< MatrixXd >((*ansList).gradient.getPtr(),(*ansList).gradient.dim()[0],(*ansList).gradient.dim()[1]);
	Map<MatrixXd> ansHessian(0,0,0);
	new (&ansHessian) Map< MatrixXd >((*ansList).hessian.getPtr(),(*ansList).hessian.dim()[0],(*ansList).hessian.dim()[1]);
	bool isDeterminisitic;
	bool isWrtLine;
	bool isCalcNodeLine;
	bool derivOutputFlag;
	int thisWrtLine;
	int thisNodeSize;
	int thisIndex;
	int thisRowLength;
	vector< MatrixXd > parentJacobians;
	vector<vector<vector<VectorXd > > > parentHessians;   // the parent Hessians are actually only used as vectors in our implementation of faa di brunos theorem, so stored as such.
	vector< MatrixXd > chainRuleJacobians(length(nodes.parentIndicesList));
	vector<vector< MatrixXd > > chainRuleHessians(length(nodes.parentIndicesList)); // chainRuleHessians index info:  first vector goes over nodes.parentIndicesList.  Second over third dimension of hessian.  Then individual matrices are first and second dims.
  
	for(int i = 0; i < length(nodes.parentIndicesList); i++){
		isDeterminisitic = 1 - nodes.stochNodeIndicators[i];
		isWrtLine = (nodes.cumulativeWrtLineNums[i] >= 0);
		isCalcNodeLine = nodes.calcNodeIndicators[i];
		thisWrtLine = nodes.cumulativeWrtLineNums[i];
		thisNodeSize = nodes.nodeLengths[i]; 			
		if(isCalcNodeLine){
			vector<int> thisWrtNodes(length(nodes.parentIndicesList[i]));
			parentJacobians.resize(length(nodes.parentIndicesList[i]));
			if(hessianFlag){
				parentHessians.resize(length(nodes.parentIndicesList[i]));
			} 
			for(int j = 0; j < length(nodes.parentIndicesList[i]); j++){
				// we can pre-calculate sumParentDims earlier.  Do this.
				if(j == 0 && isWrtLine){
					thisWrtNodes[j] = 1;
					parentJacobians[j] =  MatrixXd::Zero(nodes.wrtLineSize[thisWrtLine], nodes.totalWrtSize);
					for(int k = 0; k < nodes.wrtLineIndices[thisWrtLine].dimSize(0); k++){
						thisIndex = nodes.wrtLineIndices[thisWrtLine][k] - 1;
						parentJacobians[j](k, thisIndex) = 1;
					}
				}
				else if(nodes.parentIndicesList[i][j][0] > -1){
					thisWrtNodes[j] = 1;
					int sumParentDims = 0;
					for(int k = 0; k < nodes.parentIndicesList[i][j].dimSize(0); k++){
						sumParentDims += chainRuleJacobians[nodes.parentIndicesList[i][j][k]].rows();
					}
					parentJacobians[j] = MatrixXd::Zero(sumParentDims, nodes.totalWrtSize);
					int sumRowLengths = 0;
					for(int k = 0; k < nodes.parentIndicesList[i][j].dimSize(0); k++){
						thisRowLength = chainRuleJacobians[nodes.parentIndicesList[i][j][k]].rows();
						parentJacobians[j].block(sumRowLengths, 0, thisRowLength, nodes.totalWrtSize) = chainRuleJacobians[nodes.parentIndicesList[i][j][k]].block(0, 0, thisRowLength, nodes.totalWrtSize);
						sumRowLengths = sumRowLengths + thisRowLength;
					} 
					// if(hessianFlag){
						// parentHessians[j].resize(nodes.totalWrtSize);
						// for(int k = 0; k < nodes.totalWrtSize; k++){
							// parentHessians[j][k].resize(nodes.totalWrtSize);
							// for(int l = 0; l < nodes.totalWrtSize; l++){
								// parentHessians[j][k][l].resize(sumParentDims);
							// }
						// }
 						// for(int k = 0; k < nodes.parentIndicesList[i][j].dimSize(0); k++){
							// for(int l = 0; l <  chainRuleHessians[nodes.parentIndicesList[i][j][k]][0].rows(); l++){
								// for(int m = 0; m <  chainRuleHessians[nodes.parentIndicesList[i][j][k]][0].cols(); m++){
									// for(int n = 0; n < sumParentDims; n++){
										// parentHessians[j][l][m](n) = chainRuleHessians[nodes.parentIndicesList[i][j][k]][n](l, m);
									// }
								// }
							// }
						// } 
					// } 
				}
				else{				
					thisWrtNodes[j] = 0;
				}
			}
			if(nodes.cppWrtArgIndices[i][0] > -1){ 
				thisDerivList = instructions[i].nodeFunPtr->calculateWithArgs_derivBlock(instructions[i].operand, derivOrders, nodes.cppWrtArgIndices[i]);
				vector<Map<MatrixXd > > thisHessian; 
 				derivOutputFlag = (isDeterminisitic) ? false : true;
				if(derivOutputFlag){
					chainRuleJacobians[i] =  MatrixXd::Zero(1, nodes.totalWrtSize);
					if(hessianFlag){
						Map<MatrixXd> iHessian(0,0,0);
						new (&iHessian) EigenMapStrd((*thisDerivList).hessian.getPtr(), (*thisDerivList).hessian.dim()[0], (*thisDerivList).hessian.dim()[1],
							EigStrDyn((*thisDerivList).hessian.strides()[1], (*thisDerivList).hessian.strides()[0]));
						thisHessian.push_back(iHessian);
						chainRuleHessians[i].resize(1);
						chainRuleHessians[i][0] =  MatrixXd::Zero(nodes.totalWrtSize, nodes.totalWrtSize);
					}
				}
				else{
					chainRuleJacobians[i] =  MatrixXd::Zero(thisNodeSize, nodes.totalWrtSize);
					if(hessianFlag){
						chainRuleHessians[i].resize(thisNodeSize);
						for(int j = 0; j < thisNodeSize; j++){
							Map<MatrixXd> iHessian(0,0,0);
							new (&iHessian) EigenMapStrd((*thisDerivList).hessian.getPtr() + static_cast<int>(static_cast<int>(j * (*thisDerivList).hessian.strides()[2])),
								(*thisDerivList).hessian.dim()[0], (*thisDerivList).hessian.dim()[1], EigStrDyn((*thisDerivList).hessian.strides()[1], (*thisDerivList).hessian.strides()[0]));
							thisHessian.push_back(iHessian);
							chainRuleHessians[i][j] = MatrixXd::Zero(nodes.totalWrtSize, nodes.totalWrtSize);
						}
					}
				} 
 				for(int j = 0; j <  nodes.wrtLineNums.dimSize(0); j++){
					int thisArgIndex = 0;
					int wrtStartNode = nodes.wrtLineIndices[j][0] - 1;
					int wrtLength = nodes.wrtLineIndices[j].dimSize(0);
					int wrtToStartNode = nodes.wrtToIndices[j][0] - 1;
					int wrtToLength = nodes.wrtToIndices[j].dimSize(0);
					int wrtFromStartNode = nodes.wrtFromIndices[j][0] - 1;
					int wrtFromLength = nodes.wrtFromIndices[j].dimSize(0);
	 
					Map<MatrixXd> thisJacobian(0,0,0);
					new (&thisJacobian) Map< MatrixXd >((*thisDerivList).gradient.getPtr(),(*thisDerivList).gradient.dimSize(0),(*thisDerivList).gradient.dimSize(1));
					for(int k = 0; k < length(nodes.parentIndicesList[i]) ; k++){
						if(thisWrtNodes[k] == 1){
							//if(i < 4){
								// cout << "thisArgIndex: " <<  thisArgIndex << "\n";

								// cout << "nodes.lineWrtArgSizeInfo[i][k]: " << nodes.lineWrtArgSizeInfo[i][k] << "\n";
								// cout << "parentJacobians[k](0,0): " << parentJacobians[k](0,0) << "\n";	
								// cout << "parentJacobians[k](0,1): " << parentJacobians[k](0,1) << "\n";			
								// cout << "parentJacobians[k](0,2): " << parentJacobians[k](0,2) << "\n";			
								// cout << "parentJacobians[k](0,3): " << parentJacobians[k](0,3) << "\n";			
								
							 
								chainRuleJacobians[i].block(0, wrtStartNode, chainRuleJacobians[i].rows(), wrtLength) += (thisJacobian).block(0, thisArgIndex, (thisJacobian).rows(), nodes.lineWrtArgSizeInfo[i][k])*parentJacobians[k].block(0, wrtStartNode, parentJacobians[k].rows(), wrtLength);
							//}
						}
						thisArgIndex += nodes.lineWrtArgSizeInfo[i][k];
					}
					if(derivOutputFlag){
						// cout << "copying \n";
						// cout << "chainRuleJacobians[i](0, 0):" << chainRuleJacobians[i](0, 0) << "\n";
						// cout << "wrtToStartNode:" << wrtToStartNode << "\n";
						// cout << "wrtFromStartNode:" << wrtFromStartNode << "\n";
						// cout << "ansJacobian.rows():" << ansJacobian.rows() << "\n";
						// cout << "wrtToLength:" << wrtToLength << "\n";
						// cout << "wrtFromLength:" << wrtFromLength << "\n";
						// cout << "pre:" << i << " " << j << " " << ansJacobian.block(0, wrtToStartNode, ansJacobian.rows(), wrtToLength) << "\n";
						// if(i == 3 && j == 0){
							// cout << "wrtStartNode " << wrtStartNode << "\n";
							// cout << "wrtLength" << wrtLength << "\n";
						// }
						ansJacobian.block(0, wrtToStartNode, ansJacobian.rows(), wrtToLength) += chainRuleJacobians[i].block(0, wrtFromStartNode, ansJacobian.rows(), wrtFromLength);
						// // cout << "post:" << i << " " << j << " " << ansJacobian.block(0, wrtToStartNode, ansJacobian.rows(), wrtToLength) << "\n";

					} 
					// if(hessianFlag){
					 	// for(int j2 = 0; j2 <  nodes.wrtLineNums.dimSize(0); j2++){
							// thisArgIndex = 0;
							// int wrtStartNode2 = nodes.wrtLineIndices[j2][0] - 1;
							// int wrtLength2 = nodes.wrtLineIndices[j2].dimSize(0);
							// int wrtToStartNode2 = nodes.wrtToIndices[j2][0] - 1;
							// int wrtToLength2 = nodes.wrtToIndices[j2].dimSize(0);
							// int wrtFromStartNode2 = nodes.wrtFromIndices[j2][0] - 1;
							// int wrtFromLength2 = nodes.wrtFromIndices[j2].dimSize(0);
							// for(int k = 0; k < length(nodes.parentIndicesList[i]) ; k++){
								// if(thisWrtNodes[k] == 1){
									// for(int dim1 = 0; dim1 < nodes.wrtLineIndices[j].dimSize(0); dim1++){
										// for(int dim2 = 0; dim2 < nodes.wrtLineIndices[j2].dimSize(0); dim2++){
											// VectorXd addToHessian = (thisJacobian).block(0, thisArgIndex, (thisJacobian).rows(), nodes.lineWrtArgSizeInfo[i][k])*parentHessians[k][dim1][dim2];
											// for(int dim3 = 0; dim3 < length(chainRuleHessians[i]); dim3++){
												// chainRuleHessians[i][dim3](nodes.wrtLineIndices[j][dim1], nodes.wrtLineIndices[j2][dim2]) += addToHessian[dim3];
											// }
										// }
									// }
								// }
								// int thisArgIndex2 = 0;
								// for(int k2 = 0; k2 < length(nodes.parentIndicesList[i]) ; k2++){
									// if((thisWrtNodes[k] == 1) && (thisWrtNodes[k2] == 1)){
											// for(int dim3 = 0; dim3 < chainRuleJacobians[i].rows(); dim3++){
												// chainRuleHessians[i][dim3].block(wrtStartNode, wrtStartNode2, wrtLength, wrtLength2) += 
												// parentJacobians[k].block(0, wrtStartNode, parentJacobians[k].rows(), wrtLength).transpose()*
												// thisHessian[dim3].block(thisArgIndex, thisArgIndex2, nodes.lineWrtArgSizeInfo[i][k], nodes.lineWrtArgSizeInfo[i][k2])*
												// parentJacobians[k2].block(0, wrtStartNode2, parentJacobians[k2].rows(), wrtLength2);
											// }
									// }
									// thisArgIndex2 += nodes.lineWrtArgSizeInfo[i][k2];
								// }
								// thisArgIndex += nodes.lineWrtArgSizeInfo[i][k];
							// }
							// if(derivOutputFlag){
								// // cout << "chainRuleHessians[i][0](0, 0)" <<chainRuleHessians[i][0](0, 0) << "\n";
								// ansHessian.block(wrtToStartNode, wrtToStartNode2, wrtToLength, wrtToLength2) += chainRuleHessians[i][0].block(wrtFromStartNode, wrtFromStartNode2, wrtFromLength, wrtFromLength2);
							// }
						// }
					// }
				} 
 			}
			else{ 
				if(valueFlag){
					thisDerivList = instructions[i].nodeFunPtr->calculateWithArgs_derivBlock(instructions[i].operand, derivOrders, nodes.cppWrtArgIndices[i]);
				}
				chainRuleJacobians[i] = MatrixXd::Zero(thisNodeSize, nodes.totalWrtSize);
				if(hessianFlag){
					chainRuleHessians[i].resize(thisNodeSize);
					for(int j = 0; j < thisNodeSize; j++){
						chainRuleHessians[i][j] = MatrixXd::Zero(nodes.totalWrtSize, nodes.totalWrtSize);
					}
				} 
			}  
		} 
		if(!isDeterminisitic){ 
 			chainRuleJacobians[i] =  MatrixXd::Zero(thisNodeSize, nodes.totalWrtSize);
			if(isWrtLine){
				for(int j = 0; j < nodes.wrtLineIndices[thisWrtLine].dimSize(0); j++){
					thisIndex = nodes.wrtLineIndices[thisWrtLine][j] - 1;
					chainRuleJacobians[i](j, thisIndex) = 1;
				}
			}
			if(isWrtLine && hessianFlag){
				chainRuleHessians[i].resize(thisNodeSize);
				for(int j = 0; j < thisNodeSize; j++){
					chainRuleHessians[i][j] = MatrixXd::Zero(nodes.totalWrtSize, nodes.totalWrtSize);
				}
			}
			if(isCalcNodeLine && valueFlag){
				(*ansList).value[0] = (*ansList).value[0] + (*thisDerivList).value[0];
			} 
		}		 
	}
	if(hessianFlag){
		ansHessian.triangularView<Eigen::Lower>() = ansHessian.transpose().triangularView<Eigen::Lower>();
	}
	if(hessianFlag){
		ansHessian.triangularView<Eigen::Lower>() = ansHessian.transpose().triangularView<Eigen::Lower>();
	}

  return(ansList);  
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
