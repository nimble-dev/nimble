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
			vector<int> thisHessianNodes(length(nodes.parentIndicesList[i]));
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
					thisHessianNodes[j] = 1;
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
					if(hessianFlag){
						parentHessians[j].resize(nodes.totalWrtSize);
						for(int k = 0; k < nodes.totalWrtSize; k++){
							parentHessians[j][k].resize(nodes.totalWrtSize);
							for(int l = 0; l < nodes.totalWrtSize; l++){
								parentHessians[j][k][l].resize(sumParentDims);
							}
						}
 						for(int k = 0; k < nodes.parentIndicesList[i][j].dimSize(0); k++){
							for(int l = 0; l <  chainRuleHessians[nodes.parentIndicesList[i][j][k]][0].rows(); l++){
								for(int m = 0; m <  chainRuleHessians[nodes.parentIndicesList[i][j][k]][0].cols(); m++){
									for(int n = 0; n < sumParentDims; n++){
										parentHessians[j][l][m](n) = chainRuleHessians[nodes.parentIndicesList[i][j][k]][n](l, m);
									}
								}
							}
						} 
					} 
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
							chainRuleJacobians[i].block(0, wrtStartNode, chainRuleJacobians[i].rows(), wrtLength) += (thisJacobian).block(0, thisArgIndex, (thisJacobian).rows(), nodes.lineWrtArgSizeInfo[i][k])*parentJacobians[k].block(0, wrtStartNode, parentJacobians[k].rows(), wrtLength);
								
						}
						thisArgIndex += nodes.lineWrtArgSizeInfo[i][k];
					}
					if(derivOutputFlag){
						ansJacobian.block(0, wrtToStartNode, ansJacobian.rows(), wrtToLength) += chainRuleJacobians[i].block(0, wrtFromStartNode, ansJacobian.rows(), wrtFromLength);
					} 
					if(hessianFlag){
					 	for(int j2 = j; j2 <  nodes.wrtLineNums.dimSize(0); j2++){
							thisArgIndex = 0;
							int wrtStartNode2 = nodes.wrtLineIndices[j2][0] - 1;
							int wrtLength2 = nodes.wrtLineIndices[j2].dimSize(0);
							int wrtToStartNode2 = nodes.wrtToIndices[j2][0] - 1;
							int wrtToLength2 = nodes.wrtToIndices[j2].dimSize(0);
							int wrtFromStartNode2 = nodes.wrtFromIndices[j2][0] - 1;
							int wrtFromLength2 = nodes.wrtFromIndices[j2].dimSize(0);

							for(int k = 0; k < length(nodes.parentIndicesList[i]) ; k++){
								if(thisWrtNodes[k] == 1 && thisHessianNodes[k] == 1){
									for(int dim1 = 0; dim1 < nodes.wrtLineIndices[j].dimSize(0); dim1++){
										for(int dim2 = 0; dim2 < nodes.wrtLineIndices[j2].dimSize(0); dim2++){
											VectorXd addToHessian = (thisJacobian).block(0, thisArgIndex, (thisJacobian).rows(), nodes.lineWrtArgSizeInfo[i][k])*parentHessians[k][dim1][dim2];
											for(int dim3 = 0; dim3 < length(chainRuleHessians[i]); dim3++){
												chainRuleHessians[i][dim3](nodes.wrtLineIndices[j][dim1] - 1, nodes.wrtLineIndices[j2][dim2] - 1) += addToHessian[dim3];
											}
										}
									}
								}
								int thisArgIndex2 = 0;
								for(int k2 = 0; k2 < length(nodes.parentIndicesList[i]) ; k2++){
									if(thisWrtNodes[k] == 1 && thisWrtNodes[k2] == 1){
											for(int dim3 = 0; dim3 < chainRuleJacobians[i].rows(); dim3++){
												chainRuleHessians[i][dim3].block(wrtStartNode, wrtStartNode2, wrtLength, wrtLength2) += 
												parentJacobians[k].block(0, wrtStartNode, parentJacobians[k].rows(), wrtLength).transpose()*
												thisHessian[dim3].block(thisArgIndex, thisArgIndex2, nodes.lineWrtArgSizeInfo[i][k], nodes.lineWrtArgSizeInfo[i][k2])*
												parentJacobians[k2].block(0, wrtStartNode2, parentJacobians[k2].rows(), wrtLength2);
											}
									}
									thisArgIndex2 += nodes.lineWrtArgSizeInfo[i][k2];
								}
								thisArgIndex += nodes.lineWrtArgSizeInfo[i][k];
							}
							if(derivOutputFlag){
								ansHessian.block(wrtToStartNode, wrtToStartNode2, wrtToLength, wrtToLength2) += chainRuleHessians[i][0].block(wrtFromStartNode, wrtFromStartNode2, wrtFromLength, wrtFromLength2);
							}
						}
					}
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
  return(ansList);  
}

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
