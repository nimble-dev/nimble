#include <nimble/nimbleCppAD.h>

nimSmartPtr<NIMBLE_ADCLASS> NIM_DERIVS_CALCULATE(
    NodeVectorClassNew_derivs &nodes, NimArr<1, double> &derivOrders) {
  nimSmartPtr<NIMBLE_ADCLASS> ansList =
      new NIMBLE_ADCLASS;  // This will be returned from this funciton.
  nimSmartPtr<NIMBLE_ADCLASS>
      thisDerivList;  // Used to store derivative output from individual nodes.
  const vector<NodeInstruction> &instructions = nodes.getInstructions();
  bool hessianFlag = false;   // Are second order derivs requested?
  bool jacobianFlag = false;  // Are first order derivs requested?
  bool valueFlag = false;     // Is the function value requested?

  for (int i = 0; i < derivOrders.dimSize(0); i++) {
    if (derivOrders[i] == 0) {
      valueFlag = true;
    }
    if (derivOrders[i] == 1) {
      jacobianFlag = true;
    }
    if (derivOrders[i] == 2) {
      hessianFlag = true;
      jacobianFlag = true;  // If second order is asked for, first order must be
                            // calculated as well (although we could technically
                            // skip copying this into output)
    }
  }
  if (!valueFlag && (!jacobianFlag && !hessianFlag)) {
    printf(
        "Error: must have at least one of 0, 1, or 2 in the 'order' "
        "argument.\n");
    return (ansList);
  }
  if (valueFlag) {
    (*ansList).value.initialize(0, 1, 1);
  }
  if (valueFlag && (!jacobianFlag && !hessianFlag)) {  // If only function value
                                                       // is asked for, skip
                                                       // derivative calculation
    for (int i = 0; i < length(nodes.parentIndicesList); i++) {
      if (nodes.calcNodeIndicators[i]) {
        (*ansList).value[0] +=
            instructions[i].nodeFunPtr->calculateBlock(instructions[i].operand);
      }
    }
    return (ansList);
  }
  derivOrders.setSize(valueFlag + jacobianFlag + hessianFlag);
  if (valueFlag) {
    derivOrders[0] = 0;
  }
  if (jacobianFlag) {
    derivOrders[valueFlag] = 1;
    if (hessianFlag) {
      derivOrders[valueFlag + 1] = 2;
    }
  }
  // Initialize the Jacobian and Hessian of the returned list
  (*ansList).gradient.initialize(0, 1, 1, nodes.totalOutWrtSize);
  if (hessianFlag) {
    (*ansList).hessian.initialize(0, 1, nodes.totalOutWrtSize,
                                  nodes.totalOutWrtSize, 1);
  }
  Map<MatrixXd> ansJacobian((*ansList).gradient.getPtr(),
                            (*ansList).gradient.dim()[0],
                            (*ansList).gradient.dim()[1]);
  Map<MatrixXd> ansHessian((*ansList).hessian.getPtr(),
                           (*ansList).hessian.dim()[0],
                           (*ansList).hessian.dim()[1]);
  bool isDeterminisitic;
  bool isWrtLine;
  bool isCalcNodeLine;
  bool derivOutputFlag;
  int thisWrtLine;
  int thisNodeSize;
  int thisIndex;
  int thisRowLength;
  vector<MatrixXd> parentJacobians;  // For each node, parentJacobians is
                                     // repopulated to store the Jacobian
                                     // matrices of all of the parent nodes of
                                     // that node.

  vector<vector<vector<VectorXd> > > parentHessians;  // Similar to
                                                      // parentJacobians.  The
                                                      // parent Hessians are
                                                      // actually only used as
                                                      // vectors in our
                                                      // implementation of Faa
                                                      // Di Bruno's theorem, so
                                                      // stored as such. E.g.
                                                      // parentHessians[i][j][k]
                                                      // would be the i'th
                                                      // parent node's Hessian
                                                      // array's third dimension
                                                      // (a vector).

  vector<MatrixXd> chainRuleJacobians(length(
      nodes.parentIndicesList));  // Similar storage as parentJacobians above.
                                  // chainRuleJacobians[i] is the Jacobian
                                  // matrix for the i'th node (as determined by
                                  // nodes.parentIndicesList). Values of
                                  // chainRuleJacobians are not overwritten, but
                                  // stored for use in populating the
                                  // parentJacobians of downstream nodes.
  vector<vector<MatrixXd> > chainRuleHessians(  // Due to the linear algebra
                                                // necessary for chain rule
                                                // calculations,
      length(nodes.parentIndicesList));  // the chainRuleHessians are stored in
                                         // a different format than the parent
                                         // Hessians. The outer vector goes over
                                         // nodes.parentIndicesList. The second
                                         // vector goes over the third dimension
                                         // of the Hessian (the third dimension
                                         // has length = sum length of wrt
                                         // arguments). Then individual
                                         // MatrixXds are first and second dims
                                         // of the Hessian.

  for (int i = 0; i < length(nodes.parentIndicesList); i++) {
    isDeterminisitic =
        1 - nodes.stochNodeIndicators[i];  // Is node i deterministic?
    isWrtLine = (nodes.cumulativeWrtLineNums[i] >= 0);  // Is node i a wrt node
                                                        // (i.e. was it included
                                                        // in the wrt argument
                                                        // to nimDerivs?)
    isCalcNodeLine = nodes.calcNodeIndicators[i];  // Is node i a node included
                                                   // in the call to
                                                   // calculate(nodes) ?
    thisWrtLine = nodes.cumulativeWrtLineNums[i];  // If node i is a wrt node,
                                                   // which wrt node is it (i.e.
                                                   // where does it fall in the
                                                   // topological ordering of
                                                   // wrt nodes) ?
    thisNodeSize = nodes.nodeLengths[i];  // Total size (= product of lengths of
                                          // dimensions) of this node
    if (isCalcNodeLine) {
      // If this node was included in the call to calculate(nodes), we first
      // populate its parent derivative information.

      vector<int> thisWrtNodes(
          length(nodes.parentIndicesList[i]));  // for each parent node,  do we
                                                // need to take derivatives wrt
                                                // that node?
      vector<int> thisHessianNodes(
          length(nodes.parentIndicesList[i]));  // for each parent node,  do we
                                                // need to calculate hessians
                                                // wrt that node?
      parentJacobians.resize(length(nodes.parentIndicesList[i]));
      if (hessianFlag) {
        parentHessians.resize(length(nodes.parentIndicesList[i]));
      }
      for (int j = 0; j < length(nodes.parentIndicesList[i]);
           j++) {  // Iterate over all arguments to calcWithArgs() for node i
                   // and populate parent deriv info.
                   // we can pre-calculate sumParentDims earlier.  Do this.
        if (j == 0 && isWrtLine) {  // This is the case where this node (node i)
                                    // was included in the wrt argument to
                                    // nimDerivs.  Note that for the purposes of
                                    // the chain rule, we set the parent 1st
                                    // derivatives of this node wrt itself (its
                                    // first argument) as a matrix with 1s on
                                    // the diagonal.
          thisWrtNodes[j] = 1;
          parentJacobians[j] = MatrixXd::Zero(nodes.wrtLineSize[thisWrtLine],
                                              nodes.totalWrtSize);
          for (int k = 0; k < nodes.wrtLineIndices[thisWrtLine].dimSize(0);
               k++) {
            thisIndex = nodes.wrtLineIndices[thisWrtLine][k] - 1;
            parentJacobians[j](k, thisIndex) = 1;
          }
        } else if (nodes.parentIndicesList[i][j][0] >
                   -1) {  // This is the case where at least one of the parent
                          // nodes of the j'th argument to this node (node i)'s
                          // calculateWithArgs() function depends on a wrt node.
                          // We populate parentJacobians and parentHessians by
                          // retreiving the parent nodes' already calculated
                          // derivatives from the chainRuleJacobians and
                          // chainRuleHessians.  In  the case where the j'th
                          // argument depends on multiple nodes, we piece
                          // together the parentJacobians and parentHessians
                          // from multiple entires in chainRuleJacobians and
                          // chainRuleHessians.
          thisHessianNodes[j] = 1;
          thisWrtNodes[j] = 1;
          vector<int> parentDims(nodes.parentIndicesList[i][j].dimSize(0));
          int sumParentDims = 0;
          for (int k = 0; k < nodes.parentIndicesList[i][j].dimSize(0); k++) {
            parentDims[k] =
                chainRuleJacobians[nodes.parentIndicesList[i][j][k]].rows();
            sumParentDims +=  parentDims[k];
          }
          parentJacobians[j] =
              MatrixXd::Zero(sumParentDims, nodes.totalWrtSize);
          int sumRowLengths = 0;
          for (int k = 0; k < nodes.parentIndicesList[i][j].dimSize(0);
               k++) {  // Populating Jacobians and Hessians for parents is
                       // somewhat complicated by the fact that
                       //  argument j can come from multiple different nodes.
                       //  I.e., in BUGS code: x[1] ~ dnorm(...) x[2] ~
                       //  dnorm(...) y[1:2] ~ dmnorm(x[1:2], ...)
                       // Here the second argument to y[1:2]'s calculateWithArgs
                       // function would be x[1:2], and the parent derivative
                       // information would have to come from both x[1] and
                       // x[2].
                       //  The derivative information from x[1] and x[2] are
                       //  combined into a single (matrix / array) to enable the
                       //  linear algebra necessary to calculate the (Jacobian /
                       //  Hessian).
            thisRowLength =
                chainRuleJacobians[nodes.parentIndicesList[i][j][k]].rows();
            parentJacobians[j].block(sumRowLengths, 0, thisRowLength,
                                     nodes.totalWrtSize) =
                chainRuleJacobians[nodes.parentIndicesList[i][j][k]].block(
                    0, 0, thisRowLength, nodes.totalWrtSize);
            sumRowLengths = sumRowLengths + thisRowLength;
          }
          if (hessianFlag) {
            parentHessians[j].resize(nodes.totalWrtSize);
            for (int k = 0; k < nodes.totalWrtSize; k++) {
              parentHessians[j][k].resize(nodes.totalWrtSize);
              for (int l = 0; l < nodes.totalWrtSize; l++) {
                parentHessians[j][k][l].resize(sumParentDims);
              }
            }
            for (int k = 0; k < nodes.parentIndicesList[i][j].dimSize(0); k++) {
              for (int l = 0;
                   l < chainRuleHessians[nodes.parentIndicesList[i][j][k]][0]
                           .rows();
                   l++) {
                for (int m = 0;
                     m < chainRuleHessians[nodes.parentIndicesList[i][j][k]][0]
                             .cols();
                     m++) {
                  int indexTracker = 0;     
                  for (int n = 0; n < parentDims[k]; n++) {
                    parentHessians[j][l][m](indexTracker) = 
                        chainRuleHessians[nodes.parentIndicesList[i][j][k]][n](
                            l, m);
                    indexTracker++;
                  }
                }
              }
            }
          }
        } else {  // Otherwise, the j'th argument has no parents that depend on
                  // wrt nodes.
          thisWrtNodes[j] = 0;
        }
      }
      if (nodes.cppWrtArgIndices[i][0] >
          -1) {  // -1 is used as an indicator that node i doesn't have any
                 // parents that depend on WRT arguments, so its derivs don't
                 // need too be taken.  Otherwise, we proceed with the chain
                 // rule.
        thisDerivList =
            instructions[i].nodeFunPtr->calculateWithArgs_derivBlock(
                instructions[i].operand, derivOrders,
                nodes.cppWrtArgIndices[i]);  // Derivatives of calculate() for
                                             // node i are computed here.
        vector<Map<MatrixXd, Unaligned, EigStrDyn> > thisHessian;
        derivOutputFlag =
            (isDeterminisitic) ? false : true;  // derivOutputFlag is true if
                                                // the derivative output from
                                                // this node will be added to
                                                // the ansList.
        if (derivOutputFlag) {  // If derivOutputFlag == true, these derivatives
                                // will be taken of a function that returns a
                                // scalar (a density calculation), so the
                                // storage is sized appropriately.
          chainRuleJacobians[i] = MatrixXd::Zero(1, nodes.totalWrtSize);
          if (hessianFlag) {
            Map<MatrixXd, Unaligned, EigStrDyn> iHessian(
                (*thisDerivList).hessian.getPtr(),
                (*thisDerivList).hessian.dim()[0],
                (*thisDerivList).hessian.dim()[1],
                EigStrDyn((*thisDerivList).hessian.strides()[1],
                          (*thisDerivList).hessian.strides()[0]));
            thisHessian.push_back(iHessian);
            chainRuleHessians[i].resize(1);
            chainRuleHessians[i][0] =
                MatrixXd::Zero(nodes.totalWrtSize, nodes.totalWrtSize);
          }
        } else {  // If derivOutputFlag == false, these derivatives will be
                  // taken of a deterministic function that returns an object of
                  // thisNodeSize, so the storage is sized appropriately.
          chainRuleJacobians[i] =
              MatrixXd::Zero(thisNodeSize, nodes.totalWrtSize);
          if (hessianFlag) {
            chainRuleHessians[i].resize(thisNodeSize);
            for (int j = 0; j < thisNodeSize; j++) {
              Map<MatrixXd, Unaligned, EigStrDyn> iHessian(
                  (*thisDerivList).hessian.getPtr() +
                      static_cast<int>(static_cast<int>(
                          j * (*thisDerivList).hessian.strides()[2])),
                  (*thisDerivList).hessian.dim()[0],
                  (*thisDerivList).hessian.dim()[1],
                  EigStrDyn((*thisDerivList).hessian.strides()[1],
                            (*thisDerivList).hessian.strides()[0]));
              thisHessian.push_back(iHessian);
              chainRuleHessians[i][j] =
                  MatrixXd::Zero(nodes.totalWrtSize, nodes.totalWrtSize);
            }
          }
        }
        for (int j = 0; j < nodes.wrtLineNums.dimSize(0);
             j++) {  // Next we iterate over all the wrt nodes.  For each wrt
                     // node, we iterate over all of the nodes that are
                     // arguments to this node's calculateWithArgs() function.
                     // If argument k depends on wrt node j, we calcuate the
                     // first and second derivatives (stored in
                     // chainRuleJacobians[i] and chainRuleHessians[i]) of this
                     // node (node i) wrt wrt node j using Faa di Bruno's
                     // formula.
          int thisArgIndex = 0;
          int wrtStartNode = nodes.wrtLineIndices[j][0] - 1;
          int wrtLength = nodes.wrtLineIndices[j].dimSize(0);
          int wrtToStartNode = nodes.wrtToIndices[j][0] - 1;
          int wrtToLength = nodes.wrtToIndices[j].dimSize(0);
          int wrtFromStartNode = nodes.wrtFromIndices[j][0] - 1;
          int wrtFromLength = nodes.wrtFromIndices[j].dimSize(0);
          Map<MatrixXd> thisJacobian((*thisDerivList).gradient.getPtr(),
                                     (*thisDerivList).gradient.dimSize(0),
                                     (*thisDerivList).gradient.dimSize(1));

          for (int k = 0; k < length(nodes.parentIndicesList[i]); k++) {
            if (thisWrtNodes[k] == 1) {
              chainRuleJacobians[i].block(
                  0, wrtStartNode, chainRuleJacobians[i].rows(), wrtLength) +=
                  (thisJacobian)
                      .block(0, thisArgIndex, (thisJacobian).rows(),
                             nodes.lineWrtArgSizeInfo[i][k]) *
                  parentJacobians[k].block(
                      0, wrtStartNode, parentJacobians[k].rows(), wrtLength);
            }
            thisArgIndex += nodes.lineWrtArgSizeInfo[i][k];
          }
          if (derivOutputFlag) {  // Add this Jacobian to the output of this
                                  // function
            ansJacobian.block(0, wrtToStartNode, ansJacobian.rows(),
                              wrtToLength) +=
                chainRuleJacobians[i].block(0, wrtFromStartNode,
                                            ansJacobian.rows(), wrtFromLength);
          }
          if (hessianFlag) {
            for (int j2 = j; j2 < nodes.wrtLineNums.dimSize(0); j2++) {
              thisArgIndex = 0;
              int wrtStartNode2 = nodes.wrtLineIndices[j2][0] - 1;
              int wrtLength2 = nodes.wrtLineIndices[j2].dimSize(0);
              int wrtToStartNode2 = nodes.wrtToIndices[j2][0] - 1;
              int wrtToLength2 = nodes.wrtToIndices[j2].dimSize(0);
              int wrtFromStartNode2 = nodes.wrtFromIndices[j2][0] - 1;
              int wrtFromLength2 = nodes.wrtFromIndices[j2].dimSize(0);

              for (int k = 0; k < length(nodes.parentIndicesList[i]); k++) {
                if (thisWrtNodes[k] == 1 && thisHessianNodes[k] == 1) {
                  for (int dim1 = 0; dim1 < nodes.wrtLineIndices[j].dimSize(0);
                       dim1++) {
                    for (int dim2 = 0;
                         dim2 < nodes.wrtLineIndices[j2].dimSize(0); dim2++) {
                      VectorXd addToHessian =
                          (thisJacobian)
                              .block(0, thisArgIndex, (thisJacobian).rows(),
                                     nodes.lineWrtArgSizeInfo[i][k]) *
                          parentHessians[k][dim1][dim2];
                      for (int dim3 = 0; dim3 < length(chainRuleHessians[i]);
                           dim3++) {
                        chainRuleHessians[i][dim3](
                            nodes.wrtLineIndices[j][dim1] - 1,
                            nodes.wrtLineIndices[j2][dim2] - 1) +=
                            addToHessian[dim3];
                      }
                    }
                  }
                }
                int thisArgIndex2 = 0;
                for (int k2 = 0; k2 < length(nodes.parentIndicesList[i]);
                     k2++) {
                  if (thisWrtNodes[k] == 1 && thisWrtNodes[k2] == 1) {
                    for (int dim3 = 0; dim3 < chainRuleJacobians[i].rows();
                         dim3++) {
                      chainRuleHessians[i][dim3].block(
                          wrtStartNode, wrtStartNode2, wrtLength, wrtLength2) +=
                          parentJacobians[k]
                              .block(0, wrtStartNode, parentJacobians[k].rows(),
                                     wrtLength)
                              .transpose() *
                          thisHessian[dim3].block(
                              thisArgIndex, thisArgIndex2,
                              nodes.lineWrtArgSizeInfo[i][k],
                              nodes.lineWrtArgSizeInfo[i][k2]) *
                          parentJacobians[k2].block(0, wrtStartNode2,
                                                    parentJacobians[k2].rows(),
                                                    wrtLength2);
                    }
                  }
                  thisArgIndex2 += nodes.lineWrtArgSizeInfo[i][k2];
                }
                thisArgIndex += nodes.lineWrtArgSizeInfo[i][k];
              }
              if (derivOutputFlag) {  // Add this Hessian to the output of this
                                      // function
                ansHessian.block(wrtToStartNode, wrtToStartNode2, wrtToLength,
                                 wrtToLength2) +=
                    chainRuleHessians[i][0].block(
                        wrtFromStartNode, wrtFromStartNode2, wrtFromLength,
                        wrtFromLength2);
              }
            }
          }
        }
      } else {  // Otherwise, no arguments depend on wrt nodes, so we don't need
                // to take derivatives. Instead, just set Jacobian and Hessian
                // to 0.
        if (valueFlag) {  // Still may need to get the value (0'th order deriv).
          NimArr<1, double> valueOrder(1);
          valueOrder[0] = 0;
          thisDerivList =
              instructions[i].nodeFunPtr->calculateWithArgs_derivBlock(
                  instructions[i].operand, valueOrder,
                  nodes.cppWrtArgIndices[i]);
        }
        chainRuleJacobians[i] =
            MatrixXd::Zero(thisNodeSize, nodes.totalWrtSize);
        if (hessianFlag) {
          chainRuleHessians[i].resize(thisNodeSize);
          for (int j = 0; j < thisNodeSize; j++) {
            chainRuleHessians[i][j] =
                MatrixXd::Zero(nodes.totalWrtSize, nodes.totalWrtSize);
          }
        }
      }
    }
    if (!isDeterminisitic) {  // If this is a wrt node, we need to set the
                              // chainRule derivatives appropriately so that the
                              // chain rule will work for dependent nodes of
                              // this node.  That means taking the first and
                              // second derivs of the function f(x) = x, which
                              // will be the identity matrix and 0 respectively.
      chainRuleJacobians[i] = MatrixXd::Zero(thisNodeSize, nodes.totalWrtSize);
      if (isWrtLine) {
        for (int j = 0; j < nodes.wrtLineIndices[thisWrtLine].dimSize(0); j++) {
          thisIndex = nodes.wrtLineIndices[thisWrtLine][j] - 1;
          chainRuleJacobians[i](j, thisIndex) = 1;
        }
      }
      if (hessianFlag) {
        chainRuleHessians[i].resize(thisNodeSize);
        for (int j = 0; j < thisNodeSize; j++) {
          chainRuleHessians[i][j] =
              MatrixXd::Zero(nodes.totalWrtSize, nodes.totalWrtSize);
        }
      }
      if (isCalcNodeLine && valueFlag) {
        (*ansList).value[0] = (*ansList).value[0] + (*thisDerivList).value[0];
      }
    }
  }
  if (hessianFlag) {  // reflect Hessian across the diagonal
    ansHessian.triangularView<Eigen::Lower>() =
        ansHessian.transpose().triangularView<Eigen::Lower>();
  }
  return (ansList);
}

nimSmartPtr<NIMBLE_ADCLASS> NIM_DERIVS_CALCULATE(
    NodeVectorClassNew_derivs &nodes, int iNodeFunction,
    NimArr<1, double> &derivOrders) {
  nimSmartPtr<NIMBLE_ADCLASS> ADlist = new NIMBLE_ADCLASS;
  // if(nodes.getInstructions().size() < static_cast<unsigned
  // int>(iNodeFunction)) {
  // PRINTF("Warning in calculate: index of requested set of nodes is too
  // large\n");
  // return(0);
  // }
  // const NodeInstruction &oneUseInfo =
  // nodes.getInstructions()[iNodeFunction-1];
  // (*ADlist).value[1] +=
  // oneUseInfo.nodeFunPtr->calculateBlock(oneUseInfo.operand);
  return (ADlist);
}
