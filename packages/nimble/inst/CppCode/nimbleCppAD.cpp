#include <nimble/nimbleCppAD.h>

#ifdef _TIME_AD
ad_timer derivs_main_timer("derivs_main");
ad_timer derivs_calc_timer("derivs_calc");
ad_timer derivs_node_iteration_timer("derivs_node_iteration");
ad_timer derivs_getDerivs_timer("derivs_getDerivs");
ad_timer derivs_run_tape_timer("derivs_run_tape");
int id(0);
// Include the following line to see node by node timing timing
// results within use of getDerivs and use of CppAD tapes.  This
// generates a lot of output.
//#define _SHOW_NODE_BY_NODE 
void derivs_getDerivs_timer_start() {derivs_getDerivs_timer.start(false);}
void derivs_run_tape_timer_start() {derivs_run_tape_timer.start(false);}
#ifdef _SHOW_NODE_BY_NODE 
void derivs_getDerivs_timer_stop() {derivs_getDerivs_timer.stop(true);}
void derivs_run_tape_timer_stop() {derivs_run_tape_timer.stop(true);}
void derivs_show_id() {std::cout<<"(id "<<id<<")"<<std::endl;}
#else
void derivs_getDerivs_timer_stop() {derivs_getDerivs_timer.stop(false);}
void derivs_run_tape_timer_stop() {derivs_run_tape_timer.stop(false);}
void derivs_show_id() {}
#endif
  void derivs_tick_id() {++id;}

SEXP report_AD_timers() {
  derivs_main_timer.show_report();
  derivs_calc_timer.show_report();
  derivs_node_iteration_timer.show_report();
  derivs_getDerivs_timer.show_report();
  derivs_run_tape_timer.show_report();
  return(R_NilValue);
}

SEXP reset_AD_timers(SEXP SreportInterval) {
  derivs_main_timer.reset();
  derivs_calc_timer.reset();
  derivs_node_iteration_timer.reset();
  derivs_getDerivs_timer.reset();
  derivs_run_tape_timer.reset();
  derivs_main_timer.set_interval(INTEGER(SreportInterval)[0]);
  derivs_calc_timer.set_interval(INTEGER(SreportInterval)[0]);
  derivs_node_iteration_timer.set_interval(INTEGER(SreportInterval)[0]);
  derivs_getDerivs_timer.set_interval(INTEGER(SreportInterval)[0]);
  derivs_run_tape_timer.set_interval(INTEGER(SreportInterval)[0]);
  return(R_NilValue);
}

#endif

#define _DERIVS_FULLTAPE
#ifdef  _DERIVS_FULLTAPE

nimSmartPtr<NIMBLE_ADCLASS> NIM_DERIVS_CALCULATE(
						 NodeVectorClassNew_derivs &nodes,
						 const NimArr<1, double> &derivOrders) {
  // std::cout<<"need to propagate const-ness"<<std::endl;

  if(!nodes.tapeRecorded()) nodes.recordTape();
#ifdef _TIME_AD
  derivs_main_timer.start();
#endif

  nimSmartPtr<NIMBLE_ADCLASS> ansList = new NIMBLE_ADCLASS;
#ifdef _TIME_AD
  derivs_getDerivs_timer.start();
#endif
  std::vector<double> independentVars;
  nodes.runTape_setIndependent(independentVars);
  std::vector<double> dependentVars;

#ifdef _TIME_AD
  derivs_run_tape_timer.start();
#endif
    
  nodes.runTape_runTape(independentVars, dependentVars);

#ifdef _TIME_AD
  derivs_run_tape_timer.stop();
#endif

#ifdef _TIME_AD
  derivs_getDerivs_timer.stop();
#endif

  ansList->value.setSize(1);
  //  ansList->value[0] = nodes.runTape(derivOrders);
  ansList->value[0] = nodes.runTape_unpackDependent(dependentVars);

#ifdef _TIME_AD
  derivs_main_timer.stop();
#endif

  return ansList;
}

#else

// nimSmartPtr<NIMBLE_ADCLASS> NIM_DERIVS_CALCULATE(
// 						 const NodeVectorClassNew_derivs &nodes,
// 						 const NimArr<1, double> &derivOrders) {
// #ifdef _TIME_AD
//   derivs_main_timer.start();
// #endif
  
//   nimSmartPtr<NIMBLE_ADCLASS> ansList =
//       new NIMBLE_ADCLASS;  // This will be returned from this funciton.
//   nimSmartPtr<NIMBLE_ADCLASS> thisDerivList =
//       new NIMBLE_ADCLASS;  // Used to store derivative output from individual
//                            // nodes.
//   const vector<NodeInstruction> &instructions = nodes.getConstInstructions();
//   bool hessianFlag = false;   // Are second order derivs requested?
//   bool jacobianFlag = false;  // Are first order derivs requested?
//   bool valueFlag = false;     // Is the function value requested?

//   for (int i = 0; i < derivOrders.dimSize(0); i++) {
//     if (derivOrders[i] == 0) {
//       valueFlag = true;
//     }
//     if (derivOrders[i] == 1) {
//       jacobianFlag = true;
//     }
//     if (derivOrders[i] == 2) {
//       hessianFlag = true;
//       jacobianFlag = true;  // If second order is asked for, first order must be
//                             // calculated as well (although we could technically
//                             // skip copying this into output)
//     }
//   }

//   if (!valueFlag && (!jacobianFlag && !hessianFlag)) {
//     printf(
//         "Error: must have at least one of 0, 1, or 2 in the 'order' "
//         "argument.\n");
//     return (ansList);
//   }
//   if (valueFlag) {
//     (*ansList).value.initialize(0, 1, 1);
//   }
//   int iLength;
//   if (valueFlag && (!jacobianFlag && !hessianFlag)) {  // If only function value
//                                                        // is asked for, skip
//     iLength = nodes.parentIndicesList.size();          // derivative calculation
//     for (int i = 0; i < iLength; i++) {
//       if (nodes.calcNodeIndicators[i]) {
//         (*ansList).value[0] +=
//             instructions[i].nodeFunPtr->calculateBlock(instructions[i].operand);
//       }
//     }
//     return (ansList);
//   }
//   NimArr<1, double> newDerivOrders;
//   newDerivOrders.setSize(valueFlag + jacobianFlag + hessianFlag);
//   if (valueFlag) {
//     newDerivOrders[0] = 0;
//   }
//   if (jacobianFlag) {
//     newDerivOrders[valueFlag] = 1;
//     if (hessianFlag) {
//       newDerivOrders[valueFlag + 1] = 2;
//     }
//   }

//   // Initialize the Jacobian and Hessian of the returned list
//   (*ansList).jacobian.initialize(0, 1, 1, nodes.totalOutWrtSize);
//   if (hessianFlag) {
//     (*ansList).hessian.initialize(0, 1, nodes.totalOutWrtSize,
//                                   nodes.totalOutWrtSize, 1);
//   }
//   Map<MatrixXd> ansJacobian((*ansList).jacobian.getPtr(),
//                             (*ansList).jacobian.dim()[0],
//                             (*ansList).jacobian.dim()[1]);
//   Map<MatrixXd> ansHessian((*ansList).hessian.getPtr(),
//                            (*ansList).hessian.dim()[0],
//                            (*ansList).hessian.dim()[1]);
//   bool isDeterminisitic;
//   bool isWrtLine;
//   bool isCalcNodeLine;
//   bool derivOutputFlag;
//   int thisWrtLine;
//   int thisNodeSize;
//   int jLength;
//   int kLength;
//   int k2Length;
//   int lLength;
//   int dim3Length;
//   int lineWrtSizeIJ;
//   int wrtNodeK;
//   int thisRows;
//   int thisCols;
//   int thisArgIndex;
//   int wrtNodeJ;

//   vector<MatrixXd> chainRuleJacobians(
//       nodes.parentIndicesList
//           .size()); // chainRuleJacobians[i] is the Jacobian
//                      // matrix for the i'th node (as determined by
//                      // nodes.parentIndicesList). Values of
//                      // chainRuleJacobians are not overwritten, but
//                      // stored for use in populating the
//                      // parentJacobians of downstream nodes.

//   vector<vector<MatrixXd> > chainRuleHessians(  
//       nodes.parentIndicesList.size());  // The outer vector goes
//   //                                        over
//   //                                        // nodes.parentIndicesList. The
//   //                                        second
//   //                                        // vector goes over the third
//   //                                        dimension
//   //                                        // of the Hessian (the third
//   //                                        dimension
//   //                                        // has length = sum length of wrt
//   //                                        // arguments). Then individual
//   //                                        // MatrixXds are first and second
//   //                                        dims
//   //                                        // of the Hessian.
//   iLength = nodes.parentIndicesList.size();
//   int sumAddedScalarNodes = 0;
//   for (int i = 0; i < iLength; i++) {
// #ifdef _TIME_AD
//     derivs_node_iteration_timer.start();
// #endif
//     isDeterminisitic =
//         1 - nodes.stochNodeIndicators[i];  // Is node i deterministic?
//     isWrtLine = (nodes.cumulativeWrtLineNums[i] >= 0);  // Is node i a wrt node
//                                                         // (i.e. was it included
//                                                         // in the wrt argument
//                                                         // to nimDerivs?)
//     isCalcNodeLine = nodes.calcNodeIndicators[i];  // Is node i a node included
//                                                    // in the call to
//                                                    // calculate(nodes) ?

//     bool isAddedScalarNode = nodes.isAddedScalarNode[i];  
//     if(isAddedScalarNode){
//       sumAddedScalarNodes++;
//     }

//     thisWrtLine = nodes.cumulativeWrtLineNums[i];  // If node i is a wrt node,
//                                                    // which wrt node is it (i.e.
//                                                    // where does it fall in the
//                                                    // topological ordering of
//                                                    // wrt nodes) ?
//     thisNodeSize = nodes.nodeLengths[i];  // Total size (= product of lengths of
//                                           // dimensions) of this node
//     if (isCalcNodeLine | isAddedScalarNode) {
     
//       if ((nodes.cppWrtArgIndices[i][0] >
//           -1) | isAddedScalarNode) {  // -1 is used as an indicator that node i doesn't have any
//                  // parents that depend on WRT arguments, so its derivs don't
//                  // need too be taken.  Otherwise, we proceed with the chain
//                  // rule.
//                  // auto t2 = std::chrono::high_resolution_clock::now();

//         if(isAddedScalarNode){
//           thisDerivList->value[0] = 0;
//           thisDerivList->jacobian = nodes.thisAddedNodeJacobianList[0];
//           if(hessianFlag){
//             thisDerivList->hessian.initialize(0, 1, nodes.totalOutWrtSize,
//                                   nodes.totalOutWrtSize, 1);
//           }
//         }
//         else{
// #ifdef _TIME_AD
// 	  derivs_calc_timer.start();
// #endif
// 	  instructions[i - sumAddedScalarNodes].nodeFunPtr->calculateWithArgs_derivBlock(
// 											 instructions[i - sumAddedScalarNodes].operand, newDerivOrders, nodes.cppWrtArgIndices[i],
// 											 thisDerivList);  // Derivatives of calculate() for
// 	  // node i are computed here.
// #ifdef _TIME_AD
// 	  derivs_calc_timer.stop();
// #endif
//         }
//         thisRows = (*thisDerivList).jacobian.dimSize(0);
//         thisCols = (*thisDerivList).jacobian.dimSize(1);
//         vector<Map<MatrixXd, Unaligned, EigStrDyn> > thisHessian;
//         derivOutputFlag =
//             (isDeterminisitic) ? false : true;  // derivOutputFlag is true if
//                                                 // the derivative output from
//                                                 // this node will be added to
//                                                 // the ansList.
//         if (hessianFlag) {
//           chainRuleHessians[i].resize(thisNodeSize);
//           for (int j = 0; j < thisNodeSize; j++) {
//             Map<MatrixXd, Unaligned, EigStrDyn> iHessian(
//                 (*thisDerivList).hessian.getPtr() +
//                     static_cast<int>(static_cast<int>(
//                         j * (*thisDerivList).hessian.strides()[2])),
//                 thisCols, thisCols,
//                 EigStrDyn((*thisDerivList).hessian.strides()[1],
//                           (*thisDerivList).hessian.strides()[0]));
//             thisHessian.push_back(iHessian);
//             chainRuleHessians[i][j] =
//                 MatrixXd::Zero(nodes.totalWrtSize, nodes.totalWrtSize);
//           }
//         }
//         chainRuleJacobians[i] = MatrixXd::Zero(thisRows, nodes.totalWrtSize);
//         Map<MatrixXd> thisJacobian((*thisDerivList).jacobian.getPtr(), thisRows,
//                                    thisCols);
//         jLength =
//             nodes.topLevelWrtDeps[i].size();  /// also need to ensure that if
//                                               /// all wrtDeps == 0 for node i,
//                                               /// no derivs are taken or calcs
//                                               /// are performed.
//                                               //  may be taken care of by if
//                                               //  (nodes.cppWrtArgIndices[i][0]
//                                               //  > -1
//         thisArgIndex = 0;
//         if(i == 5){

//         }
//         for (int j = 0; j < jLength; j++) {
//           lineWrtSizeIJ = nodes.lineWrtArgSizeInfo[i][j];
//           if ((j == 0) & isWrtLine) {
//             if(!isDeterminisitic){
//               int wrtStartNode = nodes.wrtLineIndices[thisWrtLine][0] - 1;
//               int wrtLength = nodes.wrtLineIndices[thisWrtLine].dimSize(0);
//               if (thisRows == 1) {
//                 if (wrtLength == 1) {
//                   chainRuleJacobians[i](0, wrtStartNode) +=
//                       thisJacobian(0, thisArgIndex);
//                 } else {
//                   chainRuleJacobians[i].row(0).segment(wrtStartNode, wrtLength) +=
//                       (thisJacobian).row(0).segment(thisArgIndex, lineWrtSizeIJ);
//                 }
//               } else {
//                 chainRuleJacobians[i].block(0, wrtStartNode, thisRows,
//                                             wrtLength) +=
//                     (thisJacobian)
//                         .block(0, thisArgIndex, thisRows, lineWrtSizeIJ);
//               }
//               if (hessianFlag) {
//                 chainRuleHessians[i][0].block(wrtStartNode, wrtStartNode,
//                                               wrtLength, wrtLength) +=
//                     (thisHessian[0])
//                         .block(thisArgIndex, thisArgIndex, lineWrtSizeIJ,
//                               lineWrtSizeIJ);
//                 int thisArgIndex2 = lineWrtSizeIJ;
//                 for (int j2 = 1; j2 < jLength; j2++) {
//                   int lineWrtSizeIJ2 = nodes.lineWrtArgSizeInfo[i][j2];
//                   k2Length = nodes.topLevelWrtDeps[i][j2].size();
//                   int addToIndex2 = 0;
//                   for (int k2 = 0; k2 < k2Length; k2++) {
//                       int lLength2 = nodes.topLevelWrtDeps[i][j2][k2].dimSize(0);
//                       for (int l2 = 0; l2 < lLength2; l2++) {
//                         if (nodes.topLevelWrtDeps[i][j2][k2][l2] > 0) {
//                           int wrtNodeK2 =
//                               nodes.topLevelWrtDeps[i][j2][k2][l2] - 1;
//                           int wrtStartNode2 =
//                               nodes.wrtLineIndices[wrtNodeK2][0] - 1;
//                           int wrtLength2 =
//                               nodes.wrtLineIndices[wrtNodeK2].dimSize(0);
//                           bool transposeMatrices;
//                           if (wrtNodeK2 < thisWrtLine) {
//                             transposeMatrices = true;
//                           } else {
//                             transposeMatrices = false;
//                           }
//                           dim3Length = chainRuleHessians[i].size();
//                           for (int dim3 = 0; dim3 < dim3Length; dim3++) {
//                             if (nodes.parentIndicesList
//                                     [nodes.parentIndicesList[i][j2][k2]][0][0] ==
//                                 -1 * (wrtNodeK2 + 2)) {  // parent is wrt node,
//                                                         // mult by identity
//                               if(transposeMatrices){
//                                 chainRuleHessians[i][dim3].block(
//                                   wrtStartNode2, wrtStartNode,
//                                   wrtLength2, wrtLength) +=
//                                   thisHessian[dim3].block(
//                                       thisArgIndex, thisArgIndex2 + addToIndex2, lineWrtSizeIJ,
//                                       nodes.lineWrtArgSizeInfo[i][j2]).transpose();
//                               }
//                               else{
//                               chainRuleHessians[i][dim3].block(
//                                   wrtStartNode, wrtStartNode2,
//                                   wrtLength, wrtLength2) +=
//                                   thisHessian[dim3].block(
//                                       thisArgIndex, thisArgIndex2 + addToIndex2, lineWrtSizeIJ,
//                                       nodes.lineWrtArgSizeInfo[i][j2]);

//                               }
//                             } else {
//                               if(transposeMatrices){
//                                 chainRuleHessians[i][dim3].block(
//                                     wrtStartNode2, wrtStartNode,
//                                     wrtLength2, wrtLength) +=
//                                     chainRuleJacobians[nodes.parentIndicesList
//                                                           [i][j2][k2]]
//                                         .block(
//                                             0, wrtStartNode2,
//                                             chainRuleJacobians
//                                                 [nodes.parentIndicesList[i][j2][k2]]
//                                                     .rows(),
//                                             wrtLength2).transpose() * 
//                                             thisHessian[dim3].block(
//                                         thisArgIndex, thisArgIndex2 + addToIndex2, lineWrtSizeIJ,
//                                         nodes.lineWrtArgSizeInfo[i][j2]).transpose();
//                               }
//                               else{
//                                 chainRuleHessians[i][dim3].block(
//                                     wrtStartNode, wrtStartNode2,
//                                     wrtLength, wrtLength2) +=
//                                     thisHessian[dim3].block(
//                                         thisArgIndex, thisArgIndex2 + addToIndex2, lineWrtSizeIJ,
//                                         nodes.lineWrtArgSizeInfo[i][j2]) *
//                                     chainRuleJacobians[nodes.parentIndicesList
//                                                           [i][j2][k2]]
//                                         .block(
//                                             0, wrtStartNode2,
//                                             chainRuleJacobians
//                                                 [nodes.parentIndicesList[i][j2][k2]]
//                                                     .rows(),
//                                             wrtLength2);
//                               }
//                             }
//                           }
//                         }
//                     }
//                                                      addToIndex2 += nodes.nodeLengths[nodes.parentIndicesList[i][j2][k2]];;

//                   }
//                   thisArgIndex2 += lineWrtSizeIJ2;
//                 }
//               }
//             }
//           } else {
//             kLength = nodes.topLevelWrtDeps[i][j].size();
//             int addToIndex = 0;
//             for (int k = 0; k < kLength; k++) {
//                 lLength = nodes.topLevelWrtDeps[i][j][k].dimSize(0);
//                 int parentRowLength =
//                     chainRuleJacobians[nodes.parentIndicesList[i][j][k]].rows();
//                 for (int l = 0; l < lLength; l++) {
//                   if (nodes.topLevelWrtDeps[i][j][k][l] > 0) {
//                     wrtNodeK = nodes.topLevelWrtDeps[i][j][k][l] - 1;
//                     int wrtStartNode = nodes.wrtLineIndices[wrtNodeK][0] - 1;
//                     int wrtLength = nodes.wrtLineIndices[wrtNodeK].dimSize(0);
                    
//                     if (nodes
//                             .parentIndicesList[nodes.parentIndicesList[i][j][k]]
//                                               [0][0] ==
//                         -1 * (wrtNodeK +
//                               2)) {  // parent is wrt node, mult by identity
//                       chainRuleJacobians[i]
//                           .block(  // this doesn't always need to be saved
//                               0, wrtStartNode, thisRows, wrtLength) +=
//                           (thisJacobian)
//                               .block(0, thisArgIndex + addToIndex, thisRows, wrtLength);
//                     } else {
//                       chainRuleJacobians[i]
//                           .block(  // this doesn't always need to be saved
//                               0, wrtStartNode, thisRows, wrtLength) +=
//                           (thisJacobian)
//                               .block(0, thisArgIndex + addToIndex, thisRows, parentRowLength) *
//                           chainRuleJacobians[nodes.parentIndicesList[i][j][k]]
//                               .block(0, wrtStartNode, parentRowLength,
//                                      wrtLength);
//                     }
//                     if (hessianFlag) {
//                       int thisArgIndex2 = nodes.lineWrtArgSizeInfo[i][0];
//                       for (int j2 = 1; j2 <= j; j2++) { // j and j2 iterate over arguments to node i's calculate function
//                         int lineWrtSizeIJ2 = nodes.lineWrtArgSizeInfo[i][j2];
//                         int k2Length = nodes.topLevelWrtDeps[i][j2].size();
//                         int addToIndex2 = 0;
//                         for (int k2 = 0; k2 < k2Length; k2++) { // k and k2 iterate over parent nodes to argument j and j2 
//                             int lLength2 =
//                                 nodes.topLevelWrtDeps[i][j2][k2].dimSize(0);
//                             if((j2 == j) & (k2 == k)){
//                                 lLength2 = l + 1; 
//                             }
//                             for (int l2 = 0; l2 < lLength2; l2++) { // l and l2 iterate over top level wrt nodes that the parent nodes (k and k2) depend on
//                               if(nodes.topLevelWrtDeps[i][j2][k2][l2] > 0) {
//                                 int wrtNodeK2 =
//                                     nodes.topLevelWrtDeps[i][j2][k2][l2] - 1;
//                                 int wrtStartNode2 =
//                                     nodes.wrtLineIndices[wrtNodeK2][0] - 1;
//                                 int wrtLength2 =
//                                     nodes.wrtLineIndices[wrtNodeK2].dimSize(0);
//                                 bool transposeMatrices;
//                                 if (wrtNodeK2 < wrtNodeK) {
//                                   transposeMatrices = true;
//                                 } else {
//                                   transposeMatrices = false;
//                                 }
//                                 dim3Length = chainRuleHessians[i].size();
//                                 if(!nodes.isAddedScalarNode[nodes.parentIndicesList[i][j][k]]){
//                                   if (nodes.parentIndicesList
//                                           [nodes.parentIndicesList[i][j][k]][0]
//                                           [0] != -1 * (wrtNodeK + 2)) {
//                                     for (int m1 = 0; m1 < wrtLength; m1++) {
//                                       for (int m2 = 0; m2 < wrtLength2;
//                                           m2++) {
//                                         VectorXd addToHessian(parentRowLength);
//                                         for (int m3 = 0; m3 < parentRowLength;
//                                             m3++) {
//                                           addToHessian[m3] = chainRuleHessians
//                                               [nodes.parentIndicesList[i][j][k]]
//                                               [m3](nodes.wrtLineIndices
//                                                           [wrtNodeK][m1] -
//                                                       1,
//                                                   nodes.wrtLineIndices
//                                                           [wrtNodeK2][m2] -
//                                                       1);
//                                       }
//                                         addToHessian =
//                                             (thisJacobian)
//                                                 .block(0, thisArgIndex + addToIndex,
//                                                       thisRows,
//                                                       parentRowLength) *
//                                             addToHessian;

//                                         for (int m3 = 0; m3 < dim3Length; m3++) {
//                                           if(transposeMatrices){
//                                           chainRuleHessians[i][m3](
//                                               nodes.wrtLineIndices[wrtNodeK2]
//                                                                   [m2] -
//                                                   1,
//                                               nodes.wrtLineIndices[wrtNodeK]
//                                                                   [m1] -
//                                                   1) += addToHessian[m3];
//                                                 }
//                                                 else{
//                                                   chainRuleHessians[i][m3](
//                                               nodes.wrtLineIndices[wrtNodeK]
//                                                                   [m1] -
//                                                   1,
//                                               nodes.wrtLineIndices[wrtNodeK2]
//                                                                   [m2] -
//                                                   1) += addToHessian[m3];
//                                                 }
//                                         }
//                                       }
//                                     }
//                                   }
//                                 }
//                                 if(!isAddedScalarNode){
//                                   for (int dim3 = 0; dim3 < dim3Length; dim3++) {
//                                     if (nodes.parentIndicesList
//                                             [nodes.parentIndicesList[i][j][k]][0]
//                                             [0] ==
//                                         -1 * (wrtNodeK + 2)) {  // parent is wrt
//                                                                 // node, mult by
//                                                                 // identity
//                                       if (nodes.parentIndicesList
//                                               [nodes.parentIndicesList[i][j2][k2]]
//                                               [0][0] ==
//                                           -1 * (wrtNodeK2 +
//                                                 2)) {  // parent is wrt node, mult
//                                                       // by identity
//                                         if(transposeMatrices){
//                                           chainRuleHessians[i][dim3].block(
//                                             wrtStartNode2, wrtStartNode,
//                                             wrtLength2, wrtLength) +=
//                                             thisHessian[dim3].block(
//                                                 thisArgIndex + addToIndex, thisArgIndex2 + addToIndex2,
//                                                 nodes.lineWrtArgSizeInfo[i][j],
//                                                 nodes.lineWrtArgSizeInfo[i][j2]).transpose();
//                                         }
//                                         else{
//                                         chainRuleHessians[i][dim3].block(
//                                             wrtStartNode, wrtStartNode2,
//                                             wrtLength, wrtLength2) +=
//                                             thisHessian[dim3].block(
//                                                 thisArgIndex + addToIndex, thisArgIndex2 + addToIndex2,
//                                                 nodes.lineWrtArgSizeInfo[i][j],
//                                                 nodes.lineWrtArgSizeInfo[i][j2]);
//                                         }
//                                       } else {
//                                         if(transposeMatrices){
//                                         chainRuleHessians[i][dim3].block(
//                                             wrtStartNode2, wrtStartNode,
//                                             wrtLength2, wrtLength) +=
//                                             chainRuleJacobians
//                                                 [nodes.parentIndicesList[i][j2]
//                                                                         [k2]]
//                                                     .block(
//                                                         0, wrtStartNode2,
//                                                         chainRuleJacobians
//                                                             [nodes
//                                                                 .parentIndicesList
//                                                                     [i][j2][k2]]
//                                                                 .rows(),
//                                                         wrtLength2).transpose() *
//                                             thisHessian[dim3].block(
//                                                 thisArgIndex + addToIndex, thisArgIndex2 + addToIndex2,
//                                                 nodes.lineWrtArgSizeInfo[i][j],
//                                                 nodes.lineWrtArgSizeInfo[i][j2]).transpose();
//                                         }
//                                         else{
//                                           chainRuleHessians[i][dim3].block(
//                                             wrtStartNode, wrtStartNode2,
//                                             wrtLength, wrtLength2) +=
//                                             thisHessian[dim3].block(
//                                                 thisArgIndex + addToIndex, thisArgIndex2 + addToIndex2,
//                                                 nodes.lineWrtArgSizeInfo[i][j],
//                                                 nodes.lineWrtArgSizeInfo[i][j2]) *
//                                             chainRuleJacobians
//                                                 [nodes.parentIndicesList[i][j2]
//                                                                         [k2]]
//                                                     .block(
//                                                         0, wrtStartNode2,
//                                                         chainRuleJacobians
//                                                             [nodes
//                                                                 .parentIndicesList
//                                                                     [i][j2][k2]]
//                                                                 .rows(),
//                                                         wrtLength2);

//                                         }            
//                                       }
//                                     } else if (nodes.parentIndicesList
//                                                   [nodes.parentIndicesList[i][j2]
//                                                                           [k2]]
//                                                   [0][0] ==
//                                               -1 * (wrtNodeK2 + 2)) {
//                                     if(transposeMatrices){
//                                       chainRuleHessians[i][dim3].block(
//                                           wrtStartNode2, wrtStartNode, wrtLength2,
//                                           wrtLength) +=
//                                           thisHessian[dim3].block(
//                                               thisArgIndex + addToIndex, thisArgIndex2 + addToIndex2,
//                                               nodes.lineWrtArgSizeInfo[i][j],
//                                               nodes.lineWrtArgSizeInfo[i][j2]).transpose() *
//                                               chainRuleJacobians
//                                               [nodes.parentIndicesList[i][j][k]]
//                                                   .block(
//                                                       0, wrtStartNode,
//                                                       chainRuleJacobians
//                                                           [nodes.parentIndicesList
//                                                               [i][j][k]]
//                                                               .rows(),
//                                                       wrtLength);
//                                             }
//                                             else{
//                                               chainRuleHessians[i][dim3].block(
//                                           wrtStartNode, wrtStartNode2, wrtLength,
//                                           wrtLength2) +=
//                                           chainRuleJacobians
//                                               [nodes.parentIndicesList[i][j][k]]
//                                                   .block(
//                                                       0, wrtStartNode,
//                                                       chainRuleJacobians
//                                                           [nodes.parentIndicesList
//                                                               [i][j][k]]
//                                                               .rows(),
//                                                       wrtLength)
//                                                   .transpose() *
//                                           thisHessian[dim3].block(
//                                               thisArgIndex + addToIndex, thisArgIndex2 + addToIndex2,
//                                               nodes.lineWrtArgSizeInfo[i][j],
//                                               nodes.lineWrtArgSizeInfo[i][j2]);
//                                             }
//                                     } else {
//                                       if(transposeMatrices){
//                                       chainRuleHessians[i][dim3].block(
//                                           wrtStartNode2, wrtStartNode, wrtLength2,
//                                           wrtLength) +=
//                                           chainRuleJacobians
//                                               [nodes.parentIndicesList[i][j2][k2]]
//                                                   .block(
//                                                       0, wrtStartNode2,
//                                                       chainRuleJacobians
//                                                           [nodes.parentIndicesList
//                                                               [i][j2][k2]]
//                                                               .rows(),
//                                                       wrtLength2).transpose() *
//                                           thisHessian[dim3].block(
//                                               thisArgIndex + addToIndex, thisArgIndex2 + addToIndex2,
//                                               nodes.lineWrtArgSizeInfo[i][j],
//                                               nodes.lineWrtArgSizeInfo[i][j2]) *
//                                           chainRuleJacobians
//                                               [nodes.parentIndicesList[i][j][k]]
//                                                   .block(
//                                                       0, wrtStartNode,
//                                                       chainRuleJacobians
//                                                           [nodes.parentIndicesList
//                                                               [i][j][k]]
//                                                               .rows(),
//                                                       wrtLength);
//                                                     }
//                                                     else{
//                                                     chainRuleHessians[i][dim3].block(
//                                           wrtStartNode, wrtStartNode2, wrtLength,
//                                           wrtLength2) +=
//                                           chainRuleJacobians
//                                               [nodes.parentIndicesList[i][j][k]]
//                                                   .block(
//                                                       0, wrtStartNode,
//                                                       chainRuleJacobians
//                                                           [nodes.parentIndicesList
//                                                               [i][j][k]]
//                                                               .rows(),
//                                                       wrtLength)
//                                                   .transpose() *
//                                           thisHessian[dim3].block(
//                                               thisArgIndex + addToIndex, thisArgIndex2 + addToIndex2,
//                                               nodes.lineWrtArgSizeInfo[i][j],
//                                               nodes.lineWrtArgSizeInfo[i][j2])   *
//                                           chainRuleJacobians
//                                               [nodes.parentIndicesList[i][j2][k2]]
//                                                   .block(
//                                                       0, wrtStartNode2,
//                                                       chainRuleJacobians
//                                                           [nodes.parentIndicesList
//                                                               [i][j2][k2]]
//                                                               .rows(),
//                                                       wrtLength2); 
//                                                     }
//                                     }
//                                   }
//                                 }
//                             }
//                         }
//                         addToIndex2 += nodes.nodeLengths[nodes.parentIndicesList[i][j2][k2]];;
//                         }
//                         thisArgIndex2 += lineWrtSizeIJ2;
//                       }
//                     }
//                   }
//                 }
//                 addToIndex += nodes.nodeLengths[nodes.parentIndicesList[i][j][k]];
//             }
//           }
//           thisArgIndex += lineWrtSizeIJ;
//         }
//         if (derivOutputFlag) {  // Add this Jacobian to the output of this
//           jLength = nodes.allNeededWRTCopyVars[i].dimSize(0);
//           for (int j = 0; j < jLength; j++) {
//             wrtNodeJ = nodes.allNeededWRTCopyVars[i][j] - 1;
//             int wrtToStartNode = nodes.wrtToIndices[wrtNodeJ][0] - 1;
//             int wrtToLength = nodes.wrtToIndices[wrtNodeJ].dimSize(0);
//             int wrtFromStartNode = nodes.wrtFromIndices[wrtNodeJ][0] - 1;
//             int wrtFromLength = nodes.wrtFromIndices[wrtNodeJ].dimSize(0);
//             ansJacobian.row(0).segment(wrtToStartNode, wrtToLength) +=
//                 chainRuleJacobians[i].row(0).segment(wrtFromStartNode,
//                                                      wrtFromLength);
//             if(hessianFlag){
//               for (int j2 = j; j2 < jLength; j2++) {
//                 int wrtNodeJ2 = nodes.allNeededWRTCopyVars[i][j2] - 1;
//                 int wrtToStartNode2 = nodes.wrtToIndices[wrtNodeJ2][0] - 1;
//                 int wrtToLength2 = nodes.wrtToIndices[wrtNodeJ2].dimSize(0);
//                 int wrtFromStartNode2 = nodes.wrtFromIndices[wrtNodeJ2][0] - 1;
//                 int wrtFromLength2 = nodes.wrtFromIndices[wrtNodeJ2].dimSize(0);
//                 ansHessian.block(wrtToStartNode, wrtToStartNode2, wrtToLength,
//                                 wrtToLength2) +=
//                     chainRuleHessians[i][0].block(wrtFromStartNode,
//                                                   wrtFromStartNode2,
//                                                   wrtFromLength, wrtFromLength2);
//               }
//             }
//           }
//         }
//       } else {  // Otherwise, no arguments depend on wrt nodes, so we don't need
//                 // to take derivatives. Instead, just set Jacobian and Hessian
//                 // to 0.
//         if (valueFlag) {  // Still may need to get the value (0'th order deriv).
//           NimArr<1, double> valueOrder(1);
//           valueOrder[0] = 0;
//           instructions[i].nodeFunPtr->calculateWithArgs_derivBlock(
//               instructions[i].operand, valueOrder, nodes.cppWrtArgIndices[i],
//               thisDerivList);
//         }
//       }
//     }
//     if (!isDeterminisitic) {  
//       if (isCalcNodeLine && valueFlag) {
//         (*ansList).value[0] = (*ansList).value[0] + (*thisDerivList).value[0];
//       }
//     }
// #ifdef _TIME_AD
//     derivs_node_iteration_timer.stop();	
// #endif
//   }
//   if (hessianFlag) {  // reflect Hessian across the diagonal
//     ansHessian.triangularView<Eigen::Lower>() =
//         ansHessian.transpose().triangularView<Eigen::Lower>();
//   }
// #ifdef _TIME_AD
//   derivs_main_timer.stop();
//   derivs_main_timer.report();
//   derivs_node_iteration_timer.report();
//   derivs_calc_timer.report();
//   derivs_getDerivs_timer.report();
//   derivs_run_tape_timer.report();
// #endif
//   return (ansList);
// }

#endif

nimSmartPtr<NIMBLE_ADCLASS> NIM_DERIVS_CALCULATE(
    NodeVectorClassNew_derivs &nodes, const double derivOrders) {
  NimArr<1, double> orders(1);
  orders[0] = derivOrders;
  return (NIM_DERIVS_CALCULATE(nodes, orders));
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

void nimbleFunctionCppADbase::getDerivs(nimbleCppADinfoClass &ADinfo,
                                        NimArr<1, double> &derivOrders,
                                        const NimArr<1, double> &wrtVector,
                                        nimSmartPtr<NIMBLE_ADCLASS> &ansList) {
#ifdef _TIME_AD
  derivs_getDerivs_timer_start();
  derivs_tick_id();
  derivs_show_id();  
#endif
  std::size_t n = ADinfo.independentVars.size();  // dim of independent vars

  std::size_t wrt_n = wrtVector.size();            // dim of wrt vars
  if(wrt_n == 2){
    if(wrtVector[1] == -1){
      wrt_n = 1;
    }
  }
  int orderSize = derivOrders.size();
  double array_derivOrders[orderSize];

  std::memcpy(array_derivOrders, derivOrders.getPtr(),
	      orderSize * sizeof(double));

  int maxOrder =
    *std::max_element(array_derivOrders, array_derivOrders + orderSize);
  bool ordersFound[3] = {false};

  for (int i = 0; i < orderSize; i++) {
    if ((array_derivOrders[i] > 2) | (array_derivOrders[i] < 0)) {
      printf("Error: Derivative orders must be between 0 and 2.\n");
    }
    ordersFound[static_cast<int>(array_derivOrders[i])] = true;
  }
  vector<double> value_ans;
#ifdef _TIME_AD
  derivs_run_tape_timer_start();
#ifdef _SHOW_NODE_BY_NODE 
  std::cout<<"Running value "<<std::endl;
#endif
#endif
  value_ans = ADinfo.ADtape->Forward(0, ADinfo.independentVars);
#ifdef _TIME_AD
  derivs_run_tape_timer_stop();
#endif
  if (ordersFound[0] == true) {
    ansList->value = vectorDouble_2_NimArr(value_ans);
  }
  if(maxOrder > 0){
    std::size_t q = value_ans.size();
    vector<bool> infIndicators(q); 
    for(size_t inf_ind = 0; inf_ind < q; inf_ind++){
      if(((value_ans[inf_ind] == -std::numeric_limits<double>::infinity()) |
          (value_ans[inf_ind] == std::numeric_limits<double>::infinity())) | 
	 (isnan(value_ans[inf_ind]))){
	infIndicators[inf_ind] = true;
      }
      else{
	infIndicators[inf_ind] = false;
      }
    }
    if (ordersFound[1] == true) {
      ansList->jacobian.setSize(q, wrt_n, false, false); // setSize may be costly.  Possible to setSize outside of fxn, within chain rule algo, and only resize when necessary?
    }
    if (ordersFound[2] == true) {
      ansList->hessian.setSize(wrt_n, wrt_n, q, false, false);
    }
    vector<double> cppad_derivOut;
    for (size_t dy_ind = 0; dy_ind < q; dy_ind++) {
      std::vector<double> w(q, 0);
      w[dy_ind] = 1;
      if (maxOrder == 1) {   
	if(infIndicators[dy_ind] == false){
#ifdef _TIME_AD
#ifdef _SHOW_NODE_BY_NODE 
	  std::cout<<"Running deriv "<<dy_ind<<std::endl;
#endif
	  derivs_run_tape_timer_start();
#endif
	  cppad_derivOut = ADinfo.ADtape->Reverse(1, w);
#ifdef _TIME_AD
	  derivs_run_tape_timer_stop();
#endif
	}
      } else {
	for (size_t vec_ind = 0; vec_ind < wrt_n; vec_ind++) {
	  if(infIndicators[dy_ind] == false){
	    int dx1_ind = wrtVector[vec_ind] - 1;
	    std::vector<double> x1(n, 0);  // vector specifying first derivatives.
	    // first specify coeffs for first dim
	    // of s across all directions r, then
	    // second dim, ...
	    x1[dx1_ind] = 1;
#ifdef _TIME_AD
	    derivs_run_tape_timer_start();
#endif
	    ADinfo.ADtape->Forward(1, x1);
	    cppad_derivOut = ADinfo.ADtape->Reverse(2, w);
#ifdef _TIME_AD
	    derivs_run_tape_timer_stop();
#endif
	  }
	  for (size_t vec_ind2 = 0; vec_ind2 < wrt_n; vec_ind2++) {
	    if(infIndicators[dy_ind] == false){
	      int dx2_ind = wrtVector[vec_ind2] - 1;
	      ansList->hessian[wrt_n * wrt_n * dy_ind + wrt_n * vec_ind + vec_ind2] =
		cppad_derivOut[dx2_ind * 2 + 1];
	    }
	    else{
	      ansList->hessian[wrt_n * wrt_n * dy_ind + wrt_n * vec_ind + vec_ind2] = 
		CppAD::numeric_limits<double>::quiet_NaN();
	    }
	  }
	}
      }
      if (ordersFound[1] == true) {
	for (size_t vec_ind = 0; vec_ind < wrt_n; vec_ind++) {
	  if(infIndicators[dy_ind] == false){
	    int dx1_ind = wrtVector[vec_ind] - 1;
	    ansList->jacobian[vec_ind * q + dy_ind] =
	      cppad_derivOut[dx1_ind * maxOrder + 0];
	  }
	  else{
	    ansList->jacobian[vec_ind * q + dy_ind] =
	      CppAD::numeric_limits<double>::quiet_NaN();
	  }     
	}
      }
    }
  }
#ifdef _TIME_AD
  derivs_getDerivs_timer_stop();
#endif
}
