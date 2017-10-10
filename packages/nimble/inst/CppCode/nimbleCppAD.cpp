#include <nimble/nimbleCppAD.h>
//#include <nimble/accessorClasses.h>

nimSmartPtr<NIMBLE_ADCLASS>  NIM_DERIVS_CALCULATE(NodeVectorClassNew &nodes, NimArr<1, double> &derivOrders) {
  nimSmartPtr<NIMBLE_ADCLASS> ADlist = new NIMBLE_ADCLASS;
  return(ADlist);
  /* const vector<NodeInstruction> &instructions = nodes.getInstructions();
  vector<NodeInstruction>::const_iterator iNode(instructions.begin());
  vector<NodeInstruction>::const_iterator iNodeEnd(instructions.end());
  for(; iNode != iNodeEnd; iNode++)
    ans += iNode->nodeFunPtr->calculateBlock(iNode->operand);
  return(ans); */
}


nimSmartPtr<NIMBLE_ADCLASS>  NIM_DERIVS_CALCULATE(NodeVectorClassNew &nodes,  int iNodeFunction, NimArr<1, double> &derivOrders) {
  nimSmartPtr<NIMBLE_ADCLASS> ADlist = new NIMBLE_ADCLASS;
  return(ADlist);
  /* const vector<NodeInstruction> &instructions = nodes.getInstructions();
  vector<NodeInstruction>::const_iterator iNode(instructions.begin());
  vector<NodeInstruction>::const_iterator iNodeEnd(instructions.end());
  for(; iNode != iNodeEnd; iNode++)
    ans += iNode->nodeFunPtr->calculateBlock(iNode->operand);
  return(ans); */
}

