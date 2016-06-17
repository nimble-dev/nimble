#include "nimble/nimbleGraph.h"

graphNode::graphNode(int inputCgraphID, NODETYPE inputType, const string &inputName ) :
  role(UNKNOWNROLE),
  type(inputType),
  CgraphID(inputCgraphID),
  name(inputName),
  touched(false),
  numChildren(0) {
  RgraphID = CgraphID + 1;
};

void graphNode::addChild( graphNode *toNode, int childParentExpressionID) {
#ifdef _DEBUGNIMGRAPH
  PRINTF("adding child %s to parent %s\n", toNode->name.c_str(), name.c_str());
#endif
  children.push_back(toNode);
  childrenParentExpressionIDs.push_back(childParentExpressionID);
  numChildren++;
  toNode->addParent(this);
};

void graphNode::addParent(graphNode *fromNode) {
#ifdef _DEBUGNIMGRAPH
  PRINTF("adding parent %s to child %s\n", fromNode->name.c_str(), name.c_str());
#endif
  parents.push_back(fromNode);
};

void SEXP_2_nodeType(SEXP Stypes, vector<NODETYPE> &ans) {
  //  enum NODETYPE {UNKNOWNTYPE, STOCH, DETERM, RHSONLY};
  if(!isString(Stypes)) {
    PRINTF("Error:  called for SEXP that is not a string!\n"); 
    return;
  }
  int nn = LENGTH(Stypes);
  ans.resize(nn);
  string oneString;
  for(int i = 0; i < nn; i++) {
    oneString.assign(CHAR(STRING_ELT(Stypes, i)), LENGTH(STRING_ELT(Stypes, i)));
    if(oneString == "stoch")
      ans[i] = STOCH;
    else if(oneString == "determ")
      ans[i] = DETERM;
    else if(oneString == "RHSonly")
      ans[i] = RHSONLY;
    else if(oneString == "LHSinferred")
      ans[i] = LHSINFERRED;
    else if(oneString == "unknown")
      ans[i] = UNKNOWNTYPE;
    else {
      ans[i] = UNKNOWNTYPE;
      PRINTF("In SEXP_2_nodeType: unknown string type label %s\n", oneString.c_str());
    }
  }
}

SEXP setGraph(SEXP SedgesFrom, SEXP SedgesTo, SEXP SedgesFrom2ParentExprIDs, SEXP Stypes, SEXP Snames, SEXP SnumNodes) {
  vector<int> edgesFrom = SEXP_2_vectorInt(SedgesFrom, -1); // -1 subtracted here
  vector<int> edgesTo = SEXP_2_vectorInt(SedgesTo, -1); // -1 substracted here
  vector<int> edgesFrom2ParentExprIDs = SEXP_2_vectorInt(SedgesFrom2ParentExprIDs);
  vector<NODETYPE> types;
  SEXP_2_nodeType(Stypes, types);
  vector<string> names;
  STRSEXP_2_vectorString(Snames, names);
  int numNodes = SEXP_2_int(SnumNodes);
  nimbleGraph *newGraph = new nimbleGraph;
  newGraph->setNodes(edgesFrom, edgesTo, edgesFrom2ParentExprIDs, types, names, numNodes);
  SEXP SextPtrAns;
  PROTECT(SextPtrAns = R_MakeExternalPtr(newGraph, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(SextPtrAns, &nimbleGraphFinalizer, TRUE);
  UNPROTECT(1);
  return(SextPtrAns);
};

void nimbleGraphFinalizer(SEXP SgraphExtPtr) {
  nimbleGraph *graphPtr = static_cast<nimbleGraph *>(R_ExternalPtrAddr(SgraphExtPtr));
  delete graphPtr;
}

nimbleGraph::~nimbleGraph() {
  int n = graphNodeVec.size();
  for(int i = 0; i < n; i++) {
    delete graphNodeVec[i];
  }
}

SEXP anyStochDependencies(SEXP SgraphExtPtr) {
  nimbleGraph *graphPtr = static_cast<nimbleGraph *>(R_ExternalPtrAddr(SgraphExtPtr)); 
  vector<int> ans(graphPtr->anyStochDependencies());
  SEXP Sans;
  PROTECT(Sans = allocVector(LGLSXP, ans.size()));
  int *SansPtr = INTEGER(Sans);
  for(int i = 0; i < ans.size(); i++) {
    if(ans[i] == 0) PRINTF("Element %i was not processed\n", i);
    SansPtr[i] = ans[i]==2 ? 1 : 0;
  }
  UNPROTECT(1);
  return(Sans);
}

SEXP getDependencies(SEXP SgraphExtPtr, SEXP Snodes, SEXP Somit, SEXP Sdownstream) {
  nimbleGraph *graphPtr = static_cast<nimbleGraph *>(R_ExternalPtrAddr(SgraphExtPtr)); 
  vector<int> nodes = SEXP_2_vectorInt(Snodes, -1); // subtract 1 index for C
  vector<int> omit = SEXP_2_vectorInt(Somit, -1);
  bool downstream = SEXP_2_bool(Sdownstream);
  vector<int> ans = graphPtr->getDependencies(nodes, omit, downstream);
  return(vectorInt_2_SEXP(ans, 1)); // add 1 index for R
}

SEXP anyStochParents(SEXP SgraphExtPtr) {
  nimbleGraph *graphPtr = static_cast<nimbleGraph *>(R_ExternalPtrAddr(SgraphExtPtr)); 
  vector<int> ans(graphPtr->anyStochParents());
  SEXP Sans;
  PROTECT(Sans = allocVector(LGLSXP, ans.size()));
  int *SansPtr = INTEGER(Sans);
  for(int i = 0; i < ans.size(); i++) {
    if(ans[i] == 0) PRINTF("Element %i was not processed\n", i);
    SansPtr[i] = ans[i]==2 ? 1 : 0;
  }
  UNPROTECT(1);
  return(Sans);
}

void nimbleGraph::setNodes(const vector<int> &edgesFrom, const vector<int> &edgesTo,
		     const vector<int> &edgesFrom2ParentExprIDs,
		     const vector<NODETYPE> &types,
		     const vector<string> &names,
		     int inputNumNodes) {
  numNodes = inputNumNodes;
  int numEdges = edgesFrom.size();

#ifdef _DEBUGNIMGRAPH
  PRINTF("numNodes %i\n", numNodes);
  PRINTF("numEdges %i\n", numEdges);
#endif
  if(numEdges != edgesTo.size() | numEdges != edgesFrom2ParentExprIDs.size() | numNodes != types.size() | numNodes != names.size()) {
    PRINTF("Something is not the right size\n");
    return;
  }
  graphNodeVec.resize(numNodes);
  for(int iNode = 0; iNode < numNodes; iNode++) {
    graphNodeVec[iNode] = new graphNode(iNode, types[iNode], names[iNode]);
  }
  for(int iEdge = 0; iEdge < numEdges; iEdge++) {
    graphNodeVec[ edgesFrom[iEdge]]->addChild( graphNodeVec[edgesTo[iEdge]], edgesFrom2ParentExprIDs[iEdge] );
  }
}

vector<int> nimbleGraph::anyStochDependencies() {
  vector<int> ans(numNodes, 0);
  for(int i = 0; i < numNodes; i++) {
    anyStochDependenciesOneNode(ans, i); 
  }
  return(ans);
}

bool nimbleGraph::anyStochDependenciesOneNode(vector<int> &anyStochDependencies, int CgraphID) {
  // 0 = untouched, 1 = false, 2 = true
  if(anyStochDependencies[CgraphID] != 0) return(anyStochDependencies[CgraphID] == 2);
  bool thisHasAstochDep(false);
  graphNode *thisGraphNode = graphNodeVec[CgraphID];
  graphNode *thisChildNode;
  int numChildren = thisGraphNode->numChildren;
  /* If no children, answer is false */
  if(numChildren == 0) {
    anyStochDependencies[CgraphID] = 1;
    return(false);
  }
  int i(0);
  /* Check type of children without recursing.  If any are STOCH, answer is true */
  while((i < numChildren) & (!thisHasAstochDep)) {
    if(thisGraphNode->children[i]->type == STOCH) {
      thisHasAstochDep = true;
    }
    i++;
  }
  /* If answer was true, we're done */
  if(thisHasAstochDep) {
    anyStochDependencies[CgraphID] = 2;
    return(true);
  }
  /* all children were not STOCH, so now recurse through children */
  i = 0;
  while((i < numChildren) & (!thisHasAstochDep)) {
    thisChildNode = thisGraphNode->children[i];
    if(anyStochDependenciesOneNode(anyStochDependencies, thisChildNode->CgraphID)) {
      thisHasAstochDep = true;
    }
    i++;
  }
  if(thisHasAstochDep) {
    anyStochDependencies[CgraphID] = 2;
    return(true);
  }
  anyStochDependencies[CgraphID] = 1;
  return(false);
}


vector<int> nimbleGraph::anyStochParents() {
  vector<int> ans(numNodes, 0);
  for(int i = static_cast<int>(numNodes - 1); i >= 0; i--) {
    anyStochParentsOneNode(ans, i);
  }
  return(ans);
}

bool nimbleGraph::anyStochParentsOneNode(vector<int> &anyStochParents, int CgraphID) {
  // 0 = untouched, 1 = false, 2 = true
  if(anyStochParents[CgraphID] != 0) return(anyStochParents[CgraphID] == 2);
  bool thisHasAstochParent(false);
  graphNode *thisGraphNode = graphNodeVec[CgraphID];
  graphNode *thisParentNode;
  int numParents = thisGraphNode->parents.size();
  if(numParents == 0) {
    anyStochParents[CgraphID] = 1;
    return(false);
  }
  int i(0);
  while((i < numParents) & (!thisHasAstochParent)) {
    if(thisGraphNode->parents[i]->type == STOCH) {
      thisHasAstochParent = true;
    }
    i++;
  }
  if(thisHasAstochParent) {
    anyStochParents[CgraphID] = 2;
    return(true);
  }
  i = 0;
  while((i < numParents) & (!thisHasAstochParent)) {
    thisParentNode = thisGraphNode->parents[i];
    if(anyStochParentsOneNode(anyStochParents, thisParentNode->CgraphID)) {
      thisHasAstochParent = true;
    }
    i++;
  }
  if(thisHasAstochParent) {
    anyStochParents[CgraphID] = 2;
    return(true);
  }
  anyStochParents[CgraphID] = 1;
  return(false);
}

//#define _DEBUG_GETDEPS

vector<int> nimbleGraph::getDependencies(const vector<int> &Cnodes, const vector<int> &Comit, bool downstream) {
  // assume on entry that touched = false on all nodes
  // Cnodes and Comit are C-indices (meaning they start at 0)
  int n = Comit.size();
  int i;
  vector<int> ans;
  // touch omit nodes
#ifdef _DEBUG_GETDEPS
  int iDownstream = static_cast<int>(downstream);
  PRINTF("debugging output for getDependencies with %i nodes, %i omits, and downstream = %i.  C indices (graphIDs) shown are 0-based\n", Cnodes.size(), Comit.size(), iDownstream);
#endif
  for(i = 0; i < n; i++) {
    graphNodeVec[ Comit[i] ]->touched = true;
#ifdef _DEBUG_GETDEPS
    PRINTF("touching %i to omit\n", Comit[i]);
#endif
  }
  n = Cnodes.size();
  graphNode *thisGraphNode;
  int thisGraphNodeID;
  for(i = 0; i < n; i++) {
    thisGraphNodeID = Cnodes[i];
    thisGraphNode = graphNodeVec[ thisGraphNodeID ];
#ifdef _DEBUG_GETDEPS
    PRINTF("Working on input node %i\n", thisGraphNodeID);
#endif
    if(!thisGraphNode->touched) {
#ifdef _DEBUG_GETDEPS
      PRINTF("  Adding node %i to ans and recursing\n", thisGraphNodeID);
#endif
      ans.push_back(thisGraphNodeID);
      thisGraphNode->touched = true;
      getDependenciesOneNode(ans, thisGraphNodeID, downstream, 1);      
    } else {
#ifdef _DEBUG_GETDEPS
      PRINTF("  Node %i was already touched\n", thisGraphNodeID);
#endif
    }
  }

  // untouch nodes and omit
  n = Comit.size();
  for(i = 0; i < n; i++) {
    graphNodeVec[ Comit[i] ]->touched = false;
  }
  n = ans.size();
  for(i = 0; i < n; i++) {
    graphNodeVec[ ans[i] ]->touched = false;
  }
  std::sort(ans.begin(), ans.end());
  return(ans);
}

void nimbleGraph::getDependenciesOneNode(vector<int> &deps, int CgraphID, bool downstream, int recursionDepth) {
  if(recursionDepth > graphNodeVec.size()) {
    PRINTF("ERROR: getDependencies has recursed too far.  Something must be wrong.\n");
    return;
  }
#ifdef _DEBUG_GETDEPS
  PRINTF("    Entering recursion for node %i\n", CgraphID);
#endif
  graphNode *thisGraphNode = graphNodeVec[CgraphID];
  int numChildren = thisGraphNode->numChildren;
  int i(0);
  graphNode *thisChildNode;
  int thisChildCgraphID;
#ifdef _DEBUG_GETDEPS
  PRINTF("      Starting to iterate through %i children of node %i\n", numChildren, CgraphID);
#endif
  for(; i < numChildren; i++) {
    thisChildNode = thisGraphNode->children[i];
    if(thisChildNode->touched) continue;
    thisChildCgraphID = thisChildNode->CgraphID;
#ifdef _DEBUG_GETDEPS
    PRINTF("        Adding child node %i\n", thisChildCgraphID);
#endif
    deps.push_back(thisChildNode->CgraphID);
    thisChildNode->touched = true;
    if(downstream | (thisChildNode->type != STOCH)) {
#ifdef _DEBUG_GETDEPS
    PRINTF("          Recursing into child node %i\n", thisChildCgraphID);
#endif
      getDependenciesOneNode(deps, thisChildCgraphID, downstream, recursionDepth + 1);
    }
  }
#ifdef _DEBUG_GETDEPS
  PRINTF("      Done iterating through %i children of node %i\n", numChildren, CgraphID);
#endif
}




