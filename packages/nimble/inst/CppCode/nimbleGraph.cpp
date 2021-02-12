/*
 * NIMBLE: an R package for programming with BUGS models.
 * Copyright (C) 2014-2017 Perry de Valpine, Christopher Paciorek,
 * Daniel Turek, Clifford Anderson-Bergman, Nick Michaud, Fritz Obermeyer,
 * Duncan Temple Lang.
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, a copy is available at
 * https://www.R-project.org/Licenses/
 */

#include "nimble/nimbleGraph.h"

graphNode::graphNode(int inputCgraphID, NODETYPE inputType, const string &inputName ) :
  role(UNKNOWNROLE),
  type(inputType),
  CgraphID(inputCgraphID),
  name(inputName),
  touched(false),
  numChildren(0),
  numPaths(-1) {
  RgraphID = CgraphID + 1;
}

void graphNode::addChild( graphNode *toNode, int childParentExpressionID) {
#ifdef _DEBUGNIMGRAPH
  PRINTF("adding child %s to parent %s\n", toNode->name.c_str(), name.c_str());
#endif
  children.push_back(toNode);
  childrenParentExpressionIDs.push_back(childParentExpressionID);
  numChildren++;
  toNode->addParent(this);
}

void graphNode::addParent(graphNode *fromNode) {
#ifdef _DEBUGNIMGRAPH
  PRINTF("adding parent %s to child %s\n", fromNode->name.c_str(), name.c_str());
#endif
  parents.push_back(fromNode);
}

void SEXP_2_nodeType(SEXP Stypes, vector<NODETYPE> &ans) {
  //  enum NODETYPE {UNKNOWNTYPE, STOCH, DETERM, RHSONLY};
  if(!Rf_isString(Stypes)) {
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
    else if(oneString == "unknownIndex")
      ans[i] = UNKNOWNINDEX;
    else if(oneString == "unknown")
      ans[i] = UNKNOWNTYPE;
    else {
      ans[i] = UNKNOWNTYPE;
      PRINTF("In SEXP_2_nodeType: unknown string type label %s\n", oneString.c_str());
    }
  }
}

SEXP C_setGraph(SEXP SedgesFrom, SEXP SedgesTo, SEXP SedgesFrom2ParentExprIDs, SEXP SnodeFunctionIDs, SEXP Stypes, SEXP Snames, SEXP SnumNodes) {
  vector<int> edgesFrom = SEXP_2_vectorInt(SedgesFrom, -1); // -1 subtracted here
  vector<int> edgesTo = SEXP_2_vectorInt(SedgesTo, -1); // -1 substracted here
  vector<int> edgesFrom2ParentExprIDs = SEXP_2_vectorInt(SedgesFrom2ParentExprIDs);
  vector<int> nodeFunctionIDs = SEXP_2_vectorInt(SnodeFunctionIDs, -1);
  vector<NODETYPE> types;
  SEXP_2_nodeType(Stypes, types);
  vector<string> names;
  STRSEXP_2_vectorString(Snames, names);
  int numNodes = SEXP_2_int(SnumNodes);
  nimbleGraph *newGraph = new nimbleGraph;
  newGraph->setNodes(edgesFrom, edgesTo, edgesFrom2ParentExprIDs, nodeFunctionIDs, types, names, numNodes);
  SEXP SextPtrAns;
  PROTECT(SextPtrAns = R_MakeExternalPtr(newGraph, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(SextPtrAns, &nimbleGraphFinalizer, TRUE);
  UNPROTECT(1);
  return(SextPtrAns);
}

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

SEXP C_anyStochDependencies(SEXP SgraphExtPtr) {
  nimbleGraph *graphPtr = static_cast<nimbleGraph *>(R_ExternalPtrAddr(SgraphExtPtr));
  vector<int> ans(graphPtr->anyStochDependencies());
  SEXP Sans;
  PROTECT(Sans = Rf_allocVector(LGLSXP, ans.size()));
  int *SansPtr = INTEGER(Sans);
  for(unsigned int i = 0; i < ans.size(); i++) {
    if(ans[i] == 0) PRINTF("Element %i was not processed\n", i);
    SansPtr[i] = ans[i]==2 ? 1 : 0;
  }
  UNPROTECT(1);
  return(Sans);
}

SEXP C_getDependencies(SEXP SgraphExtPtr, SEXP Snodes, SEXP Somit, SEXP Sdownstream) {
  nimbleGraph *graphPtr = static_cast<nimbleGraph *>(R_ExternalPtrAddr(SgraphExtPtr));
  vector<int> nodes = SEXP_2_vectorInt(Snodes, -1); // subtract 1 index for C
  vector<int> omit = SEXP_2_vectorInt(Somit, -1);
  std::sort(omit.begin(), omit.end());
  bool downstream = SEXP_2_bool(Sdownstream);
  vector<int> ans = graphPtr->getDependencies(nodes, omit, downstream);
  return(vectorInt_2_SEXP(ans, 1)); // add 1 index for R
}

SEXP C_getParents(SEXP SgraphExtPtr, SEXP Snodes, SEXP Somit, SEXP Sdownstream) {
  nimbleGraph *graphPtr = static_cast<nimbleGraph *>(R_ExternalPtrAddr(SgraphExtPtr));
  vector<int> nodes = SEXP_2_vectorInt(Snodes, -1); // subtract 1 index for C
  vector<int> omit = SEXP_2_vectorInt(Somit, -1);
  std::sort(omit.begin(), omit.end());
  bool downstream = SEXP_2_bool(Sdownstream);
  vector<int> ans = graphPtr->getParents(nodes, omit, downstream);
  return(vectorInt_2_SEXP(ans, 1)); // add 1 index for R
}

SEXP C_getDependencyPathCountOneNode(SEXP SgraphExtPtr, SEXP Snode) {
  nimbleGraph *graphPtr = static_cast<nimbleGraph *>(R_ExternalPtrAddr(SgraphExtPtr));
  int node = SEXP_2_int(Snode, 0, -1); // subtract 1 index for C
  int result = graphPtr->getDependencyPathCountOneNode(node);
  return(int_2_SEXP(result)); 
}

SEXP C_anyStochParents(SEXP SgraphExtPtr) {
  nimbleGraph *graphPtr = static_cast<nimbleGraph *>(R_ExternalPtrAddr(SgraphExtPtr));
  vector<int> ans(graphPtr->anyStochParents());
  SEXP Sans;
  PROTECT(Sans = Rf_allocVector(LGLSXP, ans.size()));
  int *SansPtr = INTEGER(Sans);
  for(unsigned int i = 0; i < ans.size(); i++) {
    if(ans[i] == 0) PRINTF("Element %i was not processed\n", i);
    SansPtr[i] = ans[i]==2 ? 1 : 0;
  }
  UNPROTECT(1);
  return(Sans);
}

void nimbleGraph::setNodes(const vector<int> &edgesFrom, const vector<int> &edgesTo,
			   const vector<int> &edgesFrom2ParentExprIDs,
			   const vector<int> &nodeFunctionIDs,
			   const vector<NODETYPE> &types,
			   const vector<string> &names,
			   int inputNumNodes) {
  if(inputNumNodes < 0) PRINTF("Error in setNodes: inputNumNodes < 0\n");
  numNodes = static_cast<unsigned int>(inputNumNodes);
  unsigned int numEdges = edgesFrom.size();

#ifdef _DEBUGNIMGRAPH
  PRINTF("numNodes %i\n", numNodes);
  PRINTF("numEdges %i\n", numEdges);
#endif
  if((numEdges != edgesTo.size()) | (numEdges != edgesFrom2ParentExprIDs.size()) | (numNodes != types.size()) | (numNodes != names.size())) {
    PRINTF("Something is not the right size\n");
    return;
  }
  if(numNodes != nodeFunctionIDs.size()) {
    PRINTF("Wrong length for nodeFunctionIDs\n");
    return;
  }
  graphNodeVec.resize(numNodes);
  for(unsigned int iNode = 0; iNode < numNodes; iNode++) {
    graphNodeVec[iNode] = new graphNode(iNode, types[iNode], names[iNode]);
  }
  for(unsigned int iEdge = 0; iEdge < numEdges; iEdge++) {
    graphNodeVec[ edgesFrom[iEdge]]->addChild( graphNodeVec[edgesTo[iEdge]], edgesFrom2ParentExprIDs[iEdge] );
  }
  for(unsigned int iNode = 0; iNode < numNodes; iNode++) {
    graphNodeVec[iNode]->nodeFunctionNode = graphNodeVec[ nodeFunctionIDs[iNode] ];
  }

}

vector<int> nimbleGraph::anyStochDependencies() {
  vector<int> ans(numNodes, 0);
  for(unsigned int i = 0; i < numNodes; i++) {
    anyStochDependenciesOneNode(ans, i);
  }
  return(ans);
}

bool nimbleGraph::anyStochDependenciesOneNode(vector<int> &anyStochDependencies,  int CgraphID) {
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

bool nimbleGraph::anyStochParentsOneNode(vector<int> &anyStochParents,  int CgraphID) {
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

//#define _DEBUG_GETPATHS

int nimbleGraph::getDependencyPathCountOneNode(const int Cnode) {
  int result(0);
  int i(0);
  graphNode *thisGraphNode;
  graphNode *thisChildNode;

  thisGraphNode = graphNodeVec[ Cnode ];
  if(thisGraphNode->numPaths >= 0)  // already calculated
    return(thisGraphNode->numPaths);

  int numChildren = thisGraphNode->numChildren;
#ifdef _DEBUG_GETPATHS
  PRINTF("debugging output for getDependencyPathCountOneNode with node %i, which has %i children\n", Cnode, numChildren);
#endif

  if(numChildren == 0) {
    result = 0;
  } else {
    for(; i < numChildren; i++) {
      thisChildNode = thisGraphNode->children[i];
#ifdef _DEBUG_GETPATHS
      PRINTF("node %i has child %i\n", Cnode, thisChildNode->CgraphID);
#endif
      if(thisChildNode->type == STOCH) {
        result++;
      } else {
        result += getDependencyPathCountOneNode(thisChildNode->CgraphID); 
      }
    }
  }
  thisGraphNode->numPaths = result; 
  return(result);
}


//#define _DEBUG_GETDEPS

vector<int> nimbleGraph::getDependencies(const vector<int> &Cnodes, const vector<int> &Comit, bool downstream) {
  // assume on entry that touched = false on all nodes
  // Cnodes and Comit are C-indices (meaning they start at 0)
  int n = Comit.size();
  int i;
  vector<int> ans;
  vector<int> tempAns; // This will store LHSinferred nodes, which need to be tracked during recursion but not returned
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
  //  vector<int>::const_iterator omitFinder; 
  for(i = 0; i < n; i++) {
    thisGraphNodeID = Cnodes[i];

    // Need to check Comit
    // the touching of all Comit nodes still blocks them in the recursion
    // but for the input nodes, we need to check if they are in Comit because
    // being touched could also occur from another input node
    if(std::binary_search(Comit.begin(), Comit.end(), thisGraphNodeID)) continue;
    //    omitFinder = std::find(Comit.begin(), Comit.end(), thisGraphNodeID);
    //   if(omitFinder != Comit.end()) continue; // it was in omits
    
    thisGraphNode = graphNodeVec[ thisGraphNodeID ];
#ifdef _DEBUG_GETDEPS
    PRINTF("Working on input node %i\n", thisGraphNodeID);
#endif
    if(!thisGraphNode->touched) {
#ifdef _DEBUG_GETDEPS
      PRINTF("  Adding node %i to ans and recursing\n", thisGraphNodeID);
#endif

      /* LHSINFERRED means e.g. x[1:10] ~ dmnorm() and x[2:3] is used on a RHS. So x[2:3] is LHSINFERRED. IT's not a real node for calculation (no nodeFunction), but it is a vertex in the graph */
      if(thisGraphNode->type != LHSINFERRED) {
	ans.push_back(thisGraphNodeID);
	thisGraphNode->touched = true;
      } else { /* need to include nodeFunctionNode and its non-LHSINFERRED children*/
	/* the current LHSINFERRED node will not be touched or included */
	graphNode* nodeFunctionNode = thisGraphNode->nodeFunctionNode;
	if(!nodeFunctionNode->touched) {
	  int nodeFunctionNodeID = nodeFunctionNode->CgraphID;
	  ans.push_back(nodeFunctionNodeID);
	  nodeFunctionNode->touched = true;
	  getDependenciesOneNode(ans, tempAns, nodeFunctionNodeID, downstream, 1, false);
	}
      }
      getDependenciesOneNode(ans, tempAns, thisGraphNodeID, downstream, 1);
    } else {
#ifdef _DEBUG_GETDEPS
      PRINTF("  Node %i was already touched.\n", thisGraphNodeID);
#endif
      if((thisGraphNode->type == STOCH) & !downstream) {
	/* In this case the input node was already touched, so it is a dependency */
	/* of something earlier on the input list.  But since it was on the input list */
	/* we still need to get its dependencies.  But if downstream is TRUE (==1), then */
	/* its dependencies will have already been pursued so we don't need to */
#ifdef _DEBUG_GETDEPS
      PRINTF("  But is stochastic and downstream is false, so we are recursing into its dependencies.\n");
#endif
      getDependenciesOneNode(ans, tempAns, thisGraphNodeID, downstream, 1);
      }
    }
  }

  // untouch nodes and omit
  n = Comit.size();
  for(i = 0; i < n; i++) {
    graphNodeVec[ Comit[i] ]->touched = false;
  }
  /* The purpose of storing tempAns (LHSINFERRED IDs) was so the touched flags could be cleared here: */
  n = tempAns.size();
  for(i = 0; i < n; i++) {
    graphNodeVec[ tempAns[i] ]->touched = false;
  }
  n = ans.size();
  for(i = 0; i < n; i++) {
    graphNodeVec[ ans[i] ]->touched = false;
  }
  std::sort(ans.begin(), ans.end());
  return(ans);
}


void nimbleGraph::getDependenciesOneNode(vector<int> &deps,
					 vector<int> &tempDeps, /* LHSinferred */
					 int CgraphID,
					 bool downstream,
					 unsigned int recursionDepth,
					 bool followLHSinferred) {
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
    if(!followLHSinferred) {
      if(thisChildNode->type == LHSINFERRED) continue;
    }
    thisChildCgraphID = thisChildNode->CgraphID;
#ifdef _DEBUG_GETDEPS
    PRINTF("        Adding child node %i\n", thisChildCgraphID);
#endif
    if(thisChildNode->type == LHSINFERRED)
      tempDeps.push_back(thisChildNode->CgraphID);
    else
      deps.push_back(thisChildNode->CgraphID); /* LHSINFERRED nodes used to be included here and stripped in R before final return, but there was a bug stripping "split" nodes that have %.s% notation, so now LHSINFERRED are not returned. */
    thisChildNode->touched = true;
    if(downstream | (thisChildNode->type != STOCH)) {
#ifdef _DEBUG_GETDEPS
      PRINTF("          Recursing into child node %i\n", thisChildCgraphID);
#endif
      getDependenciesOneNode(deps, tempDeps, thisChildCgraphID, downstream, recursionDepth + 1);
    }
  }
#ifdef _DEBUG_GETDEPS
  PRINTF("      Done iterating through %i children of node %i\n", numChildren, CgraphID);
#endif
}

/**************/
/* getParents */
/**************/

// #define _DEBUG_GETPARENTS

vector<int> nimbleGraph::getParents(const vector<int> &Cnodes, const vector<int> &Comit, bool upstream) {
  // This is very much like getDependencies, but by default it does *not* include input nodes in output.
    // assume on entry that touched = false on all nodes
  // Cnodes and Comit are C-indices (meaning they start at 0)
  int n = Comit.size();
  int i;
  vector<int> ans;
  vector<int> tempAns; // This will store LHSinferred nodes, which need to be tracked during recursion but not returned
  // touch omit nodes
  for(i = 0; i < n; i++) {
    graphNodeVec[ Comit[i] ]->touched = true;
  }
  n = Cnodes.size();
  graphNode *thisGraphNode;
  int thisGraphNodeID;
  //  vector<int>::const_iterator omitFinder;
#ifdef _DEBUG_GETPARENTS
  std::cout<<"debugging getParents with n = "<<n<<std::endl;
#endif
  for(i = 0; i < n; i++) {
    thisGraphNodeID = Cnodes[i];
#ifdef _DEBUG_GETPARENTS
    std::cout<<"working on input node C-ID = "<<thisGraphNodeID<<std::endl;
#endif
    if(std::binary_search(Comit.begin(), Comit.end(), thisGraphNodeID)) continue;
    //   omitFinder = std::find(Comit.begin(), Comit.end(), thisGraphNodeID);
    //    if(omitFinder != Comit.end()) continue; // it was in omits
    thisGraphNode = graphNodeVec[ thisGraphNodeID ];
    if(!thisGraphNode->touched) { // It is not a parent of another input node
#ifdef _DEBUG_GETPARENTS
      std::cout<<"not touched"<<std::endl;
#endif
      if(thisGraphNode->type != LHSINFERRED) { // It is not a LHSinferred (split) node
#ifdef _DEBUG_GETPARENTS
	std::cout<<"not LHSinferred"<<std::endl;
#endif
	//      ans.push_back(thisGraphNodeID); // This algorithm does not return its input nodes.
	// tempAns.push_back(thisGraphNodeID);  // We do *not* touch and record this node because
	// thisGraphNode->touched = true;       // if it is a parent of another input node, it should be included and traced.
      } else { // It is a LHSinferred node, so behave as if its full declared node was input
#ifdef _DEBUG_GETPARENTS
	std::cout<<"LHSinferred"<<std::endl;
#endif
	graphNode* nodeFunctionNode = thisGraphNode->nodeFunctionNode;
	if(!nodeFunctionNode->touched) {
	  int nodeFunctionNodeID = nodeFunctionNode->CgraphID;
	  //ans.push_back(nodeFunctionNodeID);
	  tempAns.push_back(nodeFunctionNodeID);
	  nodeFunctionNode->touched = true;
	  getParentsOneNode(ans, tempAns, nodeFunctionNodeID, upstream, 1, false);
	  // This imitates getDependencies, but I am not clear if both this and the next call would both be needed
	}
      }
      getParentsOneNode(ans, tempAns, thisGraphNodeID, upstream, 1);
    } else { // This is a parent of another input node
#ifdef _DEBUG_GETPARENTS
      std::cout<<"touched"<<std::endl;
#endif
      if((thisGraphNode->type == STOCH) & !upstream) {
	/* In this case the input node was already touched, so it is a parent */
	/* of something earlier on the input list.  But since it was on the input list */
	/* we still need to get its parents.  But if upstream is TRUE (==1), then */
	/* its dependencies will have already been pursued so we don't need to */
	getParentsOneNode(ans, tempAns, thisGraphNodeID, upstream, 1);
      }
    }
  }
  // untouch nodes and omit
  n = Comit.size();
  for(i = 0; i < n; i++) {
    graphNodeVec[ Comit[i] ]->touched = false;
  }
  /* The purpose of storing tempAns (LHSINFERRED IDs) was so the touched flags could be cleared here: */
  n = tempAns.size();
  for(i = 0; i < n; i++) {
    graphNodeVec[ tempAns[i] ]->touched = false;
  }
  n = ans.size();
  for(i = 0; i < n; i++) {
    graphNodeVec[ ans[i] ]->touched = false;
  }
  std::sort(ans.begin(), ans.end());
  return(ans);
}

void nimbleGraph::getParentsOneNode(vector<int> &deps,
				    vector<int> &tempDeps, /* LHSinferred */
				    int CgraphID,
				    bool upstream,
				    unsigned int recursionDepth,
				    bool followLHSinferred) {
  if(recursionDepth > graphNodeVec.size()) {
    PRINTF("ERROR: getDependencies has recursed too far.  Something must be wrong.\n");
    return;
  }
  graphNode *thisGraphNode = graphNodeVec[CgraphID];
  int numParents = thisGraphNode->parents.size();
  int i(0);
  graphNode *thisParentNode;
  int thisParentCgraphID;
#ifdef _DEBUG_GETPARENTS
  std::cout<<"working on "<< numParents<<" parents"<<std::endl;
#endif
  for(; i < numParents; i++) {
    thisParentNode = thisGraphNode->parents[i];
    if(thisParentNode->touched) continue;
    if(!followLHSinferred) { // I don't think this is relevant in the parent case (only the child case) but it shouldn't hurt
      if(thisParentNode->type == LHSINFERRED) continue;
    }
    thisParentCgraphID = thisParentNode->CgraphID;
#ifdef _DEBUG_GETPARENTS
  std::cout<<"working on parent node C-ID "<< thisParentCgraphID<<std::endl;
#endif
    if(thisParentNode->type == LHSINFERRED)
      tempDeps.push_back(thisParentNode->CgraphID);
    else
      deps.push_back(thisParentNode->CgraphID); 
    thisParentNode->touched = true;
    if(upstream | (thisParentNode->type != STOCH)) {
      getParentsOneNode(deps, tempDeps, thisParentCgraphID, upstream, recursionDepth + 1);
    }
  }
}

/**********************************/
/* getConditionallyIndependentSets*/
/**********************************/

SEXP C_getConditionallyIndependentSets(SEXP SgraphExtPtr,
				       SEXP Snodes,
				       SEXP SgivenNodes,
				       SEXP Somit,
				       SEXP SstartUp,
				       SEXP SstartDown) {
  nimbleGraph *graphPtr = static_cast<nimbleGraph *>(R_ExternalPtrAddr(SgraphExtPtr));
  vector<int> nodes = SEXP_2_vectorInt(Snodes, -1); // subtract 1 index for C
  vector<int> givenNodes = SEXP_2_vectorInt(SgivenNodes, -1); 
  vector<int> omit = SEXP_2_vectorInt(Somit, -1);
  std::sort(omit.begin(), omit.end());
  bool startUp = SEXP_2_bool(SstartUp);
  bool startDown = SEXP_2_bool(SstartDown);
  vector<vector<int> > result = graphPtr->getAllCondIndSets(nodes, givenNodes, omit, startUp, startDown);
  /* sort sets by find node in each set. */
  /* I'm not sure when this would matter, but at least for testing purposes */ 
  /* it is helpful to establish a canonical order of results. */
  struct comp {
    /* compare first elements of results vectors */
    /* empty vectors if they occur go last.*/
    const vector<vector<int> > &result;
    comp(const vector<vector<int> > &result_) : result(result_) {};
    bool operator() (int i,int j) {
      if(result[j].size() == 0) return(i); // second or both empty 
      if(result[i].size() == 0) return(j); // first empty, second non-empty
      return(result[i][0] < result[j][0]);
    }
  };
  vector<int> sort_order(result.size());
  int numEmpty(0);
  for(int i = 0; i < result.size(); ++i) {
    sort_order[i] = i; // sort_order will be 0:number of results-1
    if(result[i].size()==0) ++numEmpty;
  }
  // this will put sort_order as sorting indices of first elements of result vectors
  std::sort(sort_order.begin(), sort_order.end(), comp(result));
  
  SEXP Sresult = PROTECT(Rf_allocVector(VECSXP, result.size() - numEmpty ) );
  for(int i = 0; i < result.size(); ++i) {
    if(result[sort_order[i] ].size() > 0) {
      SET_VECTOR_ELT(Sresult, i, PROTECT(vectorInt_2_SEXP(result[sort_order[i] ], 1)));
    }
  }
  UNPROTECT(1 + result.size() - numEmpty);
  return(Sresult); // add 1 index for R
}

vector<vector<int> > nimbleGraph::getAllCondIndSets(const vector<int> &Cnodes,
						    const vector<int> &CgivenNodes,
						    const vector<int> &Comit,
						    bool startUp,
						    bool startDown) {
  vector<vector<int> > results;
  if(!Cnodes.size()) return results;
  
  // Make isGivenVec.
  vector<bool> isGivenVec(numNodes, false);
  for(int i = 0; i < CgivenNodes.size(); ++i) {
    isGivenVec[ CgivenNodes[i] ] = true;
  }
  
  // Touch Comit for use by all iterations
  for(int i = 0; i < Comit.size(); i++) {
    graphNodeVec[ Comit[i] ]->touched = true;
  }

  // Get first seed node
  int iCurrentInputNode = 0;
  vector<int> inputNodes(1);
  do {
    int inputNodeID = Cnodes[iCurrentInputNode];
    inputNodes[0] = inputNodeID;
    // get conditionally independent set for that node
    results.push_back(getCondIndSet(inputNodes, isGivenVec, Comit, startUp, startDown));
    // Find next available seed node not already in a set
    bool done(false);
    while(!done) {
      ++iCurrentInputNode;
      done = iCurrentInputNode >= Cnodes.size();
      if(!done) done = !(graphNodeVec[ Cnodes[iCurrentInputNode] ]->touched);
    };
  } while ( iCurrentInputNode < Cnodes.size() );

  // untouch the entire graph
  // the book-keeping of touched nodes across multiple cond. ind. sets
  // would be potentially more burdensome than simply untouching everything.
  for(int i = 0; i < numNodes; i++) {
    graphNodeVec[ i ]->touched = false;
  }

  // return all sets
  return results;
}

vector<int> nimbleGraph::getCondIndSet(const vector<int> &Cnodes,
				       const vector<bool>  &isGivenVec,
				       const vector<int> &Comit,
				       bool startUp,
				       bool startDown) {
  // Cnodes are C (0-based) indices for stochastic nodes to seed the search
  // for a conditionally independent set.  It makes most sense for this to be
  // a single node.  If it is multiple nodes, and if they are really in different
  // sets, the result will be a union of their sets.
  //
  // Comit are C indices for nodes to skip in searching.
  //
  // isGivenVec is the length of the graph.  An element is true if we shouldn't traverse through it.
  // That means if looking up, stop recursing up but do look down, and
  // if looking down, stop recursing down but do look up.
  // As initially envisioned, isGivenVec would be true for top nodes and data nodes,
  // resulting in conditionally independent sets of latent nodes.
  //
  // Cnodes should include only stochastic nodes.  This should be checked
  // prior to entry to this function.  If it doesn't, reasonable results
  // should still be returned.
  //
  // startUp and startDown say whether recursion should start up and/or start down
  // If Cnodes are latent nodes, then startUp and startDown should both be true.
  // If e.g. Cnodes give data nodes, then startUp should be true and startDown false.
  //
  // In this algo, LHSINFERRED nodes can be treated like any others.
  int i;
  vector<int> ans;
  // omit nodes were already touched
  int n = Cnodes.size();
  graphNode *thisGraphNode;
  int thisGraphNodeID;
  //  vector<int>::const_iterator omitFinder;
  for(i = 0; i < n; i++) {
    thisGraphNodeID = Cnodes[i];
    if(std::binary_search(Comit.begin(), Comit.end(), thisGraphNodeID)) continue;
    //    omitFinder = std::find(Comit.begin(), Comit.end(), thisGraphNodeID);
    //    if(omitFinder != Comit.end()) continue; // it was in omits
    thisGraphNode = graphNodeVec[ thisGraphNodeID ];
    if(!thisGraphNode->touched) { // It has not been found starting from another input node
      bool isGiven = isGivenVec[thisGraphNodeID];
      if(thisGraphNode->type == STOCH && (!isGiven))
	ans.push_back(thisGraphNodeID);
      thisGraphNode->touched = true;
      expandCondIndSet(ans, thisGraphNodeID, startUp, startDown, isGivenVec, 1);
    }
  }
  /* The purpose of storing tempAns  was so the touched flags could be cleared here: */
  std::sort(ans.begin(), ans.end());
  return(ans);
}

void nimbleGraph::expandCondIndSet(vector<int> &deps,
				   int CgraphID,
				   bool goUp,
				   bool goDown,
				   const vector<bool> &isGivenVec,
				   unsigned int recursionDepth) {
  // We don't need a tempDeps in this algo b/c we untouch the entire graph when done.
  graphNode *thisGraphNode = graphNodeVec[CgraphID];
  for(int dir = 0; dir < 2; ++dir) { // 0 for down, 1 for up
    bool goingDown = dir == 0;
    if(goingDown && (!goDown)) continue;
    if((!goingDown) && (!goUp)) continue;
    int numRelatives;
    graphNode *thisRelNode;
    int thisRelCgraphID;
    if(goingDown)
      numRelatives = thisGraphNode->numChildren;
    else
      numRelatives = thisGraphNode->parents.size();
    int i(0);
    for(; i < numRelatives; i++) {
      if(goingDown)
	thisRelNode = thisGraphNode->children[i];
      else
	thisRelNode = thisGraphNode->parents[i];
      
      if(thisRelNode->touched) continue; // If it's already been handled, continue

      thisRelCgraphID = thisRelNode->CgraphID;
      bool isGiven = isGivenVec[thisRelCgraphID]; 
      if(thisRelNode->type == STOCH && (!isGiven)) // Record latent stochastic nodes for results
	deps.push_back(thisRelNode->CgraphID);
      
      thisRelNode->touched = true;

      // If we're looking down, we always want to recurse up.
      // If we're looking up, we always fully stop at a given node.
      // If we're looking down, we recurse down if the isn't given.
      // The whole recursion ends when every parent and child was already touched (processed).
      // Note that other parent nodes of a given node found going down (a data node) need inclusion,
      // but other child nodes of a given node found going up (a top node) do not need inclusion.
      bool recurseUp = goingDown || ( (!goingDown) && (!isGiven) );
      bool recurseDown = !isGiven;
      
      if(recurseUp || recurseDown) {
	expandCondIndSet(deps, thisRelCgraphID,
			 recurseUp, recurseDown,
			 isGivenVec, recursionDepth + 1);
      }
    }
  }
}

/**********************/
/* getDependencyPaths */
/**********************/

class depStep_class {
private:
  int nodeID_and_parentExprID[2]; // nodeID is first element.  parentExprID is second element
public:
  int& nodeID() {return nodeID_and_parentExprID[0];}
  const int& nodeID() const {return nodeID_and_parentExprID[0];}
  int& parentExprID() {return nodeID_and_parentExprID[1];}
  const int& parentExprID() const {return nodeID_and_parentExprID[1];}
  depStep_class() {};
  depStep_class(int nodeID_, int parentExprID_) {
    nodeID() = nodeID_;
    parentExprID() = parentExprID_;
  }
};

typedef vector<depStep_class> depPath;
typedef vector<depPath> multipleDepPaths;

multipleDepPaths getDependencyPaths_recurse(const graphNode *currentNode, depPath &root_path, int parentExprID) {
  multipleDepPaths result;
  depStep_class thisStep(currentNode->RgraphID, parentExprID);
  root_path.push_back(thisStep);
  if((currentNode->type == STOCH && root_path.size() > 1)) { // stop at non-first stochastic node
    result.push_back(root_path);
  } else {
    // Recurse through child nodes, if there are any.
    // If there are zero children, we have a terminal deterministic node, which should not be recorded.
    for(unsigned int i = 0; i < currentNode->numChildren; ++i) {
      multipleDepPaths result_oneChild = getDependencyPaths_recurse(currentNode->children[i],
								    root_path,
								    currentNode->childrenParentExpressionIDs[i]);      
      for(unsigned int j = 0; j < result_oneChild.size(); ++j) {
	result.push_back( result_oneChild[j] );
      }
    }
  }
  root_path.pop_back();
  return result;
}

SEXP C_getDependencyPaths(SEXP SgraphExtPtr, SEXP Snodes) {
  nimbleGraph *graphPtr = static_cast<nimbleGraph *>(R_ExternalPtrAddr(SgraphExtPtr));
  vector<int> nodes = SEXP_2_vectorInt(Snodes, -1); // subtract 1 index for C
  if(nodes.size() != 1) {
    PRINTF("Input to C_getDependencyPaths should be one and only one nodeID.");
    return R_NilValue;
  }
  if(nodes[0] >= static_cast<int>(graphPtr->graphNodeVec.size())) {
    PRINTF("Input to C_getDependencyPaths has a nodeID that is too large.");
    return R_NilValue;

  }
  if(graphPtr->graphNodeVec[ nodes[0] ]->numChildren == 0) {
    return R_NilValue;
  }
  depPath root_path;
  multipleDepPaths result = getDependencyPaths_recurse(graphPtr->graphNodeVec[ nodes[0] ],
  						       root_path,
  						       INT_MIN); //R_defines.h says INT_MIN is NA_INTEGER
  SEXP Sresult = PROTECT(Rf_allocVector(VECSXP, result.size()));
  SEXP Sdim;
  for(unsigned int i = 0; i < result.size(); ++i) {
    SEXP SdepPath = PROTECT(Rf_allocVector(INTSXP, 2*result[i].size()));
    int *depPathPtr = INTEGER(SdepPath);
    int thisResultSize = result[i].size();
    for(int j = 0; j < thisResultSize; ++j) {
      depPathPtr[j] = result[i][j].nodeID();
      depPathPtr[j + thisResultSize] = result[i][j].parentExprID();
    }
    Sdim = PROTECT(Rf_allocVector(INTSXP, 2));
    INTEGER(Sdim)[0] = thisResultSize;
    INTEGER(Sdim)[1] = 2;
    Rf_setAttrib(SdepPath, R_DimSymbol, Sdim);
    SET_VECTOR_ELT(Sresult, i, SdepPath);
    UNPROTECT(2);
  }
  UNPROTECT(1);
  return Sresult;
}
