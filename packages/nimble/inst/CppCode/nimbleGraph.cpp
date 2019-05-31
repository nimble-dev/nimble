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
  bool downstream = SEXP_2_bool(Sdownstream);
  vector<int> ans = graphPtr->getDependencies(nodes, omit, downstream);
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
  vector<int>::const_iterator omitFinder; 
  for(i = 0; i < n; i++) {
    thisGraphNodeID = Cnodes[i];

    // Need to check Comit
    // the touching of all Comit nodes still blocks them in the recursion
    // but for the input nodes, we need to check if they are in Comit because
    // being touched could also occur from another input node
    omitFinder = std::find(Comit.begin(), Comit.end(), thisGraphNodeID);
    if(omitFinder != Comit.end()) continue; // it was in omits
    
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
      for(int j = 0; j < result_oneChild.size(); ++j) {
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
  if(nodes[0] >= graphPtr->graphNodeVec.size()) {
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
  for(int i = 0; i < result.size(); ++i) {
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
