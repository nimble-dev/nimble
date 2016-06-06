#ifndef _NIMBLEGRAPH
#define _NIMBLEGRAPH
#include "RcppUtils.h"
#include "R.h"
#include<Rinternals.h>
#include<vector>
#include<algorithm>
#include<string>
using std::vector;
using std::string;

/* NODETYPE labels the type of node, regardless of where it is in the graph */
enum NODETYPE {UNKNOWNTYPE, STOCH, DETERM, RHSONLY, LHSINFERRED};
/* NODEROLE labels how a node fits in the graph */
enum NODEROLE {UNKNOWNROLE, TOP, LATENT, END, DATA};

struct graphNode {
 public:
  NODEROLE role;
  NODETYPE type;
  int RgraphID;
  int CgraphID; // always RgraphID-1
  string name;
  bool touched; /* This is for arbitrary use by graph-traversing algorithms.  By convention it should be left at false for all nodes after completion of an algorithm. */
  int numChildren;
  vector<graphNode*> children; /* pointers to child nodes */
  vector<int> childrenParentExpressionIDs; /* integer labels of how this node is used by child nodes. */
  vector<graphNode*> parents; /* pointers to parent nodes*/
  graphNode(int inputCgraphID, NODETYPE inputType, const string &inputName);
  void addChild(graphNode *toNode, int childParentExpressionID);
  void addParent(graphNode *fromNode);
};

struct nimbleGraph {
public:
  vector<graphNode*> graphNodeVec;
  int numNodes;
  void setNodes(const vector<int> &edgesFrom, const vector<int> &edgesTo,
		const vector<int> &edgesFrom2ParentExprIDs,
		const vector<NODETYPE> &types,
		const vector<string> &names,
		int inputNumNodes);
  vector<int> anyStochDependencies();
  bool anyStochDependenciesOneNode(vector<int> &anyStochDependencies, int CgraphID);
  vector<int> anyStochParents();
  bool anyStochParentsOneNode(vector<int> &anyStochParents, int CgraphID);
  vector<int> getDependencies(const vector<int> &Cnodes, const vector<int> &Comit, bool downstream);
  void getDependenciesOneNode(vector<int> &deps, int CgraphID, bool downstream, int recursionDepth);
  ~nimbleGraph();
};

void nimbleGraphFinalizer(SEXP SgraphExtPtr);

extern "C" {
  SEXP setGraph(SEXP SedgesFrom, SEXP SedgesTo, SEXP SedgesFrom2ParentExprIDs, SEXP Stypes, SEXP Snames, SEXP SnumNodes); 
  SEXP anyStochDependencies(SEXP SextPtr);
  SEXP anyStochParents(SEXP SextPtr);
  SEXP getDependencies(SEXP SextPtr, SEXP Snodes, SEXP Somit, SEXP Sdownstream);
}

#endif
