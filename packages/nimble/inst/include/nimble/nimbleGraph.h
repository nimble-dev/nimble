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
enum NODETYPE {UNKNOWNTYPE, STOCH, DETERM, RHSONLY, LHSINFERRED, UNKNOWNINDEX};
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
  unsigned int numChildren;
  graphNode *nodeFunctionNode; /* If this is a split vertex, point to declared node it is part of.   Otherwise point to self.*/
  vector<graphNode*> children; /* pointers to child nodes */
  vector<int> childrenParentExpressionIDs; /* integer labels of how this node is used by child nodes. */
  vector<graphNode*> parents; /* pointers to parent nodes*/
  int numPaths;  /* for use in path counting -- number of paths starting at the node and terminating at stochastic nodes */
  graphNode(int inputCgraphID, NODETYPE inputType, const string &inputName);
  void addChild(graphNode *toNode, int childParentExpressionID);
  void addParent(graphNode *fromNode);
};

struct nimbleGraph {
public:
  vector<graphNode*> graphNodeVec;
  unsigned int numNodes;
  void setNodes(const vector<int> &edgesFrom, const vector<int> &edgesTo,
		const vector<int> &edgesFrom2ParentExprIDs,
		const vector<int> &nodeFunctionIDs,
		const vector<NODETYPE> &types,
		const vector<string> &names,
		int inputNumNodes);
  vector<int> anyStochDependencies();
  bool anyStochDependenciesOneNode(vector<int> &anyStochDependencies, int CgraphID);
  vector<int> anyStochParents();
  bool anyStochParentsOneNode(vector<int> &anyStochParents, int CgraphID);
  vector<int> getDependencies(const vector<int> &Cnodes, const vector<int> &Comit, bool downstream);
  void getDependenciesOneNode(vector<int> &deps, vector<int> &tempDeps, int CgraphID, bool downstream, unsigned int recursionDepth, bool followLHSinferred = true);
  vector<int> getParents(const vector<int> &Cnodes, const vector<int> &Comit, bool upstream, bool oneStep);
  void getParentsOneNode(vector<int> &deps, vector<int> &tempDeps, int CgraphID, bool upstream, unsigned int recursionDepth, bool recurse = true, bool followLHSinferred = true);
  int getDependencyPathCountOneNode(const int Cnode);
  
  vector<vector<int> > getAllCondIndSets(const vector<int> &Cnodes,
                                         const vector<int> &CgivenNodes,
                                         const vector<int> &Comit,
                                         bool startUp,
                                         bool startDown,
                                         bool unknownsAsGiven);
  vector<int> getCondIndSet(const vector<int> &Cnodes,
                            const vector<bool>  &isGivenVec,
                            const vector<bool>  &isLatentVec,
                            const vector<int> &Comit,
                            bool startUp,
                            bool startDown,
                            bool unknownsAsGiven);
  void expandCondIndSet(vector<int> &deps,
                        int CgraphID,
                        bool goUp,
                        bool goDown,
                        const vector<bool> &isGivenVec,
                        const vector<bool> &isLatentVec,
                        bool unknownsAsGiven,
                        unsigned int recursionDepth);
  ~nimbleGraph();
};

void nimbleGraphFinalizer(SEXP SgraphExtPtr);

extern "C" {
  SEXP C_setGraph(SEXP SedgesFrom, SEXP SedgesTo, SEXP SedgesFrom2ParentExprIDs, SEXP SnodeFunctionIDs, SEXP Stypes, SEXP Snames, SEXP SnumNodes);
  SEXP C_anyStochDependencies(SEXP SextPtr);
  SEXP C_anyStochParents(SEXP SextPtr);
  SEXP C_getDependencies(SEXP SextPtr, SEXP Snodes, SEXP Somit, SEXP Sdownstream);
  SEXP C_getParents(SEXP SextPtr, SEXP Snodes, SEXP Somit, SEXP upstream, SEXP SoneStep);
  SEXP C_getDependencyPathCountOneNode(SEXP SgraphExtPtr, SEXP Snode);
  SEXP C_getConditionallyIndependentSets(SEXP SgraphExtPtr,
                                         SEXP Snodes,
                                         SEXP SgivenNodes,
                                         SEXP Somit,
                                         SEXP SstartUp,
                                         SEXP SstartDown,
                                         SEXP SunknownsAsGiven);
}

/**********************/
/* getDependencyPaths */
/**********************/

extern "C" {
  SEXP C_getDependencyPaths(SEXP SgraphExtPtr, SEXP Snodes);
}

#endif
