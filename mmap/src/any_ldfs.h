/*
 * any_lrtao.h
 *
 *  Created on: 4 Mar 2016
 *      Author: radu
 */

#ifndef SRC_ANY_LDFS_H_
#define SRC_ANY_LDFS_H_

#include "aobf2.h"

/**
 * Learning Depth-First AO search (LDFS)
 * Expand depth-first the current best partial solution tree until a complete
 * solution is found. Then, propagate the values of the leaf nodes to the
 * ancestors, and update the best markings on the OR-to-AND arcs. If the
 * heuristic evaluation function on the current best partial solution tree
 * yields a lower bound that is larger than the current best upper bound,
 * then expand any of its tip nodes, and update the values of the ancestors
 * while redoing the markings. Continue the search until the lower bound
 * equals the upper bound. Notice that LDFS works fine if the pruning is
 * disabled. In that case the lower bound will increase slower than before
 * because it will expand nodes that don't actually contribute to improving
 * neither the best upper bound nor the best lower bound.
 *
 */

class AnyLDFS : public AOBF2 {

protected:

	// Leaf nodes of the current solution tree (ie. SUM subproblems)
	std::vector<AobfSearchNode*> m_leaves;

protected:

	// AO* search
	int aostar(bool verbose = true);

	// expand a node
	bool expand(AobfSearchNode *node);

	// revise the q-value of a node (also marks the best successor)
	bool revise(AobfSearchNode *node);

	// expand and revise the q-value of a node
	void expandAndRevise(AobfSearchNode *node);

	// expand and update the q-values of the node and its ancestors
	void expandAndUpdate(AobfSearchNode* node);

	// propagate the q-values and u-values (should be triggered by a leaf node)
	void update(AobfSearchNode* node);

	// mark the best AND child of an OR node (ie, min q-value)
	bool setBestChild(AobfSearchNode *node);

	// compute the current best partial solution tree by tracing down markings
	bool findBestPartialTree();

	// check if a tip node is deadend
	bool isDeadend(AobfSearchNode* node);

public:

	// solve
	int solve();

public:
	// constructor/destructor
	AnyLDFS(ProgramOptions* opt);
	~AnyLDFS();

};

inline AnyLDFS::AnyLDFS(ProgramOptions* opt) : AOBF2(opt) {

}

inline AnyLDFS::~AnyLDFS() {

}

#endif /* SRC_ANY_LDFS_H_ */
