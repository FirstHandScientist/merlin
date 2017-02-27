/*
 * any_anaobf.h
 *
 *  Created on: 9 Oct 2015
 *      Author: radu
 */

#ifndef SRC_ANY_AAOBF_H_
#define SRC_ANY_AAOBF_H_


#include "aobf2.h"

#define MARK_FEASIBLE_TIGHTEST 1
#define MARK_FEASIBLE_EPSILON 2

/**
 * This is an implementation of the anytime non-parametric AOBF algorithm
 * which is a hybrid of depth-first and best-first search on the AO graph.
 */

class AnyAAOBF : public AOBF2 {

protected:

	// Tip nodes of the current feasible partial solution tree
	std::vector<AobfSearchNode*> m_tipsF;

	// Current assignment of FST
	std::vector<int> m_assignmentF;

	// Number of full repairs
	size_t m_numRepairs;

	// Number of expansions from FST (depth-first)
	size_t m_numExpansionsFST;

	// Number of expansions from BST (best-first)
	size_t m_numExpansionsBST;

protected:

	// AO* search
	int aostar(bool verbose = true);

	// expand a node
	bool expand(AobfSearchNode *node);

	// revise the q-value of a node
	bool revise(AobfSearchNode *node);

	// expand and revise a node (from the BST)
	void expandAndRevise(AobfSearchNode *node);

	// expand and revise a node (from the FST)
	void expandAndReviseF(AobfSearchNode* node);

	// arrange tip nodes
	void arrangeTipNodesF();

	// select a tip node
	AobfSearchNode* chooseTipNodeF();

	// propagate the q-values
	void update(AobfSearchNode* node);

	// propagate the f-values following the expansion of a terminal node (no markings)
	void updateF(AobfSearchNode* node);

	// revise the f-value of a node (without 'feasible' marking)
	bool reviseF(AobfSearchNode *node);

	// recompute q-values starting from leaf nodes (correct markings as well)
	bool repair();

	// redo feasible markings to get the most promissing partial solution tree
	bool repairF();

	// mark the best AND child of an OR node (ie, the max e-value)
	bool markFeasibleChild(AobfSearchNode *node);

	// compute the best partial solution tree by tracing down the max e-values
	bool findFeasiblePartialTree();

public:

	// solve
	int solve();

public:
	AnyAAOBF(ProgramOptions *opt);
	~AnyAAOBF();

};

inline AnyAAOBF::AnyAAOBF(ProgramOptions *opt) : AOBF2(opt) {
	m_numRepairs = 0;
	m_numExpansionsFST = 0;
	m_numExpansionsBST = 0;
}

inline AnyAAOBF::~AnyAAOBF() {

}


#endif /* SRC_ANY_ANAOBF_H_ */
