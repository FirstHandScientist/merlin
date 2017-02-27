/*
 * any_bbmb.h
 *
 *  Created on: 28 Oct 2015
 *      Author: radu
 */

#ifndef SRC_ANY_BBMB_BOUNDED_H_
#define SRC_ANY_BBMB_BOUNDED_H_

#include "search.h"
#include "aobb_search_node.h"
#include "aobb_search_space.h"
#include "bound_propagator.h"


/**
 * Anytime depth-first OR Branch-and-Bound over the MAP variables, guided
 * by weighted mini-bucket heuristics. The conditioned summation is lower bounded
 * namely it is not solved exactly anymore.
 */

class AnyBBMB_Bounded : public Search {
protected:

	// The DFS stack of nodes
	std::stack<AobbSearchNode*> m_stack;

	// Reusable vector for node expansions (to avoid repeated (de)allocation of memory)
	std::vector<AobbSearchNode*> m_expand;

	// The context minimal AND/OR search graph (OR caching)
	scoped_ptr<AobbSearchSpace> m_space;

	// cache hits
	size_t m_cacheHits;

	// Bound propagator
	scoped_ptr<BoundPropagator> m_propagator;

	// Last MAP variable
	int m_lastMapVar;

protected:

	// calculate the h-value of an OR node; cache the h-values of its AND children
	double heuristicOR(AobbSearchNode* node);

	// Initialize the search
	AobbSearchNode* initSearchSpace();

	// Checks if search is done
	virtual bool isDone() const;

	// Returns the next node to expand
	virtual AobbSearchNode* nextNode();

	// Expands a node
	virtual bool doExpand(AobbSearchNode* n);

	// Processes a node
	bool doProcess(AobbSearchNode* n);

	// Prunes a node
	bool doPruning(AobbSearchNode* n);

	// Caches a node
	bool doCaching(AobbSearchNode* n);

	// generates the children of an AND node, writes them into argument vector,
	// returns true if no children
	bool generateChildrenAND(AobbSearchNode*, std::vector<AobbSearchNode*>&);
	// generates the children of an OR node, writes them into argument vector,
	// returns true if no children
	bool generateChildrenOR(AobbSearchNode*, std::vector<AobbSearchNode*>&);

	// checks if the node can be pruned (only meant for AND nodes)
	bool canBePruned(AobbSearchNode*);

	// get the cache context of a node (ie, assignment)
	void addCacheContext(AobbSearchNode* node, const set<int>& ctxt) const;

	// returns the next leaf node
	AobbSearchNode* nextLeaf();

public:

	// initialize
	void init();

	// Main solver routine
	virtual int solve();

public:

	AnyBBMB_Bounded(ProgramOptions* opt);
	virtual ~AnyBBMB_Bounded();

protected:

	// Time and node checkpoints
	double m_prevTimeCheckpoint;
	double m_prevNodeCheckpoint;

};

/* inline definitions */

inline AnyBBMB_Bounded::AnyBBMB_Bounded(ProgramOptions* opt) : Search(opt) {

	// Preallocate space for expansion vector. 128 should be plenty.
	m_expand.reserve(128);
	m_cacheHits = 0;
	m_lastMapVar = UNKNOWN;
	m_prevNodeCheckpoint = 0;
	m_prevTimeCheckpoint = 0;
}

inline AnyBBMB_Bounded::~AnyBBMB_Bounded() {

}

inline bool AnyBBMB_Bounded::isDone() const {
	return ( m_stack.empty() );
}


#endif /* SRC_ANY_BBMB_BOUNDED_H_ */
