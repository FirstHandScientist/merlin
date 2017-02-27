/*
 * aobb.h
 *
 *  Created on: Jun 7, 2013
 *      Author: radu
 */

#ifndef IBM_MAP_AOBB_H_
#define IBM_MAP_AOBB_H_

#include "search.h"
#include "aobb_search_node.h"
#include "aobb_search_space.h"
#include "bound_propagator.h"

// AND/OR Branch and Bound (AOBB) with mini-bucket heuristics
// (assumes a constrained pseudo tree, ie. all MAP variables at the top)

class AOBB : public Search {

protected:

	// The DFS stack of nodes
	std::stack<AobbSearchNode*> m_stack;

	// Reusable vector for node expansions (to avoid repeated (de)allocation of memory)
	std::vector<AobbSearchNode*> m_expand;

	// The context minimal AND/OR search graph (OR caching)
	scoped_ptr<AobbSearchSpace> m_space;

	// cache hits
	size_t m_cacheHits;
	size_t m_deadendsCP; // due to determinism
	size_t m_deadendsUB; // due to upper bounds

	// Bound propagator
	scoped_ptr<BoundPropagator> m_propagator;

	std::vector<val_t> m_prevAssignment;

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

	void addCacheContext(AobbSearchNode* node, const set<int>& ctxt) const;

	// returns the next leaf node
	AobbSearchNode* nextLeaf();

	// get the path assignment (as SAT variables)
	void getPathAssignment(AobbSearchNode* node, std::vector<int>& vars);

public:

	// initialize
	void init();

	// Main solver routine
	virtual int solve();

public:

	AOBB(ProgramOptions* opt);
	virtual ~AOBB();

protected:

	// Time and node checkpoints
	double m_prevTimeCheckpoint;
	double m_prevNodeCheckpoint;
};

/* inline definitions */

inline AOBB::AOBB(ProgramOptions* opt) : Search(opt) {

	// Preallocate space for expansion vector. 128 should be plenty.
	m_expand.reserve(128);
	m_cacheHits = 0;
	m_deadendsCP = 0;
	m_deadendsUB = 0;
	m_prevNodeCheckpoint = 0;
	m_prevTimeCheckpoint = 0;
}

inline AOBB::~AOBB() {

}

inline bool AOBB::isDone() const {
	return ( m_stack.empty() );
}

#endif /* AOBB_H_ */
