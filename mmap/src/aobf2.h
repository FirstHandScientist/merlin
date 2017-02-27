/*
 * aobf.h
 *
 *  Created on: Apr 17, 2013
 *      Author: radu
 */

#ifndef IBM_MAP_AOBF2_H_
#define IBM_MAP_AOBF2_H_

#include "search.h"
#include "aobf_search_node.h"
#include "aobf_search_space.h"
#include "aobb_search_node.h"
#include "aobb_search_space.h"
#include "bound_propagator.h"

/* Standard Best-First AND/OR Search solver (AOBF) - OR caching */

class AOBF2 : public Search {

protected:

	// The context-minimal AND/OR graph
	scoped_ptr<AobfSearchSpace2> m_space;

	// Tip nodes of the current best partial solution tree
	std::vector<AobfSearchNode*> m_tips;

	// global search index
	size_t m_globalSearchIndex;

	// weight that inflates the heuristic (default is 1.0)
	double m_epsilon;

	// number of deadends due to constraint propagation
	size_t m_deadendsCP;

	// number of SUM evaluations
	size_t m_numSumEvals;

	// number of cache hits (MAP)
	size_t m_numCacheHits;

public:

	// initialize
	void init();

	// solver
	virtual int solve();

protected:

	// AO* search over the AND/OR graph
	virtual int aostar(bool verbose = false);

	// expand a node
	virtual bool expand(AobfSearchNode* node);

	// revise the value of a node
	virtual bool revise(AobfSearchNode* node);

	// calculate the heuristic value of an OR node
	double heuristicOR(AobfSearchNode* n);

	// expand and revise
	virtual void expandAndRevise(AobfSearchNode* node);

	// find best partial solution tree
	virtual bool findBestPartialTree();

	// arrange tip nodes
	void arrangeTipNodes();

	// select a tip node
	AobfSearchNode* chooseTipNode();

	// get the assigment along the current path from the root
	void getPathAssignment(AobfSearchNode* node, std::vector<int>& vars);

	// check if variable assignment is a deadend (constraint propagation)
	virtual bool isDeadend(AobfSearchNode* node);
	bool isDeadend(AobfSearchNode* node, int var, int val);

	// returns a string of the context instantiation (AND context)
	std::string context(int node_type, const std::set<int>& C);

	// comparison function for sorting tip nodes (ascending)
	struct CompNodeHeurAsc : public std::binary_function<AobfSearchNode*, AobfSearchNode*, bool> {
		bool operator()(const AobfSearchNode *x, const AobfSearchNode *y) const {
			return (x->getHeur() < y->getHeur());
		}
	};

	// comparison function for sorting tip nodes (descending)
	struct CompNodeHeurDesc : public std::binary_function<AobfSearchNode*, AobfSearchNode*, bool> {
		bool operator()(const AobfSearchNode *x, const AobfSearchNode *y) const {
			return (x->getHeur() > y->getHeur());
		}
	};

	// comparison function (ascending)
	struct CompNodeIndexAsc : public std::binary_function<AobfSearchNode*, AobfSearchNode*, bool> {
		bool operator()(const AobfSearchNode *x, const AobfSearchNode *y) const {
			return (x->getIndex() < y->getIndex());
		}
	};

	// comparison function (descending)
	struct CompNodeIndexDesc : public std::binary_function<AobfSearchNode*, AobfSearchNode*, bool> {
		bool operator()(const AobfSearchNode *x, const AobfSearchNode *y) const {
			return (x->getIndex() > y->getIndex());
		}
	};

	// comparison function (ascending)
	struct CompNodeDepthAsc : public std::binary_function<AobfSearchNode*, AobfSearchNode*, bool> {
		bool operator()(const AobfSearchNode *x, const AobfSearchNode *y) const {
			return (x->getDepth() < y->getDepth());
		}
	};

	// comparison function (descending)
	struct CompNodeDepthDesc : public std::binary_function<AobfSearchNode*, AobfSearchNode*, bool> {
		bool operator()(const AobfSearchNode *x, const AobfSearchNode *y) const {
			return (x->getDepth() > y->getDepth());
		}
	};

	// Write solution to output file
#ifdef UAI_COMPETITION
	void writeSolution();
#endif

protected:

	// AO search (for the conditioned SUM problem)
	std::stack<AobbSearchNode*> ao_stack;
	std::vector<AobbSearchNode*> ao_expand;
	scoped_ptr<AobbSearchSpace> ao_space;
	scoped_ptr<BoundPropagator> ao_propagator;

	// Solve the conditioned SUM subproblem by AO search
	double _ao(int var, const std::vector<int>& assignment);

	void _initAoSearch(int var);
	void _heuristicOR(AobbSearchNode*);
	bool _doProcess(AobbSearchNode*);
	void _addCacheContext(AobbSearchNode*, const set<int>&) const;
	bool _doCaching(AobbSearchNode* node);
	AobbSearchNode* _nextLeaf();
	bool _generateChildrenAND(AobbSearchNode*, vector<AobbSearchNode*>&);
	bool _generateChildrenOR(AobbSearchNode* n, vector<AobbSearchNode*>& chi);
	AobbSearchNode* _nextNode();
	bool _doExpand(AobbSearchNode* n);
	bool _doPruning(AobbSearchNode* n);
	bool _canBePruned(AobbSearchNode* n);
	void _getPathAssignment(AobbSearchNode* node, std::vector<int>& vars);

public:
	AOBF2(ProgramOptions* o);
	virtual ~AOBF2();

protected:
	double m_prevTimeCheckpoint;
	double m_prevNodeCheckpoint;

};

/* inline functions */

inline AOBF2::AOBF2(ProgramOptions* o) : Search(o) {
	m_prevTimeCheckpoint = 0;
	m_prevNodeCheckpoint = 0;
	m_globalSearchIndex = 0;
	m_epsilon = o->weight; // default
	m_deadendsCP = 0;
	m_numSumEvals = 0;
	m_numCacheHits = 0;
}

inline AOBF2::~AOBF2() {

}

#endif /* AOBF_H_ */
