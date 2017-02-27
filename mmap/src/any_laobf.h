/*
 * any_laobf.h
 *
 *  Created on: 29 Sep 2015
 *      Author: radu
 */

#ifndef SRC_ANY_LAOBF_H_
#define SRC_ANY_LAOBF_H_


#include "search.h"
#include "aobf_search_node.h"
#include "aobf_search_space.h"
#include "aobb_search_node.h"
#include "aobb_search_space.h"
#include "bound_propagator.h"


/**
 * Anytime AOBF with depth-first lookaheads (OR contexts)
 *
 */
class AnyLAOBF : public Search {

protected:

	// The context-minimal AND/OR graph
	scoped_ptr<AobfSearchSpace2> m_spaceBFS;

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

	// g-value of the best partial tree
	double m_g;

	// number of lookaheads
	size_t m_numLookaheads;

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
	//bool isDeadend(AobfSearchNode* node, int var, int val);

	// check if the AND node (var, val) can be pruned below the current OR node
	bool canBePruned(AobfSearchNode* n, int var, int val);

	// Update ancestors values (q-values) along marked paths
	void update(AobfSearchNode* node);

	// Update all ancestors values (q-values)
	void updateAll(AobfSearchNode* node);

	bool check();

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

private:

	// DFS data structures (for BnB probes as well as the conditioned SUMs)
	std::stack<AobbSearchNode*> m_stack;
	std::vector<AobbSearchNode*> m_expand;
	scoped_ptr<AobbSearchSpace> m_spaceDFS;
	scoped_ptr<BoundPropagator> m_propagator;

	// BRAOBB rotating stacks
	size_t m_stackCount;        // counter for current stack
	size_t m_stackLimit;        // expansion limit for stack rotation
	MyStack* m_rootStack;       // the root stack
	queue<MyStack*> m_stacks;   // the queue of active stacks

	// Solve the conditioned SUM subproblem by depth-first AO search
	double _sum(int var, const std::vector<int>& assignment);
	double _probe(int var, const std::vector<int>& assignment);

	AobbSearchNode* _initSearchSpace(int var);
	double _heuristicOR(AobbSearchNode*);
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
	void _reset(AobbSearchNode* p);

public:
	AnyLAOBF(ProgramOptions* o);
	virtual ~AnyLAOBF();

protected:
	double m_prevTimeCheckpoint;
	double m_prevNodeCheckpoint;

};


/* inline functions */

inline AnyLAOBF::AnyLAOBF(ProgramOptions* o) : Search(o) {
	m_prevTimeCheckpoint = 0;
	m_prevNodeCheckpoint = 0;
	m_globalSearchIndex = 0;
	m_epsilon = o->weight; // default
	m_deadendsCP = 0;
	m_numSumEvals = 0;
	m_numCacheHits = 0;
	m_stackCount = 0;
	m_stackLimit = o->rotateLimit; // default: 1000 node expansions per stack
	m_rootStack = NULL;
	m_g = 0;
	m_numLookaheads = 0;
}

inline AnyLAOBF::~AnyLAOBF() {

}

inline void AnyLAOBF::_reset(AobbSearchNode* p) {
	assert(p);
	while (m_stacks.size()) {
		delete m_stacks.front();
		m_stacks.pop();
	}
	m_rootStack = new MyStack(NULL);
	m_rootStack->push(p);
	m_stacks.push(m_rootStack);
}


#endif /* SRC_ANY_LAOBF_H_ */
