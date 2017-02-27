/*
 * any_laobf.h
 *
 *  Created on: 29 Sep 2015
 *      Author: radu
 */

#ifndef SRC_ANY_SAOBF_H_
#define SRC_ANY_SAOBF_H_


#include "search.h"
#include "aobf_search_node.h"
#include "aobf_search_space.h"
#include "aobb_search_node.h"
#include "aobb_search_space.h"
#include "bound_propagator.h"


/**
 * Anytime Stochastic AOBF (OR contexts)
 *
 */
class AnySAOBF : public Search {

	// A priority queue for outside tip nodes
	class priority_queue {
	protected:

		// comparison function for sorting tip nodes in a priority queue
		struct CompNodeHeur : public std::binary_function<AobfSearchNode*, AobfSearchNode*, bool> {
			bool operator()(const AobfSearchNode* x, const AobfSearchNode *y) const {
				return x->getHeur() < y->getHeur();
			}
		};

		typedef std::vector<std::multiset<AobfSearchNode*, CompNodeHeur> > index_t;

		// tip nodes indexed by depth and sorted by heuristic value
		index_t m_nodes;

		// number of elements
		size_t m_size;

	public:
		// constructor and destructor
		priority_queue() : m_size(0) {};
		~priority_queue() {};

		// init the data structure
		inline void init(int sz) {
			m_nodes.resize(sz);
		}

		// add a tip node to the structure
		inline void push(AobfSearchNode* n) {
			int d = n->getDepth();
			assert(d >= 0 && d < (int) m_nodes.size());
			m_nodes[d].insert(n);
			++m_size;
		}

		// erase a tip node from the structure
		inline void erase(AobfSearchNode* n) {
			int d = n->getDepth();
			std::multiset<AobfSearchNode*, CompNodeHeur>::iterator si;
			for (si = m_nodes[d].begin(); si != m_nodes[d].end(); ++si) {
				if (n == *si) {
					m_nodes[d].erase(si);
					if (m_size > 0) --m_size;
					break;
				}
			}
		}

		// pop the deepest tip node having the smallest heuristic value
		inline AobfSearchNode* pop() {

			AobfSearchNode* tip = NULL;
			index_t::reverse_iterator ri = m_nodes.rbegin();
			for (; ri != m_nodes.rend(); ++ri) {
				bool found = false;
				std::multiset<AobfSearchNode*, CompNodeHeur>& tips = *ri;
				if (tips.empty()) {
					continue;
				} else {
					std::multiset<AobfSearchNode*, CompNodeHeur>::iterator si;
					for (si = tips.begin(); si != tips.end(); ++si) {
						tip = *si;
						assert(tip->isFringe() == true);
						tips.erase(si);
						found = true;
						if (m_size > 0) --m_size;
						break;
					}
				}

				if (found) break;
			}

			return tip;
		}

		inline size_t size() {
			return m_size;
		}
	};

protected:
	// comparison function for sorting tip nodes in a priority queue
	struct CompNodePriority : public std::binary_function<AobfSearchNode*, AobfSearchNode*, bool> {
		bool operator()(const AobfSearchNode* x, const AobfSearchNode *y) const {
			if (x->getHeur() > y->getHeur()) {
				return true;
			} else if (x->getHeur() == y->getHeur()) {
				return (x->getDepth() < y->getDepth());
			} else {
				return false;
			}
		}
	};

protected:

	// The context-minimal AND/OR graph
	scoped_ptr<AobfSearchSpace2> m_spaceBFS;

	// Tip nodes of the current best partial solution tree
	std::vector<AobfSearchNode*> m_tips;

	// Tip nodes outside the current best partial solution tree
	priority_queue m_tipsOut;

	// global search index
	size_t m_globalSearchIndex;

	// weight that inflates the heuristic (default is 1.0)
	double m_epsilon;

	// number of SUM evaluations
	size_t m_numSumEvals;

	// number of cache hits (MAP)
	size_t m_numCacheHits;

	// Number of expansions outside BST
	size_t m_numExpansionsOut;

	// Number of expansions inside BST (best-first)
	size_t m_numExpansionsIn;

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

	// revise the q-value of a node
	virtual bool revise(AobfSearchNode* node);

	// revise the f-value of a node
	virtual bool reviseF(AobfSearchNode* node);

	// calculate the heuristic value of an OR node
	double heuristicOR(AobfSearchNode* n);

	// expand and revise
	virtual void expandAndRevise(AobfSearchNode* node);

	// find best partial solution tree
	virtual bool findBestPartialTree();

	// arrange tip nodes
	void arrangeTipNodes();

	// select a tip node
	AobfSearchNode* chooseTipNode(double threshold);

	// get the assigment along the current path from the root
	void getPathAssignment(AobfSearchNode* node, std::vector<int>& vars);

	// Update ancestors values (q-values) along marked paths
	void update(AobfSearchNode* node);

	// Update all ancestors values (q-values)
	void updateAll(AobfSearchNode* node);

	// Update f-values
	void updateF(AobfSearchNode* node);

	// returns a string of the context instantiation (AND context)
	std::string context(int node_type, const std::set<int>& C);

	void check();

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
	AnySAOBF(ProgramOptions* o);
	virtual ~AnySAOBF();

protected:
	double m_prevTimeCheckpoint;
	double m_prevNodeCheckpoint;

};


/* inline functions */

inline AnySAOBF::AnySAOBF(ProgramOptions* o) : Search(o) {
	m_prevTimeCheckpoint = 0;
	m_prevNodeCheckpoint = 0;
	m_globalSearchIndex = 0;
	m_epsilon = o->weight; // default
	m_numSumEvals = 0;
	m_numCacheHits = 0;
	m_stackCount = 0;
	m_stackLimit = o->rotateLimit; // default: 1000 node expansions per stack
	m_rootStack = NULL;
	m_numExpansionsOut = 0;
	m_numExpansionsIn = 0;
}

inline AnySAOBF::~AnySAOBF() {

}

inline void AnySAOBF::_reset(AobbSearchNode* p) {
	assert(p);
	while (m_stacks.size()) {
		delete m_stacks.front();
		m_stacks.pop();
	}
	m_rootStack = new MyStack(NULL);
	m_rootStack->push(p);
	m_stacks.push(m_rootStack);
}


#endif /* SRC_ANY_SAOBF_H_ */
