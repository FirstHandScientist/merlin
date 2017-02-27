/*
 * pseudo_tree.h
 *
 *  Created on: Mar 22, 2013
 *      Author: radu
 */

#ifndef IBM_ANYTIME_PSEUDO_TREE_H_
#define IBM_ANYTIME_PSEUDO_TREE_H_

#include "base.h"
#include "graph.h"
#include "function.h"
#include "problem.h"

/* forward declaration */
class PseudotreeNode;

/* Represents the guiding pseudo tree */
class Pseudotree {

	friend class PseudotreeNode;

protected:

	// Subproblem ordering heuristic
	int m_subOrder;

	// Height of the pseudo tree
	int m_height;

	// Height of a conditioned pseudo tree
//	int m_heightConditioned;

	// Induced width of the pseudo tree
	int m_width;

	// Width of a conditioned pseudo tree
//	int m_widthConditioned;

	// Number of disconnected components
	int m_components;

	// Size of the pseudo tree
	int m_size;

	// Size of a conditioned pseudo tree
//	int m_sizeConditioned;

	// Problem instance
	Problem* m_problem;

	// Root node of the pseudo tree
	PseudotreeNode* m_root;

	// for minfill
//	std::vector<nCost> m_initialScores;

	// Nodes of the pseudo tree
	std::vector<PseudotreeNode*> m_nodes;

	// Elimination order
	std::vector<int> m_elimOrder;

	// Keeps track of SUM roots
	std::vector<int> m_sumRoots;

public:
	/* builds the pseudo tree according to constrained order 'elim' */
	void build(Graph G, const std::vector<int>& elim, const int cachelimit = NONE);

	/* builds the pseudo tree according to constrained order 'elim', such that
	 * the MAP variables form a chain at the top of the pseudo tree */
	void forceMapChain(const std::list<int>& mapSearchOrder, const int cachelimit = NONE);

	bool isSumRoot(int var) {
		return m_sumRoots[var];
	}

	/* returns the width of the (sub) pseudo tree */
	int getWidth() const {
		return m_width;
	}

	/* returns the height of the (sub) pseudo tree */
	int getHeight() const {
		return m_height;
	}

	/* returns the number of nodes in the (sub) pseudo tree */
	int getSize() const {
		return m_size;
	}

	/* returns the number of nodes in the *full* pseudo tree */
	size_t getN() const {
		return m_nodes.size();
	}

	/* returns the number of components */
	int getComponents() const {
		return m_components;
	}

	/* returns the root node */
	PseudotreeNode* getRoot() const {
		return m_root;
	}

	/* returns the pseudotree node for variable i */
	PseudotreeNode* getNode(int i) const;

	/* returns the elimination order */
	const vector<int>& getElimOrder() const {
		return m_elimOrder;
	}

	/* augments the pseudo tree with function information */
	void resetFunctionInfo(const vector<Function*>& fns);

	/* returns the list of functions associated with a node */
	const list<Function*>& getFunctions(int i) const;

	/* select top M variables (start pseudo tree) */
	void bfs(std::vector<int>& start);

	/* return a dfs order of the variables in the pseudo tree */
	void dfs(std::list<int>& order, const bool mapOnly = true);

private:
	/* recursive helper function for outputToFile(...) below */
	void outputToFileNode(const PseudotreeNode*, ostringstream&) const;

public:
	/* outputs the pseudo tree structure to a file in ASCII format, to
	 * be used as input to a plotting script, for instance */
	void outputToFile(string basename) const;
	void outputToDot(const char *dotfile );

private:
	void drawNodesForDot( std::ofstream& outfile );
	void drawLinesForDot( std::ofstream& outfile );

protected:
	/* creates a new node in the PT for variable i, with context N. Also makes sure
	 * existing roots are checked and connected appropriately. */
	void insertNewNode(const int i, const std::set<int>& N,
			std::list<PseudotreeNode*>& roots);

	// "clears out" the pseudo tree by resetting parameters and removing the node structure
	void reset();

public:
	// constructors
	Pseudotree(Problem* p, int subOrder);
	Pseudotree(const Pseudotree& pt); // copy constructor
	~Pseudotree();
};

/* Represents a single problem variable in the pseudotree */
class PseudotreeNode {

protected:
	int m_var; // The node variable
	int m_depth; // The node's depth in the tree
	int m_subHeight; // The node's height in the tree (distance from farthest child)
	int m_subWidth; // max. width in subproblem
	PseudotreeNode* m_parent; // The parent node
	Pseudotree* m_tree;
	std::set<int> m_subproblemVars; // The variables in the subproblem (including self)
	std::vector<int> m_subproblemVarMap; // Maps variables to their index in subprob assignment
	std::set<int> m_context; // The node's full OR context (!doesn't include own variable!)
	std::set<int> m_contextAnd; // The node's full AND context (includes own variable)
	std::list<Function*> m_functions; // The functions that will be fully instantiated at this point
	std::vector<PseudotreeNode*> m_children; // The node's children

	std::set<int> m_cacheContext; // The (possibly smaller) context for (adaptive) caching
	std::list<int> m_cacheResetList; // List of var's whose cache tables need to be reset when this
	                              // var's search node is expanded (for adaptive caching)

public:

	void setParent(PseudotreeNode* p) {
		m_parent = p;
	}
	PseudotreeNode* getParent() const {
		return m_parent;
	}

	void addChild(PseudotreeNode* p) {
		m_children.push_back(p);
	}
	void setChild(PseudotreeNode* p) {
		m_children.clear();
		m_children.push_back(p);
	}
	const vector<PseudotreeNode*>& getChildren() const {
		return m_children;
	}
	void orderChildren(int subOrder);

	static bool compGreater(PseudotreeNode* a, PseudotreeNode* b);
	static bool compLess(PseudotreeNode* a, PseudotreeNode* b);

	void setFullContext(const set<int>& c) {
		m_context = c;
	}
	void addFullContext(int v) {
		m_context.insert(v);
	}
	const std::set<int>& getFullContext() const {
		return m_context;
	}

	void setAndContext(const std::set<int>& c) {
		m_contextAnd = c;
	}
	void addAndContext(int v) {
		m_contextAnd.insert(v);
	}
	const std::set<int>& getAndContext() const {
		return m_contextAnd;
	}

	void setCacheContext(const set<int>& c) { m_cacheContext = c; }
	void addCacheContext(int i) { m_cacheContext.insert(i); }
	const set<int>& getCacheContext() const { return m_cacheContext; }

	void setCacheReset(const list<int>& l) { m_cacheResetList = l; }
	void addCacheReset(int i) { m_cacheResetList.push_back(i); }
	const list<int>& getCacheReset() const { return m_cacheResetList; }


	void addFunction(Function* f) {
		m_functions.push_back(f);
	}
	void setFunctions(const std::list<Function*>& l) {
		m_functions = l;
	}
	void resetFunctions() {
		m_functions.clear();
	}
	const std::list<Function*>& getFunctions() const {
		return m_functions;
	}

	int getVar() const {
		return m_var;
	}
	int getDepth() const {
		return m_depth;
	}
	int getSubHeight() const {
		return m_subHeight;
	}
	int getSubWidth() const {
		return m_subWidth;
	}

	size_t getSubprobSize() const {
		return m_subproblemVars.size();
	}
	const std::set<int>& getSubprobVars() const {
		return m_subproblemVars;
	}
	const std::vector<int>& getSubprobVarMap() const {
		return m_subproblemVarMap;
	}
	void setSubprobVarMap(const std::vector<int>& map) {
		m_subproblemVarMap = map;
	}

public:
	const std::set<int>& updateSubprobVars(int numVars);
	int updateDepthHeight(int d);
	int updateSubWidth();

	void dumpFullContext(std::ostream& os);
	void dumpAndContext(std::ostream& os);

public:
	/* creates new pseudo tree node for variable v with context s */
	PseudotreeNode(Pseudotree* t, int v, const std::set<int>& s);
	~PseudotreeNode() {
	}
};

/* Inline definitions */

inline PseudotreeNode* Pseudotree::getNode(int i) const {
	assert(i<(int)m_nodes.size());
	return m_nodes[i];
}

inline const list<Function*>& Pseudotree::getFunctions(int i) const {
	assert(i < (int) m_nodes.size());
	return m_nodes[i]->getFunctions();
}

inline void Pseudotree::reset() {
	assert(m_problem);
	m_height = UNKNOWN;
	m_width = UNKNOWN;
	m_components = 0;
	m_size = UNKNOWN;
	m_root = NULL;
	for (vector<PseudotreeNode*>::iterator it = m_nodes.begin();
			it != m_nodes.end(); ++it) {
		delete *it;
		*it = NULL;
	}
	m_nodes.clear();
	m_nodes.reserve(m_problem->getN() + 1); // +1 for dummy variable
	m_nodes.resize(m_problem->getN(), NULL);
	m_size = m_problem->getN();
}

inline Pseudotree::Pseudotree(Problem* p, int subOrder) :
		m_subOrder(subOrder), m_height(UNKNOWN), /*m_heightConditioned(UNKNOWN),*/
		m_width(UNKNOWN), /*m_widthConditioned(UNKNOWN),*/ m_components(0),
		m_size(UNKNOWN), /*m_sizeConditioned(UNKNOWN),*/ m_problem(p), m_root(NULL) {
	assert(p);
	m_nodes.reserve(p->getN() + 1); // +1 for bogus node
	m_nodes.resize(p->getN(), NULL);
	m_size = p->getN();
}

inline Pseudotree::~Pseudotree() {
	for (std::vector<PseudotreeNode*>::iterator it = m_nodes.begin();
			it != m_nodes.end(); ++it) {
		if (*it) {
			delete *it;
		}
	}
}

/* PseudotreeNode inlines */

/* Constructor */
inline PseudotreeNode::PseudotreeNode(Pseudotree* t, int v, const set<int>& s) :
		m_var(v), m_depth(UNKNOWN), m_subHeight(UNKNOWN), m_subWidth(UNKNOWN),
		m_parent(NULL), m_tree(t), m_context(s) {
}

// dump OR context
inline void PseudotreeNode::dumpFullContext(std::ostream& os) {
	os << m_var << ": [ ";
	std::copy(m_context.begin(), m_context.end(),
			std::ostream_iterator<int>(os, " "));
	os << "]" << std::endl;
}

// dump AND context
inline void PseudotreeNode::dumpAndContext(std::ostream& os) {
	os << m_var << ": [ ";
	std::copy(m_contextAnd.begin(), m_contextAnd.end(),
			std::ostream_iterator<int>(os, " "));
	os << "]" << std::endl;
}

#endif /* PSEUDO_TREE_H_ */
