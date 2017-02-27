/*
 * join_tree.h
 *
 *  Created on: Sep 12, 2013
 *      Author: radu
 */



#ifndef IBM_MAP_JOIN_TREE2_H_
#define IBM_MAP_JOIN_TREE2_H_

#include "base.h"
#include "graph.h"
#include "problem.h"
#include "program_options.h"

#include "mex/Factor.h"

// We assume a Shenoy-Shafer architecture

class JTreeEdge;
class Pseudotree;

// A join tree node
class JTreeNode {

protected:

	int m_id;								// cluster id (from 0)
	std::set<int> m_vars; 					// cluster scope
	std::set<int> m_varsMax;				// MAX variables
	std::set<int> m_varsSum;				// SUM variables
	std::vector<JTreeEdge*> m_edges;		// all incoming edges
	mex::Factor m_factor;					// cluster's original potential
	std::vector<mex::Factor> m_factors;		// list of original factors (MC only)
	std::vector<JTreeNode*> m_children;		// child clusters (links)
	JTreeNode* m_parent;					// parent

	Problem* m_problem;
	ProgramOptions* m_options;

public:

	// initialize
	void init();
	// compute the marginal (new function)
	mex::Factor marginal(int var, const std::vector<int>& assignment);
	mex::Factor marginalIncr(int var, const std::vector<int>& assignment);
	// eliminate all variables in the current node (assume full MAP assignment)
	double upperBound(const std::vector<int>& assignment);
	double upperBoundIncr(const std::vector<int>& assignment);
	// compute the marginal (new function)
	mex::Factor marginal(int var, int ibound, const std::vector<int>& assignment);
	mex::Factor marginalIncr(int var, int ibound, const std::vector<int>& assignment);
	// eliminate all variables in the current node (assume full MAP assignment)
	double upperBound(int ibound, const std::vector<int>& assignment);
	double upperBoundIncr(int ibound, const std::vector<int>& assignment);
	// add a factor
	inline void add(const mex::Factor& f) { m_factors.push_back(f); }
	// get the edges to this node
	inline std::vector<JTreeEdge*>& edges() { return m_edges; }
	// get the node variables
	inline const std::set<int>& vars() const { return m_vars; }
	// check if the cluster contains
	inline bool contains(int var) { return (m_vars.find(var) != m_vars.end()); }
	// get id
	inline int id() { return m_id; }
	// add MAP var
	inline void addMaxVar(int var) { m_varsMax.insert(var); }
	// add SUM var
	inline void addSumVar(int var) { m_varsSum.insert(var); }
	inline const std::set<int>& varsMax() const { return m_varsMax; }
	inline const std::set<int>& varsSum() const { return m_varsSum; }
	inline void setParent(JTreeNode* par) { m_parent = par; }
	inline JTreeNode* parent() { return m_parent; }
	inline std::vector<JTreeNode*>& children() { return m_children; }
	inline void addChild(JTreeNode* ch) { m_children.push_back(ch); }
	inline mex::Factor& factor() { return m_factor; }
	inline const mex::Factor& factor() const { return m_factor; }
	inline std::vector<mex::Factor>& factors() { return m_factors; }
	inline void setProblem(Problem* p, ProgramOptions* po) { m_problem = p; m_options = po; }

	// convert to string
	std::string toString() {
		std::string s;
		s += "( ";
		for (std::set<int>::iterator si = m_vars.begin();
				si != m_vars.end(); ++si) {
			int var = (*si);
			s += static_cast<ostringstream*>( &(ostringstream() << var) )->str();
			s += " ";
		}
		s += ")";
		return s;
	}

public:

	// default constructor
	JTreeNode(int id = NONE);
	JTreeNode(int id, const std::set<int>& vars);

	// destructor
	~JTreeNode();
};

// A join tree edge

class JTreeEdge {

protected:

	JTreeNode* m_node1; 				// the 'from' node
	JTreeNode* m_node2;					// the 'to' node

	std::set<int> m_separator;			// separator

	mex::Factor m_message1to2;				// (forward) message from 1 to 2
	mex::Factor m_message2to1;				// (backward) message from 2 to 1

	std::vector<mex::Factor> m_mcMessage1to2;
	std::vector<mex::Factor> m_mcMessage2to1;
	bool m_mcApproximation;				// use  mini-clustering

	Problem* m_problem;					// problem instance
	ProgramOptions *m_options;

public:

	// send message 1 to 2 (forward)
	void sendMessage1to2(const std::vector<int>& assignment);
	// send message 2 to 1 (backward)
	void sendMessage2to1(const std::vector<int>& assignment);
	// send message 1 to 2 (forward) -- used by the incremental version
	void updateMessage1to2(const std::vector<int>& assignment);
	// send message 2 to 1 (backward) -- used by the incremental version
	void updateMessage2to1(const std::vector<int>& assignment);

	// send message 1 to 2 (forward)
	void sendMessage1to2(int ibound, const std::vector<int>& assignment);
	// send message 2 to 1 (backward)
	void sendMessage2to1(int ibound, const std::vector<int>& assignment);
	// send message 1 to 2 (forward) -- used by the incremental version
	void updateMessage1to2(int ibound, const std::vector<int>& assignment);
	// send message 2 to 1 (backward) -- used by the incremental version
	void updateMessage2to1(int ibound, const std::vector<int>& assignment);

	// returns the separator
	inline const std::set<int>& separator() const { return m_separator; }
	// get node 1
	inline JTreeNode* getNode1() { return m_node1; }
	// get node 2
	inline JTreeNode* getNode2() { return m_node2; }
	// get message 1 to 2
	inline mex::Factor& getMessage1() { return m_message1to2; }
	// get message 2 to 1
	inline mex::Factor& getMessage2() { return m_message2to1; }
	// set message 1 to 2
	inline void setMessage1(mex::Factor& f) { m_message1to2 = f; }
	// set message 2 to 1
	inline void setMessage2(mex::Factor& f) { m_message2to1 = f; }
	// get message 1 to 2
	inline std::vector<mex::Factor>& getVectorMessage1() { return m_mcMessage1to2; }
	// get message 2 to 1
	inline std::vector<mex::Factor>& getVectorMessage2() { return m_mcMessage2to1; }
	// set message 1 to 2
	inline void setVectorMessage1(std::vector<mex::Factor>& f) { m_mcMessage1to2 = f; }
	// set message 2 to 1
	inline void setVectorMessage2(std::vector<mex::Factor>& f) { m_mcMessage2to1 = f; }

	// clean the edge
	void clean();

	// set problem instance
	inline void setProblem(Problem* p, ProgramOptions* po, bool mc) {
		m_problem = p; m_options = po; m_mcApproximation = mc;
	}

public:

	// constructor
	JTreeEdge(JTreeNode* s1, JTreeNode* s2, const std::set<int>& sep);

	// destructor
	~JTreeEdge();
};


struct PathEdge {
	int from;				// 'from' cluster id
	int to;					// 'to' cluster id
	JTreeEdge* edge;			// corresponnding join tree edge
	int activeMessage;		// 1 or 2

	PathEdge(int f, int t) : from(f), to(t), edge(NULL), activeMessage(NONE) {};
	void updateMessage(const std::vector<int>& assignment);
	void updateMessage(int ibound, const std::vector<int>& assignment);
};

struct Change {
	JTreeEdge* edge;
	int message;
	mex::Factor potential;
	std::vector<mex::Factor> potentialApprox;
	Change(JTreeEdge* e, int m, mex::Factor& f) : edge(e), message(m), potential(f) {};
	Change(JTreeEdge* e, int m, std::vector<mex::Factor>& f) : edge(e), message(m), potentialApprox(f) {};
};


// The join tree
class JTree {

protected:

	JTreeNode* m_root;							// root of the join tree
	std::vector<JTreeNode*> m_clusters;			// nodes in the join tree
	std::vector<JTreeEdge*> m_edges;			// edges in the join tree
	std::vector<JTreeNode*> m_var2cluster;		// maps variables to nodes in the join tree
	std::map<int, JTreeNode*> m_id2cluster;		// maps node ids to nodes in the join tree
	std::map<int, JTreeNode*> m_map2cluster;	// maps a MAP variable to its search cluster
	std::vector<JTreeEdge*> m_messageOrder;		// message order (full)
	std::set<int> m_varsMax;					// map variables
	std::set<int> m_sumVars;					// sum variable
	double m_upperBound;						// upper bound of (partial) MAP
	int m_treewidth;							// treewidth
	size_t m_clusterSize;						// max cluster size
	size_t m_separatorSize;						// max separator size
	std::list<int> m_searchOrder; 				// MAP search order (static)
	std::vector<JTreeNode*> m_searchClusters;	// MAP search clusters (static)
	std::map<int, std::vector<PathEdge*> > m_paths;// tree paths between nodes containing MAP variables
	Problem* m_problem;							// Problem instance
	ProgramOptions* m_options;					// Program options
	std::map<int, mex::Factor> m_maxmarginals;	// MAP marginals

protected:

	// find the post order for searching the MAP variables
	void findPostOrder();

	// find a path between two nodes in the join tree
	void findPath(JTreeNode* from, JTreeNode* to, std::vector<PathEdge*>& path);

public:

	// build a join tree (from an elimination order)
	void build(const Graph& G, const std::vector<int>& elim, const bool mc = false);

	// message passing given a partial assignment to the MAP variables
	void propagate(const std::vector<int>& assignment);

	// initial message passing for mini-clustering
	void propagate(int ibound, const std::vector<int>& assignment);

	// clean previous messages
	void clean();

	// get the marginal of a MAP variable
	mex::Factor& marginal(int var) {
		return m_maxmarginals[var];
	}

	// get current MAP assignment value (if any)
	double upperBound() {
		return m_upperBound;
	}

	// set the MAP search order (static)
	void setSearchOrder(const std::list<int>& o) {
		m_searchOrder = o;
	}

	// get the MAP search order (static)
	const std::list<int>& getSearchOrder() const {
		return m_searchOrder;
	}

	// get the MAP search clusters (static)
	const std::vector<JTreeNode*>& getSearchClusters() const {
		return m_searchClusters;
	}

	// returns the search cluster of a MAP variable
	JTreeNode* getSearchCluster(int var) {
		std::map<int, JTreeNode*>::iterator mi = m_map2cluster.find(var);
		assert(mi != m_map2cluster.end());
		return mi->second;
	}

	// get the path associated with a MAP variable
	std::vector<PathEdge*>& getPath(int var) {
		std::map<int, std::vector<PathEdge*> >::iterator mi = m_paths.find(var);
		assert( mi != m_paths.end() );
		return mi->second;
	}

	inline void setUpperBound(double upbo) {
		m_upperBound = upbo;
	}

	// dump into Graphviz file format
	void toDot(const std::string& fileName);

private:

	bool dominates(const std::set<int>& A, const std::set<int>& B);
	void buildJunctionTree(bool mc = false); // call Kruskal's

public:

	JTree(Problem* p, ProgramOptions* po);
	virtual ~JTree();
};


/* inline definitions */

inline JTree::JTree(Problem* p, ProgramOptions* po) {
	m_root = NULL;
	m_treewidth = NONE;
	m_problem = p;
	m_options = po;
	m_upperBound = ELEM_NAN;
	m_clusterSize = 0;
	m_separatorSize = 0;
}

inline JTree::~JTree() {

	// clean up
	for (size_t i = 0; i < m_edges.size(); ++i) {
		if (m_edges[i]) delete m_edges[i];
	}

	for (size_t i = 0; i < m_clusters.size(); ++i) {
		if (m_clusters[i]) delete m_clusters[i];
	}

	for (std::map<int, std::vector<PathEdge*> >::iterator mi = m_paths.begin();
			mi != m_paths.end(); ++mi) {
		std::vector<PathEdge*>& p = mi->second;
		for (std::vector<PathEdge*>::iterator it = p.begin(); it != p.end(); ++it) {
			delete (*it);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
inline JTreeNode::JTreeNode(int id) {
	m_id = id;
	m_parent = NULL;
	m_problem = NULL;
	m_options = NULL;
}

inline JTreeNode::JTreeNode(int id, const std::set<int>& vars) {
	m_id = id;
	m_vars = vars;
	m_parent = NULL;
	m_problem = NULL;
	m_options = NULL;
}

inline JTreeNode::~JTreeNode() {

}

// initialize a join tree node
inline void JTreeNode::init() {

}

///////////////////////////////////////////////////////////////////////////////
inline JTreeEdge::JTreeEdge(JTreeNode* s1, JTreeNode* s2,
		const std::set<int>& sep) {

	m_node1 = s1;
	m_node2 = s2;
	m_separator = sep;
	m_problem = NULL;
	m_options = NULL;
	m_mcApproximation = false;
}

inline void JTreeEdge::clean() {
	mex::Factor empty;
	m_message1to2.swap(empty);
	m_message2to1.swap(empty);
	m_mcMessage1to2.clear();
	m_mcMessage2to1.clear();
}

inline JTreeEdge::~JTreeEdge() {

}



#endif /* JOIN_TREE_H_ */
