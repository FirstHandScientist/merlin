/*
 * search_node.h
 *
 *  Created on: Mar 28, 2013
 *      Author: radu
 */

#ifndef IBM_MAP_SEARCH_NODE_H_
#define IBM_MAP_SEARCH_NODE_H_

#include "base.h"

#include <boost/variant.hpp>

class AobfSearchNode;

/* some constants for aggregating the boolean flags */
#define FLAG_FRINGE 1 // node is on the fringe (ie, not yet expanded)
#define FLAG_SOLVED 2 // node is solved
#define FLAG_LOOKAHEAD 4 // node is looked ahead
#define FLAG_TERMINAL 8 // node is terminal
#define FLAG_EXPANDED 16 // node is expanded
#define FLAG_VISITED 32 // node is visited
#define FLAG_REPAIRED 64 // node is repaired
#define FLAG_DEADEND 128 // node is deadend (infinity cost or arc weight)
#define FLAG_FEASIBLE 256 // node is feasible

/* A search node in the context-minimal AND/OR graph */
class AobfSearchNode {
protected:

	// Flags
	unsigned int m_flags;

	// Depth in the search space
	int m_depth;

	// Node value (as in cost): lower bound
	double m_nodeValue;

	// Heuristic estimate of the node's value
	double m_heurValue;

	// Upper bound (a solution found below this node)
	double m_nodeUpperBound;

	// Pointers to parents
	std::list<AobfSearchNode*> m_parents;

	// Pointers to children
	std::list<AobfSearchNode*> m_children;

	// Node variable
	int m_var;

	// Node index
	int m_index;

	// Current parent
	AobfSearchNode* m_currentParent;

	// Hash key (context)
	size_t m_key;

	// Path assignment
	std::vector<std::pair<int, int> > m_path;

public:

	// returns the node type (NODE_AND, NODE_OR)
	virtual int getType() const = 0;

	// returns the node variable (as in X_i)
	inline int getVar() const {
		return m_var;
	}

	// returns the node value assignment (as in <X_i,a>)
	virtual int getVal() const = 0;

	// set the node value (as in cost)
	inline void setValue(double f) {
		m_nodeValue = f;
	}

	// returns the node value (as in cost)
	inline double getValue() const {
		return m_nodeValue;
	}

	// set the node heuristic value
	inline void setHeur(double h) {
		m_heurValue = h;
	}

	// returns the node heuristic value
	inline double getHeur() const {
		return m_heurValue;
	}

	// set the node value (as in cost)
	inline void setUpperBound(double ub) {
		m_nodeUpperBound = ub;
	}

	// returns the node value (as in cost)
	inline double getUpperBound() const {
		return m_nodeUpperBound;
	}

	// sets the current parent in the current partial solution tree
	inline void setCurrentParent(AobfSearchNode* node) {
		m_currentParent = node;
	}

	// returns the current parent in the current partial solution tree
	inline AobfSearchNode* getCurrentParent() {
		return m_currentParent;
	}

	// set the index
	inline void setIndex(int i) {
		m_index = i;
	}

	// returns the node index
	inline int getIndex() const {
		return m_index;
	}

	// increments the index
	inline void incIndex() {
		++m_index;
	}

	// decrements the index
	inline void decIndex() {
		--m_index;
	}

	// returns the depth of the node in the search graph
	inline int getDepth() const {
		return m_depth;
	}

	// add a child to the list
	inline void addChild(AobfSearchNode* node) {
		m_children.push_back(node);
	}

	// add a parent to the list
	inline void addParent(AobfSearchNode* node) {
		m_parents.push_back(node);
	}

	// returns the children list (const)
	inline const std::list<AobfSearchNode*>& getChildren() const {
		return m_children;
	}

	// returns the children list
	inline std::list<AobfSearchNode*>& getChildren() {
		return m_children;
	}

	// returns the parents list (const)
	inline const std::list<AobfSearchNode*>& getParents() const {
		return m_parents;
	}

	// returns the parents list
	inline std::list<AobfSearchNode*>& getParents() {
		return m_parents;
	}

	// return the path
	inline std::vector<std::pair<int, int> >& getPath() {
		return m_path;
	}

	// set best child
	virtual void setBestChild(AobfSearchNode* node) = 0;

	// returns the best child
	virtual AobfSearchNode* getBestChild() = 0;

	// set feasible child
	virtual void setFeasibleChild(AobfSearchNode* node) = 0;

	// returns the feasible child
	virtual AobfSearchNode* getFeasibleChild() = 0;

	// check if 'node' belongs to the children list
	bool hasChild(AobfSearchNode* node) const;

	// erase 'node' from the children list
	void eraseChild(AobfSearchNode* node);

	// check if 'node' belongs to the parents list
	bool hasParent(AobfSearchNode* node) const;

	// erase 'node' from the parents list
	void eraseParent(AobfSearchNode* node);

	// set the heuristic cache
	virtual void setHeurCache(double*) = 0;

	// returns the heuristic cache at OR nodes
	virtual double* getHeurCache() const = 0;

	// clear the cache
	virtual void clearHeurCache() = 0;

	// set the arc weight
	virtual void setWeight(int, double) = 0;

	// returns an arc weight
	virtual double getWeight(int) const = 0;

	// set the hash key
	inline void setKey(size_t k) {
		m_key = k;
	}

	// returns the hash key
	inline size_t getKey() const {
		return m_key;
	}

	// set/return boolean flags
	inline void setFringe(bool flag) {
		if (flag == true) { // set the bit to 1
			m_flags |= FLAG_FRINGE;
		} else { // set the bit to 0
			m_flags &= ~(FLAG_FRINGE);
		}
	}
	inline bool isFringe() const {
		return m_flags & FLAG_FRINGE;
	}
	inline void setSolved(bool flag) {
		if (flag == true) { // set the bit to 1
			m_flags |= FLAG_SOLVED;
		} else { // set the bit to 0
			m_flags &= ~(FLAG_SOLVED);
		}
	}
	inline bool isSolved() const {
		return m_flags & FLAG_SOLVED;
	}
	inline void setLookedAhead(bool flag) {
		if (flag == true) { // set the bit to 1
			m_flags |= FLAG_LOOKAHEAD;
		} else { // set the bit to 0
			m_flags &= ~(FLAG_LOOKAHEAD);
		}
	}
	inline bool isLookedAhead() const {
		return m_flags & FLAG_LOOKAHEAD;
	}
	inline void setTerminal(bool flag) {
		if (flag == true) { // set the bit to 1
			m_flags |= FLAG_TERMINAL;
		} else { // set the bit to 0
			m_flags &= ~(FLAG_TERMINAL);
		}
	}
	inline bool isTerminal() const {
		return m_flags & FLAG_TERMINAL;
	}
	inline void setExpanded(bool flag) {
		if (flag == true) { // set the bit to 1
			m_flags |= FLAG_EXPANDED;
		} else { // set the bit to 0
			m_flags &= ~(FLAG_EXPANDED);
		}
	}
	inline bool isExpanded() const {
		return m_flags & FLAG_EXPANDED;
	}
	inline void setVisited(bool flag) {
		if (flag == true) { // set the bit to 1
			m_flags |= FLAG_VISITED;
		} else { // set the bit to 0
			m_flags &= ~(FLAG_VISITED);
		}
	}
	inline bool isVisited() const {
		return m_flags & FLAG_VISITED;
	}
	inline void setRepaired(bool flag) {
		if (flag == true) { // set the bit to 1
			m_flags |= FLAG_REPAIRED;
		} else { // set the bit to 0
			m_flags &= ~(FLAG_REPAIRED);
		}
	}
	inline bool isRepaired() const {
		return m_flags & FLAG_REPAIRED;
	}
	inline void setDeadend(bool flag) {
		if (flag == true) { // set the bit to 1
			m_flags |= FLAG_DEADEND;
		} else { // set the bit to 0
			m_flags &= ~(FLAG_DEADEND);
		}
	}
	inline bool isDeadend() const {
		return m_flags & FLAG_DEADEND;
	}
	inline void setFeasible(bool flag) {
		if (flag == true) { // set the bit to 1
			m_flags |= FLAG_FEASIBLE;
		} else { // set the bit to 0
			m_flags &= ~(FLAG_FEASIBLE);
		}
	}
	inline bool isFeasible() const {
		return m_flags & FLAG_FEASIBLE;
	}

	// Converts the node info to a string
	virtual std::string toString() = 0;

public:

	AobfSearchNode();

	virtual ~AobfSearchNode();
};


class AobfSearchNodeAND : public AobfSearchNode {
protected:

	// Value assignment
	int m_val;

public:

	// Converts the node info to a string
	std::string toString();

	inline int getType() const {
		return NODE_AND;
	}

	inline int getVal() const {
		return m_val;
	}

	// set the arc weight
	void setWeight(int, double) { assert(false); };
	// returns an arc weight
	double getWeight(int) const { assert(false); return ELEM_NAN; };

	// cache manipulation
	void setHeurCache(double* d) { assert(false); }
	double* getHeurCache() const { assert(false); return NULL; }
	void clearHeurCache() { assert(false); };

	// get the best child
	inline AobfSearchNode* getBestChild() {
		assert(false); return NULL; // not used for AND nodes
	}

	// set the best child
	inline void setBestChild(AobfSearchNode* best) {
		assert(false); // not used for AND nodes
	}

	// set feasible child
	inline void setFeasibleChild(AobfSearchNode* node) {
		assert(false); // not used for AND nodes
	}

	// returns the feasible child
	inline AobfSearchNode* getFeasibleChild() {
		assert(false); return NULL; // not used for AND nodes
	}

public:

	AobfSearchNodeAND(int var, int val, int depth);
	~AobfSearchNodeAND();
};

class AobfSearchNodeOR : public AobfSearchNode {

protected:

	// Best child
	AobfSearchNode* m_bestChild;

	// Feasible child
	AobfSearchNode* m_feasibleChild;

	// Stores the precomputed heuristic values of the AND children
	double* m_heurCache;

public:

	// cache manipulation
	void setHeurCache(double* d) {
		m_heurCache = d;
	}
	double* getHeurCache() const {
		return m_heurCache;
	}
	void clearHeurCache() {
		if (m_heurCache) delete[] m_heurCache;
	}

	// Converts the node info to a string
	std::string toString();

	// returns node type
	inline int getType() const {
		return NODE_OR;
	}

	// returns node value assignment
	inline int getVal() const {
		assert(false); // no value for OR nodes
		return NONE;
	}

	// set the arc weight
	inline void setWeight(int val, double w) {
		assert(m_heurCache != NULL);
		m_heurCache[2*val+1] = w;
	}

	// returns an arc weight
	inline double getWeight(int val) const {
		assert(m_heurCache != NULL);
		return m_heurCache[2*val+1];
	}

	// get the best child
	inline AobfSearchNode* getBestChild() {
		return m_bestChild;
	}

	// set the best child
	inline void setBestChild(AobfSearchNode* best) {
		m_bestChild = best;
	}

	// get the feasible child
	inline AobfSearchNode* getFeasibleChild() {
		return m_feasibleChild;
	}

	// set the best child
	inline void setFeasibleChild(AobfSearchNode* feas) {
		m_feasibleChild = feas;
	}

public:

	AobfSearchNodeOR(int var, int depth);
	~AobfSearchNodeOR();
};

/* inline functions*/

// constructor
inline AobfSearchNode::AobfSearchNode() :
		m_flags(0),
		m_depth(NONE),
		m_nodeValue(ELEM_NAN),
		m_heurValue(ELEM_NAN),
		m_nodeUpperBound(ELEM_INF),
		m_var(NONE),
		m_index(0),
		m_currentParent(NULL),
		m_key(0) {

}

// destructor
inline AobfSearchNode::~AobfSearchNode() {
	m_parents.clear();
	m_children.clear();
}

// check if 'node' belongs to the children list
inline bool AobfSearchNode::hasChild(AobfSearchNode* node) const {
	std::list<AobfSearchNode*>::const_iterator it = m_children.begin();
	for (; it != m_children.end(); ++it) {
		AobfSearchNode* child = (*it);
		if (child == node) {
			return true;
		}
	}

	return false;
}

// erase 'node' from the children list
inline void AobfSearchNode::eraseChild(AobfSearchNode* node) {
	std::list<AobfSearchNode*>::iterator it = m_children.begin();
	for (; it != m_children.end();) {
		AobfSearchNode* child = (*it);
		if (child == node) {
			m_children.erase(it);
			break;
		} else {
			++it;
		}
	}
}

// check if 'node' belongs to the parents list
inline bool AobfSearchNode::hasParent(AobfSearchNode* node) const {
	std::list<AobfSearchNode*>::const_iterator it = m_parents.begin();
	for (; it != m_parents.end(); ++it) {
		AobfSearchNode* parent = (*it);
		if (parent == node) {
			return true;
		}
	}

	return false;
}

// erase 'node' from the parents list
inline void AobfSearchNode::eraseParent(AobfSearchNode* node) {
	std::list<AobfSearchNode*>::iterator it = m_parents.begin();
	for (; it != m_parents.end(); ) {
		AobfSearchNode* parent = (*it);
		if (parent == node) {
			m_parents.erase(it);
			break;
		} else {
			++it;
		}
	}
}

inline AobfSearchNodeAND::AobfSearchNodeAND(int var, int val, int depth) : AobfSearchNode() {
	m_var = var;
	m_val = val;
	m_depth = depth;
}

inline AobfSearchNodeAND::~AobfSearchNodeAND() {

}

inline AobfSearchNodeOR::AobfSearchNodeOR(int var, int depth) : AobfSearchNode() {
	m_var = var;
	m_depth = depth;
	m_bestChild = NULL;
	m_feasibleChild = NULL;
	m_heurCache = NULL;
}

inline AobfSearchNodeOR::~AobfSearchNodeOR() {
	clearHeurCache();
}

// convert a node to a string
inline std::string AobfSearchNodeOR::toString() {

	std::ostringstream oss;
	oss << "OR node: (x" << getVar() << ")"
		<< ", h = " << getHeur()
		<< ", q = " << getValue()
		<< ", ub = " << getUpperBound()
		<< ", depth = " << getDepth()
		<< ", weights { "
		;
	std::list<AobfSearchNode*>::const_iterator it = getChildren().begin();
	for (; it != getChildren().end(); ++it) {
		int val = (*it)->getVal();
		oss << val << ":" << getWeight(val) << " ";
	}

	oss << "}, children { ";
	std::list<AobfSearchNode*>::const_iterator ci = getChildren().begin();
	for (; ci != getChildren().end(); ++ci) {
		int val = (*ci)->getVal();
		oss << val << ":" << (*ci)->getHeur() << " ";
	}

	oss << "}";
	oss << ", expanded = " << (isExpanded() ? "YES" : "NO");
	oss << ", solved = " << (isSolved() ? "YES" : "NO");
	//oss << ", looked = " << (isLookedAhead() ? "YES" : "NO");

	return oss.str();
}

// convert a node to a string
inline std::string AobfSearchNodeAND::toString() {

	std::ostringstream oss;

	oss << "AND node: (x" << getVar() << "," << getVal() << ")"
		<< ", h = " << getHeur()
		<< ", q = " << getValue()
		<< ", ub = " << getUpperBound()
		<< ", depth = " << getDepth()
		<< ", children { "
		;
	std::list<AobfSearchNode*>::const_iterator it = getChildren().begin();
	for (; it != getChildren().end(); ++it) {
		int var = (*it)->getVar();
		oss << var << ":" << (*it)->getHeur() << " ";
	}

	oss << "}";
	oss << ", expanded = " << (isExpanded() ? "YES" : "NO");
	oss << ", solved = " << (isSolved() ? "YES" : "NO");

	return oss.str();
}

#endif /* SEARCH_NODE_H_ */
