/*
 * SearchNode.h
 *
 *  Copyright (C) 2008-2012 Lars Otten
 *  This file is part of DAOOPT.
 *
 *  DAOOPT is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  DAOOPT is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with DAOOPT.  If not, see <http://www.gnu.org/licenses/>.
 *  
 *  Created on: Oct 9, 2008
 *      Author: Lars Otten <lotten@ics.uci.edu>
 */

#ifndef IBM_MAP_SEARCHNODE_H_
#define IBM_MAP_SEARCHNODE_H_

#include "base.h"
#include "utils.h"

class Problem;
class PseudotreeNode;

/* some constants for aggregating the boolean flags */
#define FLAG_LEAF 1 // node is a leaf node
#define FLAG_CACHABLE 2 // node is candidate for caching
#define FLAG_PRUNED 4 // subproblem below was pruned
#define FLAG_NOTOPT 8 // subproblem possibly not optimally solved (-> don't cache)
#define FLAG_EXPANDED 16 // node is expanded

class AobbSearchNode;

/* data types to store child pointers */
typedef AobbSearchNode* NodeP;
typedef NodeP* CHILDLIST;

class AobbSearchNode {
protected:
	unsigned char m_flags;             // for the boolean flags
	int m_depth;                       // depth in search space
	AobbSearchNode* m_parent;          // pointer to the parent
	double m_nodeValue;                // node value (as in cost); upper bound
	double m_heurValue;                // heuristic estimate of the node's value; h-value
	double m_heurUpdated;				// updated heuristic value; q-value

	CHILDLIST m_children;              // Child nodes
	size_t m_childCountFull;      	   // Number of total child nodes (initial count)
	size_t m_childCountAct;            // Number of remaining active child nodes

	std::list<std::pair<int,int> > m_changes;

#ifndef NO_ASSIGNMENT
	vector<val_t> m_optAssignment;     // stores the optimal solution to the subproblem
#endif

public:
	virtual int getType() const = 0;
	virtual int getVar() const = 0;
	virtual val_t getVal() const = 0;

	virtual void setValue(double) = 0;
	virtual double getValue() const = 0;
	virtual double getLabel() const = 0;
	virtual void addSubSolved(double) = 0;
	virtual double getSubSolved() const = 0;
	virtual void setSubSolved(double) = 0;
	void setHeur(double d) {
		m_heurValue = d;
	}
	// the first one is overridden in SearchNodeOR, the second one isn't
	virtual double getHeur() const {
		return m_heurValue;
	}
	void setHeurUpdated(double d) {
		m_heurUpdated = d;
	}
	double getHeurUpdated() const {
		return m_heurUpdated;
	}

	virtual void setCacheContext(const context_t&) = 0;
	virtual const context_t& getCacheContext() const = 0;

	int getDepth() const {
		return m_depth;
	}
	void setDepth(int d) {
		m_depth = d;
	}

	AobbSearchNode* getParent() const {
		return m_parent;
	}

	virtual double getPathAssignment(std::vector<int>&) = 0;
	void setChild(AobbSearchNode*);
	void addChildren(const vector<AobbSearchNode*>&);
	NodeP* getChildren() const {
		return m_children;
	}
	size_t getChildCountFull() const {
		return m_childCountFull;
	}
	size_t getChildCountAct() const {
		return m_childCountAct;
	}

	bool hasChild(AobbSearchNode* node) const;
	void eraseChild(AobbSearchNode* node);
	void clearChildren();

	void setLeaf() {
		m_flags |= FLAG_LEAF;
	}
	bool isLeaf() const {
		return m_flags & FLAG_LEAF;
	}
	void setCachable() {
		m_flags |= FLAG_CACHABLE;
	}
	bool isCachable() const {
		return m_flags & FLAG_CACHABLE;
	}
	void setPruned() {
		m_flags |= FLAG_PRUNED;
	}
	bool isPruned() const {
		return m_flags & FLAG_PRUNED;
	}
	void setNotOpt() {
		m_flags |= FLAG_NOTOPT;
	}
	bool isNotOpt() const {
		return m_flags & FLAG_NOTOPT;
	}
	void setExpanded() {
		m_flags |= FLAG_EXPANDED;
	}
	bool isExpanded() const {
		return m_flags & FLAG_EXPANDED;
	}

	virtual void setHeurCache(double* d) = 0;
	virtual double* getHeurCache() const = 0;
	virtual void clearHeurCache() = 0;

	inline std::list<std::pair<int,int> >& changes(void) {
		return m_changes;
	};

#ifndef NO_ASSIGNMENT
	vector<val_t>& getOptAssig() {return m_optAssignment;}
	void setOptAssig(const std::vector<val_t>& assign) {m_optAssignment = assign;}
	void clearOptAssig() {vector<val_t> v; m_optAssignment.swap(v);}
#endif

protected:
	AobbSearchNode(AobbSearchNode* parent);

public:
	static bool heurLess(const AobbSearchNode* a, const AobbSearchNode* b);
	static bool heurGreater(const AobbSearchNode* a, const AobbSearchNode* b);
	virtual string toString() = 0;

public:
	virtual ~AobbSearchNode();
};

class AobbSearchNodeAND: public AobbSearchNode {
protected:
	val_t m_val;          // Node value, assignment to OR parent variable
	double m_nodeLabel;   // Label of arc <X_i,a>, i.e. instantiated function costs
	double m_subSolved;   // Saves solutions of optimally solved subproblems, so that
						  // their nodes can be deleted
	static context_t emptyCtxt;
	static std::list<std::pair<double, double> > emptyPSTList;

public:
	int getType() const {
		return NODE_AND;
	}
	int getVar() const {
		assert(m_parent);
		return m_parent->getVar();
	}
	val_t getVal() const {
		return m_val;
	}

	void setValue(double d) {
		m_nodeValue = d;
	}
	double getValue() const {
		return m_nodeValue;
	}
	double getLabel() const {
		return m_nodeLabel;
	}
	void addSubSolved(double d) {
		m_subSolved *= d;
	}
	double getSubSolved() const {
		return m_subSolved;
	}
	void setSubSolved(double d) {
		m_subSolved = d;
	}
//	int getDepth() {
//		//return m_parent->getDepth();
//		return m_depth;
//	}

	double getPathAssignment(std::vector<int>&);

	/* empty implementations, functions meaningless for AND nodes */
	void setCacheContext(const context_t& c) {
		assert(false);
	}
	const context_t& getCacheContext() const {
		assert(false);
		return emptyCtxt;
	}

	/* empty implementations, functions meaningless for AND nodes */
	void setHeurCache(double* d) {
	}
	double* getHeurCache() const {
		return NULL;
	}
	void clearHeurCache() {
	}

	string toString();
public:
	AobbSearchNodeAND(AobbSearchNode* p, val_t val, double label = ELEM_ONE);
	virtual ~AobbSearchNodeAND() { /* empty */
	}
};

class AobbSearchNodeOR: public AobbSearchNode {
protected:
	int m_var;             // Node variable
	int m_depth;           // Depth of corresponding variable in pseudo tree
	double* m_heurCache;   // Stores the precomputed heuristic values of the AND children
	context_t m_cacheContext; // Stores the context (for caching)

public:
	int getType() const {
		return NODE_OR;
	}
	int getVar() const {
		return m_var;
	}
	val_t getVal() const {
		assert(false);
		return NONE;
	} // no val for OR nodes!

	void setValue(double d) {
		m_nodeValue = d;
	}
	double getValue() const {
		return m_nodeValue;
	}
	double getLabel() const {
		assert(false);
		return 0;
	} // no label for OR nodes!
	void addSubSolved(double d) {
		assert(false);
	} // not applicable for OR nodes
	double getSubSolved() const {
		assert(false);
		return 0;
	} // not applicable for OR nodes
	void setSubSolved(double d) {
		assert(false);
	}

	double getPathAssignment(std::vector<int>&);

	void setCacheContext(const context_t& t) {
		m_cacheContext = t;
	}
	const context_t& getCacheContext() const {
		return m_cacheContext;
	}

	void setHeurCache(double* d) {
		m_heurCache = d;
	}
	double* getHeurCache() const {
		return m_heurCache;
	}
	void clearHeurCache();

	string toString();
public:
	AobbSearchNodeOR(AobbSearchNode* parent, int var, int depth);
	virtual ~AobbSearchNodeOR();
};

/* Inline definitions */

inline AobbSearchNode::AobbSearchNode(AobbSearchNode* parent) :
		m_flags(0), m_depth(0), m_parent(parent), m_nodeValue(ELEM_NAN),
		m_heurValue(ELEM_NAN), m_heurUpdated(ELEM_NAN),
		m_children(NULL), m_childCountFull(0), m_childCountAct(0)  {
	/* intentionally empty */
}

inline AobbSearchNode::~AobbSearchNode() {
	this->clearChildren();
}

inline void AobbSearchNode::setChild(AobbSearchNode* node) {
	m_children = new NodeP[1];
	m_children[0] = node;
	m_childCountFull = m_childCountAct = 1;
}

inline void AobbSearchNode::addChildren(
		const std::vector<AobbSearchNode*>& nodes) {
	m_children = new NodeP[nodes.size()];
	for (size_t i = 0; i < nodes.size(); ++i) {
		m_children[i] = nodes[i];
	}
	m_childCountFull = m_childCountAct = nodes.size();
}

inline bool AobbSearchNode::hasChild(AobbSearchNode* node) const {
	for (size_t i = 0; i < m_childCountFull; ++i) {
		if (m_children[i] == node)
			return true;
	}
	return false;
}

inline void AobbSearchNode::eraseChild(AobbSearchNode* node) {
	for (size_t i = 0; i < m_childCountFull; ++i) {
		if (m_children[i] == node) {
			delete m_children[i];
			m_children[i] = NULL;
			--m_childCountAct;
			return;
		}
	}
}

inline void AobbSearchNode::clearChildren() {
	for (size_t i = 0; i < m_childCountFull; ++i) {
		if (m_children[i]) {
			delete m_children[i];
			m_children[i] = NULL;
		}
	}
	m_childCountAct = 0;
	delete[] m_children;
	m_children = NULL;
}

inline void AobbSearchNodeOR::clearHeurCache() {
	if (m_heurCache) {
		delete[] m_heurCache;
		m_heurCache = NULL;
	}
}

inline AobbSearchNodeAND::AobbSearchNodeAND(AobbSearchNode* parent,
		val_t val, double label) :
		AobbSearchNode(parent), m_val(val), m_nodeLabel(label),
		m_subSolved(ELEM_ONE) {
}

inline AobbSearchNodeOR::AobbSearchNodeOR(AobbSearchNode* parent, int var,
		int depth) :
		AobbSearchNode(parent), m_var(var), m_depth(depth), m_heurCache(NULL) //, m_cacheContext(NULL)
{ /* empty */
}

inline AobbSearchNodeOR::~AobbSearchNodeOR() {
	this->clearHeurCache();
}

inline bool AobbSearchNode::heurLess(const AobbSearchNode* a,
		const AobbSearchNode* b) {
	assert(a && b);
	return a->getHeur() < b->getHeur();
}

inline bool AobbSearchNode::heurGreater(const AobbSearchNode* a,
		const AobbSearchNode* b) {
	assert(a && b);
	return a->getHeur() > b->getHeur();
}

inline double AobbSearchNodeAND::getPathAssignment(std::vector<int>& assignment) {

	std::vector<std::pair<int, int> > path;
	AobbSearchNode* c = this;
	AobbSearchNode* p = c->getParent();
	double cost = 1.0;

	while (true) {

		int var = p->getVar(); // OR parent
		int val = c->getVal(); // AND node
		cost *= c->getLabel();
		path.push_back(std::make_pair(var, val));

		c = p->getParent(); // prev AND node
		if (c == NULL) break;
		p = c->getParent(); // prev OR node
	}

	// assumes the assignment vector is already allocated
	for (std::vector<std::pair<int, int> >::iterator pi = path.begin();
			pi != path.end(); ++pi) {
		assignment[pi->first] = pi->second;
	}

	return cost;
}

inline double AobbSearchNodeOR::getPathAssignment(std::vector<int>& assignment) {

	// check if root node
	if (this->getParent() == NULL) {
		return 1.0;
	}

	std::vector<std::pair<int, int> > path;
	AobbSearchNode* c = this->getParent();
	AobbSearchNode* p = c->getParent();
	double cost = 1.0;

	while (true) {

		int var = p->getVar(); // OR parent
		int val = c->getVal(); // AND node
		cost *= c->getLabel();
		path.push_back(std::make_pair(var, val));

		c = p->getParent();
		if (c == NULL) break;
		p = c->getParent();
	}

	// assumes that the assignment vector is already allocated
	for (std::vector<std::pair<int, int> >::iterator pi = path.begin();
			pi != path.end(); ++pi) {
		assignment[pi->first] = pi->second;
	}

	return cost;
}


inline string AobbSearchNodeAND::toString() {

	std::ostringstream oss;

	oss << "AND node: (x" << getVar() << "," << getVal() << ")"
		<< ", l = " << getLabel()
		<< ", h = " << getHeur()
		<< ", q = " << getHeurUpdated()
		<< ", v = " << getValue()
		<< ", d = " << getDepth()
		;

	return oss.str();

}

inline string AobbSearchNodeOR::toString() {

	std::ostringstream oss;
	oss << "OR node: (x" << getVar() << ")"
		<< ", h = " << getHeur()
		<< ", q = " << getHeurUpdated()
		<< ", v = " << getValue()
		<< ", d = " << getDepth()
		;

	return oss.str();
}
 #endif /* SEARCHNODE_H_ */
