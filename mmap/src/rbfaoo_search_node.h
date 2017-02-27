/*
 * DfpnplusSearchNode.h
 *
 */

#ifndef DFPNPLUSSEARCHNODE_H_
#define DFPNPLUSSEARCHNODE_H_

#include "base.h"
#include "utils.h"

#include "hash_zobrist.h"

enum {
	DNINF = 1 << 30, DNLARGE = 1 << 28, UNDETERMINED = -1
};

typedef struct {
	double* m_heurCache; // Stores the precomputed heuristic values of the AND children

} or_node_union_t;

typedef struct {
	double m_nodeLabel; // Label of arc <X_i,a>, i.e. instantiated function costs
	val_t m_val;          // Node value, assignment to OR parent variable
} and_node_union_t;

class RbfaooSearchNode {
protected:
	double m_nodeValue;       // node value (as in cost)
	double m_heurValue;       // heuristic estimate of the node's value
	context_t m_cacheContext; // Stores the context (for caching)
	Zobrist m_zobrist;
	int m_var;				// Node variable
	int m_depth;			// Depth of corresponding variable in pseudo tree
	int m_nodeType;
	int m_dn;
	union {
		or_node_union_t m_orNodeInfo;
		and_node_union_t m_andNodeInfo;
	};

public:

	// constructor
	RbfaooSearchNode(int node_type, int var, val_t val, int depth, double label = ELEM_ZERO) :
			m_nodeValue(ELEM_NAN), m_heurValue(ELEM_NAN), m_var(var),
			m_depth(depth), m_nodeType(node_type), m_dn(0) {
    ; assert(node_type == NODE_AND);  
    ; assert(m_cacheContext.size() == 0);
    m_andNodeInfo.m_val = val;
    m_andNodeInfo.m_nodeLabel = label;
	}

	RbfaooSearchNode(int node_type, int var, int depth) :
		m_nodeValue(ELEM_NAN), m_heurValue(ELEM_NAN), m_var(var),
		m_depth(depth), m_nodeType(node_type), m_dn(0) {
		; assert(node_type == NODE_OR);
		; assert(m_cacheContext.size() == 0);
		m_orNodeInfo.m_heurCache = NULL;
	}

	~RbfaooSearchNode() {
		if (m_nodeType == NODE_OR && m_orNodeInfo.m_heurCache)
		delete[] m_orNodeInfo.m_heurCache;
	}

	inline int getType() const {return m_nodeType;}

	inline int getVar() const {return m_var;}

	inline val_t getVal() const {return m_andNodeInfo.m_val;}

	string toString(const RbfaooSearchNode* a);

	inline void setCacheContext(const context_t& t) {m_cacheContext = t;}

	inline void addCacheContext(val_t v) {m_cacheContext.push_back(v);}

	inline context_t& getCacheContext() {return m_cacheContext;}

	inline void setNodeValue(double d) {m_nodeValue = d;}

	inline double getNodeValue() const {return m_nodeValue;}

	inline int getDN() const {return m_dn;}

	inline void setDN(int dn) {m_dn = dn;}

	inline int getDepth() const {return m_depth;}

	inline double getLabel() const {
		; assert(m_nodeType == NODE_AND);
		return m_andNodeInfo.m_nodeLabel;
	}

	inline void setLabel(double d) {
		; assert(m_nodeType == NODE_AND);
		m_andNodeInfo.m_nodeLabel = d;
	}

	inline void setHeurCache(double* d) {
		; assert(m_nodeType == NODE_OR);
		m_orNodeInfo.m_heurCache = d;
	}

	inline double* getHeurCache() const {
		; assert(m_nodeType == NODE_OR);
		return m_orNodeInfo.m_heurCache;
	}

	inline void setHeur(double d) {m_heurValue = d;}

	inline double getHeur() const {return m_heurValue;}

	inline void clearHeurCache();

	inline Zobrist &getZobrist() {return m_zobrist;}

	friend ostream& operator <<(ostream& out, const RbfaooSearchNode& n) {

		if (n.getType() == NODE_OR) {
			out << "OR (x" << n.getVar() << "): v=" << n.getNodeValue() << ", h=" << n.getHeur();
			double *cache = n.getHeurCache();
			out << " w0=" << cache[2*0+1] << ", w1=" << cache[2*1+1];
		} else {
			out << "AND (x"<< n.getVar() << "=" << n.getVal() << "): v=" << n.getNodeValue() << ", h=" << n.getHeur();
		}
		return out;
	}

};



#endif /* DFPNPLUSSEARCHNODE_H_ */
