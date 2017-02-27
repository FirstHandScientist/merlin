/*
 * parallel_rbfaoo_search_node.h
 *
 */

#ifndef PARALLEL_RBFAOO_SEARCH_NODE_H_
#define PARALLEL_RBFAOO_SEARCH_NODE_H_

#include "rbfaoo_search_node.h"

class ParallelRbfaooSearchNode : public RbfaooSearchNode {
protected:
	double m_upperBound;       // upperBound
	double m_pseudoValue;  // heuristic estimate of the node's value

public:

	// constructor
	ParallelRbfaooSearchNode(int node_type, int var, val_t val, int depth, double label = ELEM_ZERO) : 
		RbfaooSearchNode(node_type, var, val, depth, label), m_upperBound(ELEM_NAN), m_pseudoValue(ELEM_NAN) { }

	ParallelRbfaooSearchNode(int node_type, int var, int depth) : 
		RbfaooSearchNode(node_type, var, depth), m_upperBound(ELEM_NAN),m_pseudoValue(ELEM_NAN) { }

	~ParallelRbfaooSearchNode() { }

	void setPseudoValue(double d) { m_pseudoValue = d; }

	double getPseudoValue() const { return m_pseudoValue; }

	void setUpperbound(double d) { m_upperBound = d; }

	double getUpperbound() const { return m_upperBound; }

};

#endif /* PARALLEL_RBFAOO_SEARCH_NODE_H_ */
