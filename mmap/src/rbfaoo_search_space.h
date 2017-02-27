/*
 * DfpnpLusSearchSpace.h
 * Inherited from SearchSpace.h 
 */

#ifndef DFPNPLUSSEARCHSPACE_H_
#define DFPNPLUSSEARCHSPACE_H_

#include "rbfaoo_cache_table.h"

struct ProgramOptions;

/* forward declarations */
class Pseudotree;

/* main search space structure for worker nodes. */
class RbfaooSearchSpace {

protected:
	Pseudotree* m_pseudotree;
	ProgramOptions* m_options;

	size_t m_andNodes;
	size_t m_orNodes;

public:
	RbfaooCacheTable* dfpncache; // Cache for dfpn+
	RbfaooSearchNode* root;

public:
	// returns the number of AND nodes expanded
	inline size_t getAndNodes() const {
		return m_andNodes;
	}

	// returns the number of OR nodes expanded
	inline size_t getOrNodes() const {
		return m_orNodes;
	}

	// increment the number of nodes expanded
	inline void incNodesExpanded(int nodeType) {
		assert(nodeType == NODE_AND || nodeType == NODE_OR);
		if (nodeType == NODE_AND) {
			m_andNodes += 1;
		} else {
			m_orNodes += 1;
		}
	}

public:

	RbfaooSearchSpace(Pseudotree* pt, ProgramOptions* opt) :
		m_pseudotree(pt), m_options(opt) {
		dfpncache = NULL;
		root = NULL;
		m_andNodes = 0;
		m_orNodes = 0;
	}

	~RbfaooSearchSpace() {
		if (dfpncache) {
			delete dfpncache;
		}
	}


};

#endif /* DFPNSEARCHSPACE_H_ */
