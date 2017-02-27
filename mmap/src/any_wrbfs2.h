/*
 * any_wrbfs2.h
 *
 *  Created on: 18 Nov 2015
 *      Author: radu
 */

#ifndef SRC_ANY_WRBFS2_H_
#define SRC_ANY_WRBFS2_H_

#include "search.h"
#include "rbfaoo_search_node.h"
#include "rbfaoo_search_space.h"

/**
 * Anytime weighted RBFS with caching. It assumes an OR search graph over the
 * MAP variables (context based caching). The summation subproblems are solved
 * by depth-first AND/OR search with caching (up to the --cache-size limit).
 */

class AnyWRBFS2 : public Search {

protected:

	// The context-minimal AND/OR graph (incl. cache table)
	scoped_ptr<RbfaooSearchSpace> m_space;

	// Reusable vector for node expansions (to avoid repeated (de)allocation of memory)
	std::vector<RbfaooSearchNode*> m_expand;

	// Overestimation
	double m_overestimation;

	// Node expansions internal counters
	size_t m_num_expanded;
	size_t m_num_expanded_or;
	size_t m_num_expanded_and;
	size_t m_num_expanded_or_map;
	size_t m_num_expanded_and_map;

	// Checkpoints
	double m_prevTimeCheckpoint;
	double m_prevNodeCheckpoint;

	// Weight (default 1.0)
	double m_epsilon;

	size_t m_numSumEvals;

	// Last MAP variable
	int m_lastMapVar;

protected:


  /* processes the current node (value instantiation etc.)
	 * if trackHeur==true, caches lower/upper bounds for first node at each depth */
	bool doProcess(RbfaooSearchNode*);

//	void multipleIterativeDeepeningOR(RbfaooSearchNode &n, size_t table_index,
//			double th_value);
	void multipleIterativeDeepeningOR(RbfaooSearchNode &n,
			size_t table_index, double th_value, int dn_threshold, double g_value);

//	void multipleIterativeDeepeningAND(RbfaooSearchNode &n, size_t table_index,
//			double th_value);
	void multipleIterativeDeepeningAND(RbfaooSearchNode &n,
			size_t table_index, double th_value, int dn_threshold, double g_value);

//	int generateChildrenAND(RbfaooSearchNode &n, int &best_index, int &solved_flag);
	int generateChildrenAND(RbfaooSearchNode &n, int &second_best_dn,
			int &best_index, int &solved_flag);

//	void calculateNodeValueAND(RbfaooSearchNode &n, size_t table_index,
//			int num_children, int &best_index, int &solved_flag);
	void calculateNodeValueandDNAND(RbfaooSearchNode &n, size_t table_index,
			int num_children, int &second_best_dn, int &best_index,
			int &solved_flag);

	int generateChildrenOR(RbfaooSearchNode &n, double &second_best_value,
			int &best_index, int &solved_flag, double &cutoff_value);

//	void calculateNodeValueOR(RbfaooSearchNode &n, size_t table_index,
//			int num_children, double &second_best_value, int &best_index,
//			int &solved_flag, double &cutoff_value);
	void calculateNodeValueandDNOR(RbfaooSearchNode &n, size_t table_index,
			int num_children, double &second_best_value, int &best_index,
			int &solved_flag, double &cutoff_value);

	/* computes the heuristic of a new OR node, which requires precomputing
	 * its child AND nodes' heuristic and label values, which are cached
	 * for their explicit generation */
	double heuristicOR(RbfaooSearchNode &n);

	// Recursive best-first search
	int rbfs();

	int sumDN(int dn1, int dn2);

public:

	// initialize
	void init();

	// solver
	virtual int solve();

public:

	AnyWRBFS2(ProgramOptions* o);
	virtual ~AnyWRBFS2() {};

};


/* inline functions */

inline int AnyWRBFS2::sumDN(int dn1, int dn2) {
	if (dn1 == DNINF || dn2 == DNINF)
		return DNINF;
	int d = dn1 + dn2;
	return min(static_cast<int>(DNLARGE), d);
}



#endif /* SRC_ANY_WRBFS2_H_ */
