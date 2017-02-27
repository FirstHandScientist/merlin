/*
 * rbfaoo.h
 *
 *  Created on: Apr 2, 2014
 *      Author: radu
 */

#ifndef IBM_MAP_RBFAOO_H_
#define IBM_MAP_RBFAOO_H_

#include "search.h"

#include "rbfaoo_search_node.h"
#include "rbfaoo_search_space.h"

#include "aobb_search_node.h"
#include "aobb_search_space.h"
#include "bound_propagator.h"

class RBFAOO : public Search {

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

protected:


  /* processes the current node (value instantiation etc.)
	 * if trackHeur==true, caches lower/upper bounds for first node at each depth */
	bool doProcess(RbfaooSearchNode*);

//	void multipleIterativeDeepeningOR(RbfaooSearchNode &n, size_t table_index,
//			double th_value);
	void multipleIterativeDeepeningOR(RbfaooSearchNode &n,
			size_t table_index, double th_value, int dn_threshold);

//	void multipleIterativeDeepeningAND(RbfaooSearchNode &n, size_t table_index,
//			double th_value);
	void multipleIterativeDeepeningAND(RbfaooSearchNode &n,
			size_t table_index, double th_value, int dn_threshold);

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

//private:
//
//	// AO search (for the conditioned SUM problem)
//	std::stack<AobbSearchNode*> ao_stack;
//	std::vector<AobbSearchNode*> ao_expand;
//	scoped_ptr<AobbSearchSpace> ao_space;
//	scoped_ptr<BoundPropagator> ao_propagator;
//
//	// Solve the conditioned SUM subproblem by AO search
//	double _ao(int var, const std::vector<int>& assignment);
//
//	void _initAoSearch(int var);
//	void _heuristicOR(AobbSearchNode*);
//	bool _doProcess(AobbSearchNode*);
//	void _addCacheContext(AobbSearchNode*, const set<int>&) const;
//	bool _doCaching(AobbSearchNode* node);
//	AobbSearchNode* _nextLeaf();
//	bool _generateChildrenAND(AobbSearchNode*, vector<AobbSearchNode*>&);
//	bool _generateChildrenOR(AobbSearchNode* n, vector<AobbSearchNode*>& chi);
//	AobbSearchNode* _nextNode();
//	bool _doExpand(AobbSearchNode* n);
//	bool _doPruning(AobbSearchNode* n);
//	bool _canBePruned(AobbSearchNode* n);

public:

	RBFAOO(ProgramOptions* o);
	virtual ~RBFAOO() {};

};


/* inline functions */

inline int RBFAOO::sumDN(int dn1, int dn2) {
	if (dn1 == DNINF || dn2 == DNINF)
		return DNINF;
	int d = dn1 + dn2;
	return min(static_cast<int>(DNLARGE), d);
}


#endif /* RBFAOO_H_ */
