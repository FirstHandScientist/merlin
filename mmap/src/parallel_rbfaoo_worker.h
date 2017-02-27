/*
 * parallel_rbfaoo_worker.h
 * Programmed by Akihiro Kishimoto
 *
 */  


#ifndef PARALLEL_RBFAOO_WORKER_H
#define PARALLEL_RBFAOO_WORKER_H

#include "parallel_rbfaoo_search_node.h"
#include "heuristic.h"
#include "problem.h"
#include "program_options.h"
#include "pseudo_tree.h"
#include "parallel_rbfaoo_cache_table.h"
#include "timer.h"
#include "utils.h"

/* All search algorithms should inherit from this */
class ParallelRBFAOOWorker {

protected:
	// Thread it
	int m_thread_id;

	size_t m_depth;

	unsigned int m_num_threads;

	//Information that must be passeed from ParallelRBFAOO
	Problem *m_problem;
	Pseudotree *m_pseudotree;
	ProgramOptions *m_options;
	ParallelRbfaooCacheTable *m_cache;
	Heuristic *m_heuristic;
	// assignment
	std::vector<val_t> m_assignment;

	Timer m_timer;

	// Reusable vector for node expansions (to avoid repeated (de)allocation of memory)
	std::vector<ParallelRbfaooSearchNode*> m_expand;

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

	// Weight
	double m_weight;

	double m_solutionCost;
	
	int m_result;

public:
	void start();
	
	int getSearchResult() const { return m_result; }

	double getSolutionCost() const { return m_solutionCost; }

	size_t getNumExpandedOR() const { return m_num_expanded_or; }

	size_t getNumExpandedAND() const { return m_num_expanded_and; }

	int getThreadId() const { return m_thread_id; }

protected:

	/* processes the current node (value instantiation etc.) */
	bool doProcess(ParallelRbfaooSearchNode*);

	void multipleIterativeDeepeningOR(ParallelRbfaooSearchNode &n,
					  size_t table_index, double th_value, double upper_bound);


	void multipleIterativeDeepeningAND(ParallelRbfaooSearchNode &n,
					   size_t table_index, double th_value, double upper_bound);

	int generateChildrenAND(ParallelRbfaooSearchNode &n, int &second_best_dn,
				int &best_index, unsigned int &solved_flag);

	void calculateNodeValueandDNAND(ParallelRbfaooSearchNode &n, size_t table_index,
					int num_children, int &second_best_dn, int &best_index,
					unsigned int &solved_flag);

	int generateChildrenOR(ParallelRbfaooSearchNode &n, double &pseudo_second_best_value, 
			       int &pseudo_best_index, unsigned int &solved_flag);


	void calculateNodeValueandDNOR(ParallelRbfaooSearchNode &n, size_t table_index, 
				       int num_children, double &pseudo_second_best_value, 
				       int &pseudo_best_index, unsigned int &solved_flag);



	/* computes the heuristic of a new OR node, which requires precomputing
	 * its child AND nodes' heuristic and label values, which are cached
	 * for their explicit generation */
	double heuristicOR(ParallelRbfaooSearchNode &n);

	int sumDN(int dn1, int dn2);

	// Recursive best-first search
	int rbfs();

public:
  
	// Constructor
	ParallelRBFAOOWorker(int thread_id, Problem *problem, Pseudotree *pseudotree, ProgramOptions* opt, ParallelRbfaooCacheTable *cache, Heuristic *heuristic); 

	// Destructor
	virtual ~ParallelRBFAOOWorker() {}
};

/* Inline definitions */

inline int ParallelRBFAOOWorker::sumDN(int dn1, int dn2) {
	if (dn1 == DNINF || dn2 == DNINF)
		return DNINF;
	int d = dn1 + dn2;
	return min(static_cast<int>(DNLARGE), d);
}

#endif /* PARALLEL_RBFAOO_WORKER_H */
