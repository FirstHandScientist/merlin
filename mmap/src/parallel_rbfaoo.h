/*
 * parallel_rbfaoo.h
 *
 */  


#ifndef PARALLEL_RBFAOO_H_
#define PARALLEL_RBFAOO_H_

#include "search.h"
#include "heuristic.h"
#include "parallel_rbfaoo_cache_table.h"
#include "problem.h"
#include "program_options.h"
#include "pseudo_tree.h"
#include "utils.h"

/* All search algorithms should inherit from this */
class ParallelRBFAOO : public Search {

protected:

	// The context-minimal AND/OR graph (incl. cache table)
	ParallelRbfaooCacheTable *m_cache;

	// Number of threads
	size_t m_num_threads;

	// Checkpoints
	double m_prevTimeCheckpoint;
	double m_prevNodeCheckpoint;

	// Weight
	double m_weight;

public:

	// initialize
	void init();

	// Main solver routine
	virtual int solve();

public:
  
	// Constructor
	ParallelRBFAOO(ProgramOptions* opt);

	// Destructor
	virtual ~ParallelRBFAOO(); 
};

#endif /* PARALLEL_RBFAOO_H_ */
