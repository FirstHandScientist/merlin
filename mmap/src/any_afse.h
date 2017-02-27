/*
 * any_afse.h
 *
 *  Created on: 14 Aug 2016
 *      Author: radu
 */

#ifndef SRC_ANY_AFSE_H_
#define SRC_ANY_AFSE_H_


#include "solver.h"
#include "timer.h"
#include "mex/afse.h"
#include "pseudo_tree.h"

class AnyAFSE : public Solver {

protected:
	// solver timer
	Timer m_timer;
	std::vector<int> m_elim; // elimination order
	size_t m_iterations;
	mex::afse m_afse;

	// pseudo tree that guids the search
	scoped_ptr<Pseudotree> m_pseudotree;

	// Helper functions to move between DAO & mex factor formats
	mex::vector<mex::Factor> copyFactors(void);

public:
	void init();
	int solve();

public:

	AnyAFSE(ProgramOptions* po) : Solver(po), m_iterations(0) {};
	~AnyAFSE() {};
};


#endif /* SRC_ANY_AFSE_H_ */
