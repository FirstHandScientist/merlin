/*
 * search.h
 *
 *  Created on: Mar 28, 2013
 *      Author: radu
 */

#ifndef IBM_ANYTIME_SEARCH_H_
#define IBM_ANYTIME_SEARCH_H_

/* A generic search based solver */

#include "solver.h"
#include "pseudo_tree.h"
#include "heuristic.h"
#include "timer.h"

#include "minisat/Solver.h"
#include "zchaff/zchaff_solver.h"

/* Anstract search solver */
class Search : public Solver {

protected:

	// pseudo tree that guids the search
	scoped_ptr<Pseudotree> m_pseudotree;

	// heuristic that guids the search
	scoped_ptr<Heuristic> m_heuristic;

	// search timer
	Timer m_timer;

	// best solution found so far (lower bound)
	double m_lowerBound;

	// upper bound on optimal MAP marginal solution
	double m_upperBound;

	// for constraint propagation
	minisat::Solver m_minisat;
	zchaff::CSolver m_zchaff;
	std::vector<std::vector<int> > m_var2sat;
	std::vector<std::pair<int, int> > m_sat2var;
	std::vector<std::vector<bool> > m_currentDomains; // current domains
	std::vector<int> m_elimOrder;

public:

	// initialize
	virtual void init();

	// encode determinism
	void encodeDeterminism();

	// propagate a partial assignment (full SAT solving)
	bool propagate(const std::vector<int>&);

	// propagate a single variable-value assignment
	bool lookahead(int var, int val, list<pair<int, int> >& changes,
			const std::set<int>& subtree);

public:

	Search(ProgramOptions *o);
	virtual ~Search();
};

inline Search::Search(ProgramOptions *o) : Solver(o) {
	m_lowerBound = ELEM_NAN;
	m_upperBound = ELEM_NAN;
}

inline Search::~Search() {

}

#endif /* SEARCH_H_ */
