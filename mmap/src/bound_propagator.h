/*
 * bound_propagator.h
 *
 *  Created on: Jun 7, 2013
 *      Author: radu
 */

#ifndef IBM_MAP_BOUND_PROPAGATOR_H_
#define IBM_MAP_BOUND_PROPAGATOR_H_

#include "aobb_search_space.h"
#include "aobb_search_node.h"
#include "pseudo_tree.h"
#include "program_options.h"
#include "heuristic.h"

#include "zchaff/zchaff_solver.h"

class Timer;

class BoundPropagator {

protected:

	bool m_doCaching;
	Problem* m_problem;
	AobbSearchSpace* m_space;
	Heuristic* m_heuristic;
	Pseudotree* m_pseudotree;
	ProgramOptions* m_options;
	zchaff::CSolver* m_zchaff;

	double m_lowerBound; // best solution found so far (MAX-PROD)
//	double m_upperBound; // best solution found so far (MIN-SUM), log space

#ifndef NO_ASSIGNMENT
	std::vector<std::vector<val_t> > m_solutions;
#endif

public:

	/*
	 * propagates the value of the specified search node and removes unneeded nodes
	 * returns a pointer to the parent of the highest deleted node
	 * @n: the search node to be propagated
	 * @reportSolution: should root updates be reported to problem instance?
	 */
	AobbSearchNode* propagate(AobbSearchNode* n, Timer& tm,
			std::vector<std::vector<bool> >& currentDomains,
			bool reportSolution = false, AobbSearchNode* upperLimit = NULL);

	AobbSearchNode* propagate(AobbSearchNode* n, Timer& tm,
			bool reportSolution = false, AobbSearchNode* upperLimit = NULL);

//	AobbSearchNode* propagateLog(AobbSearchNode* n, Timer& tm,
//			bool reportSolution, AobbSearchNode* upperLimit = NULL);
//
//	AobbSearchNode* updateLog(AobbSearchNode* n);

#ifndef NO_ASSIGNMENT
	void propagateTuple(AobbSearchNode* start, AobbSearchNode* end);
	std::vector<std::vector<val_t> > getSolutions() {
		return m_solutions;
	}
#endif

	inline double getLowerBound() {
		return m_lowerBound;
	}

//	inline double getUpperBound() {
//		return m_upperBound;
//	}

	inline void setSatSolver(zchaff::CSolver* zc) {
		m_zchaff = zc;
	}

	inline void setHeuristic(Heuristic* h) {
		m_heuristic = h;
	}

private:

#ifndef NO_ASSIGNMENT
	void updateSolution(double timestamp, double cost,
			std::vector<val_t>& sol,
			std::pair<size_t, size_t> nodesMap,
			std::pair<size_t, size_t> nodesAll);
#else
	void updateSolution(double timestamp, double cost,
			std::pair<size_t, size_t> nodesMap,
			std::pair<size_t, size_t> nodesAll);
#endif

//	void updateSolutionLog(double timestamp, double cost, double lowerBound,
//			std::pair<size_t, size_t> nodesMap);
public:
	BoundPropagator(Problem* p, AobbSearchSpace* s, Pseudotree* pt,
			ProgramOptions* opt, bool doCaching = true, Heuristic* h = NULL) :
			m_doCaching(doCaching), m_problem(p), m_space(s), m_heuristic(h),
			m_pseudotree(pt), m_options(opt),
			m_zchaff(NULL), m_lowerBound(ELEM_NAN) { /* empty */
	}
	virtual ~BoundPropagator() {
	}
};

#endif /* BOUND_PROPAGATOR_H_ */
