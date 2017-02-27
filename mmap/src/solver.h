// solver.h -- Abstract solver class.

/*
 * Author: Radu Marinescu
 *
 * Copyright (c) IBM Research, 2012
 *
 * This program is under IBM license. Do not distribute.
 *
 */

/*
 * NOTE: This is an internal header file.
 * You should not attempt to use it directly.
 */

#ifndef IBM_ANYTIME_SOLVER_H_
#define IBM_ANYTIME_SOLVER_H_

#include "base.h"
#include "problem.h"
#include "program_options.h"

/* Abstract solver */
class Solver {

protected:

	// problem instance
	scoped_ptr<Problem> m_problem;
	//Problem *m_problem;

	// program options
	ProgramOptions *m_options;

	// assignment
	std::vector<val_t> m_assignment;

	// load time
	double m_tmLoad;

	// heuristic time
	double m_tmHeuristic;

	// solving time
	double m_tmSolve;

	// flag indicating the problem is solved
	bool m_solved;

	// Optimal solution cost/assignment
	double m_solutionCost;
	std::vector<val_t> m_solution;

public:

	// solve the problem instance
	virtual int solve() = 0;

	// initialize the solver
	virtual void init() = 0;

	// get the assignment (const reference)
	const std::vector<int>& assignment() const;

	// get the assignment (reference)
	std::vector<int>& assignment();

	// returns the preprocessing time (will be set in subclasses)
	inline double getLoadTime() {
		return m_tmLoad;
	}

	// returns the heuristic compilation time (will be set in subclasses)
	inline double getHeursticTime() {
		return m_tmHeuristic;
	}

	// returns the solving time (will be set in subclasses)
	inline double getSolveTime() {
		return m_tmSolve;
	}

public:

	// default constructor
	Solver(ProgramOptions* o = NULL);

	// destructor
	virtual ~Solver();
};

/* inline declarations */


// default constructor
inline Solver::Solver(ProgramOptions *o) :
		m_options(o),
		m_tmLoad(0),
		m_tmHeuristic(0),
		m_tmSolve(0),
		m_solved(false),
		m_solutionCost(ELEM_NAN) {
	/* empty */
}

// destructor
inline Solver::~Solver() {

}

// get the assignment (const reference)
inline const std::vector<int>& Solver::assignment() const {
	return m_assignment;
}

// get the assignment (reference)
inline std::vector<int>& Solver::assignment() {
	return m_assignment;
}

#endif /* SOLVER_H_ */
