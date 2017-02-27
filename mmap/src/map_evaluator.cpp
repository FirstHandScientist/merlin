/*
 * evaluator.cpp
 *
 *  Created on: Oct 22, 2014
 *      Author: radu
 */

#include "map_evaluator.h"
#include "ao.h"
#include "mex/wmb.h"

int MapEvaluator::solve() {

	std::cout << "--- Starting evaluation ---" << std::endl;
	int res = SEARCH_SUCCESS;
	char* EmergencyMemory = new char[16384]; // in order to handle bad alloc

	if (m_options->evaluation == EVAL_VE) {

		try {

			m_tmSolve = 0;
			m_timer.reset(m_options->timeLimit);
			m_timer.start();
			m_solutionCost = ELEM_NAN;

			// init weighted mini-buckets
			mex::vector<mex::Factor> fs(m_problem->getC());
			for (int i = 0; i < m_problem->getC(); ++i)
				fs[i] = m_problem->getFunctions()[i]->asFactor();

			mex::wmb ve(fs);
			ve.setIBound(6);
//			ve.init();
			ve.initJG();
			ve.tighten(100);
			m_solutionCost = ve.lb() + ELEM_ENCODE(m_problem->getGlobalConstant());

			res = SEARCH_SUCCESS;

		}  catch (std::bad_alloc& ba) {
			delete[] EmergencyMemory;
			EmergencyMemory = NULL;
			std::cerr << "out of memory";
			res = SEARCH_FAILURE;
		}

		// stop timer
		m_timer.stop();
		m_tmSolve = m_timer.elapsed();

		// output solution (if found)
		std::cout <<  std::endl << std::endl;
		std::cout << "----------------- Evaluation done ------------------" << std::endl;
		std::cout << "Problem name:\t\t" << m_problem->getName() << std::endl;
		std::cout << "Status:\t\t\t" << solver_status[res] << std::endl;
		std::cout << "Time elapsed:\t\t" << (m_tmLoad + m_tmSolve) << " seconds" << std::endl;
		std::cout << "Solution:\t\t" << ELEM_DECODE(m_solutionCost) << " (" << m_solutionCost << ")" << std::endl;
		std::cout << "----------------------------------------------------" << std::endl;

	} else {

		// init search (timer, search space, assignment)
		m_tmSolve = 0;
		m_timer.reset(m_options->timeLimit);
		m_timer.start();
		m_space.reset(new AobbSearchSpace());
		m_solutionCost = ELEM_NAN;
		m_assignment.resize(m_problem->getN(), UNKNOWN);
		m_propagator.reset(new BoundPropagator(m_problem.get(), m_space.get(), m_pseudotree.get(), m_options));
		if (m_options->propagation == CP_UNIT) {
			m_propagator->setSatSolver(&m_zchaff);
		} else {
			m_propagator->setSatSolver(NULL);
		}

		try {

			// init cache table
			if (!m_space->cache) {
				m_space->cache = new AobbCacheTable(m_problem->getN());
			}

			// init root of the search space
			AobbSearchNode* first = initSearch();
			if (first) {
				m_stack.push(first);
			}

			AobbSearchNode* n = nextLeaf();
			while (n) { // throws timeout
				m_propagator->propagate(n, m_timer, m_currentDomains, true); // true = report solutions
				m_lowerBound = m_propagator->getLowerBound();
				n = nextLeaf();
			}

			// Proved optimality
			m_solved = true;
		} catch (std::bad_alloc& ba) {
			delete[] EmergencyMemory;
			EmergencyMemory = NULL;
			std::cerr << "out of memory";
			res = SEARCH_FAILURE;
		} catch (int e) {
			if (e == SEARCH_TIMEOUT) {
				res = SEARCH_TIMEOUT;
			}

		}

		// stop timer
		m_timer.stop();
		m_tmSolve = m_timer.elapsed();
		m_solutionCost = m_propagator->getLowerBound();

		// output solution (if found)
		std::cout <<  std::endl << std::endl;
		std::cout << "----------------- Evaluation done ------------------" << std::endl;
		std::cout << "Problem name:\t\t" << m_problem->getName() << std::endl;
		std::cout << "Status:\t\t\t" << solver_status[res] << std::endl;
		std::cout << "OR nodes:\t\t" << m_space->nodesOR << std::endl;
		std::cout << "AND nodes:\t\t" << m_space->nodesAND << std::endl;
		std::cout << "Time elapsed:\t\t" << (m_tmLoad + m_tmSolve) << " seconds" << std::endl;
		std::cout << "Solution:\t\t" << m_solutionCost << " (" << ELEM_ENCODE(m_solutionCost) << ")" << std::endl;
		std::cout << "----------------------------------------------------" << std::endl;

	}

	// clean up
	if (EmergencyMemory) delete[] EmergencyMemory;

	return res;
}
