/*
 * any_wrbfaoo.cpp
 *
 *  Created on: Oct 29, 2014
 *      Author: radu
 */



#include "any_wrbfaoo.h"

// Anytime Weighted RBFAOO with Restarts
int AnyWRBFAOO::solve() {

	// check if solved during initialization
	if (m_solved) {
		std::cout << "--------- Solved during initialization ---------" << std::endl;
		std::cout << "Problem name:        " << m_problem->getName() << std::endl;
		std::cout << "Status:              " << solver_status[0] << std::endl;
		std::cout << "OR nodes:            " << 0 << std::endl;
		std::cout << "AND nodes:           " << 0 << std::endl;
		std::cout << "OR nodes (MAP):      " << 0 << std::endl;
		std::cout << "AND nodes (MAP):     " << 0 << std::endl;
		std::cout << "Time elapsed:        " << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
		std::cout << "Preprocessing:       " << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
		std::cout << "Solution:            " << m_solutionCost << " (" << ELEM_ENCODE(m_solutionCost) << ")" << std::endl;
		std::cout << "-------------------------------" << std::endl;
	}

	Zobrist::initOnce(*(m_problem.get()));

	// init search space
	Timer tm;
	tm.start();
	m_space.reset(new RbfaooSearchSpace(m_pseudotree.get(), m_options));
	assert( m_options->cacheSize > 0 ); // in kilo bytes
	size_t cache_kilobytes = (m_options->cacheSize);
	if (!m_space->dfpncache) {
		m_space->dfpncache = new RbfaooCacheTable;
		m_space->dfpncache->init(cache_kilobytes);
	}
	tm.stop();
	std::cout << "Cache allocation complete: " << tm.elapsed() << " seconds" << std::endl;

	// prologue
	std::cout << "--- Starting search ---" << std::endl;
	int res = SEARCH_SUCCESS;
	m_tmSolve = 0;
	m_timer.reset(m_options->timeLimit);
	m_timer.start();
	size_t totalAndNodes = 0, totalOrNodes = 0;

	// initialize weight
	m_epsilon = m_schedule->getInitialWeight();
	while (m_epsilon >= 1.0) {

		res = rbfs();

		if (res == SEARCH_SUCCESS) { // new solution found
			if ( ISNAN(m_upperBound) ) {
				m_upperBound = m_solutionCost;
			} else {
				m_upperBound = std::min(m_upperBound, m_solutionCost);
			}
			double timestamp = m_timer.elapsed();
			totalAndNodes += m_space->getAndNodes();
			totalOrNodes += m_space->getOrNodes();

			// output solution
			std::cout << "["
				<< std::setw(9) << timestamp << "] w "
				<< std::setw(8) << m_epsilon << " "
				<< std::setw(12) << m_num_expanded_or_map << " "
				<< std::setw(12) << m_num_expanded_and_map << " "
				<< std::setw(12) << m_upperBound << " (" << 1/ELEM_DECODE(m_upperBound) << ")"
				<< std::endl;

			if (m_epsilon == 1.0) { // proved optimality
				m_solved = true;
				break;
			}
		} else { // out of time or memory
			break;
		}

		// update the weight and continue
		m_epsilon = m_schedule->getNextWeight();

		// reset the cache table
		m_space->dfpncache->clear();
	}

	// stop timer
	m_timer.stop();
	m_tmSolve = m_timer.elapsed();

	// output solution (if found)
	std::cout << std::endl << std::endl;
	std::cout << "--- Search done ---" << std::endl;
	std::cout << "Problem name:        " << m_problem->getName() << std::endl;
	std::cout << "Status:              " << solver_status[res] << std::endl;
	std::cout << "OR nodes:            " << m_num_expanded_or << std::endl;
	std::cout << "AND nodes:           " << m_num_expanded_and << std::endl;
	std::cout << "OR nodes (MAP):      " << m_num_expanded_or_map << std::endl;
	std::cout << "AND nodes (MAP):     " << m_num_expanded_and_map << std::endl;
	std::cout << "Time elapsed:        " << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
	std::cout << "Preprocessing:       " << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
//	std::cout << "Solution:            " << m_solutionCost << " (" << ELEM_ENCODE(m_solutionCost) << ")" << std::endl;
	std::cout << "Solution:            " << 1/ELEM_DECODE(m_solutionCost) << " (" << -m_solutionCost  << ")" << std::endl;
	std::cout << "-------------------------------" << std::endl;

	// clean up
	Zobrist::finishOnce();

	return res;
}


