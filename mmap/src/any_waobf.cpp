/*
 * any_waobf.cpp
 *
 *  Created on: Oct 29, 2014
 *      Author: radu
 */


#include "any_waobf.h"

int AnyWAOBF::solve() {

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

		return SEARCH_SUCCESS;
	}

	// prologue
	std::cout << "--- Starting search ---" << std::endl;
	int res = SEARCH_SUCCESS;
	char* EmergencyMemory = new char[16384]; // in order to handle bad alloc

	// init search (timer, search space, assignment)
	m_tmSolve = 0;
	m_timer.reset(m_options->timeLimit);
	m_timer.start();
	size_t totalAndNodes = 0, totalOrNodes = 0;
	size_t totalAndNodesMAP = 0, totalOrNodesMAP = 0;

	try {

		// initialize weight
		m_epsilon = m_schedule->getInitialWeight();
		while (m_epsilon >= 1.0) {

			// restart the solver
			m_solutionCost = ELEM_NAN;
			m_assignment.resize(m_problem->getN(), UNKNOWN);
			m_space.reset(new AobfSearchSpace());

			// init the AO search space
			ao_space.reset(new AobbSearchSpace());
			ao_propagator.reset(new BoundPropagator(m_problem.get(), ao_space.get(), m_pseudotree.get(), m_options));
			if (m_options->propagation == CP_UNIT) {
				ao_propagator->setSatSolver(&m_zchaff);
			} else {
				ao_propagator->setSatSolver(NULL);
			}
			ao_space->cache = new AobbCacheTable(m_problem->getN());

			// weighted AO* search
			res = aostar(m_options->verbose);
			m_solutionCost += -ELEM_ENCODE(m_problem->getGlobalConstant());

			if (res == SEARCH_SUCCESS) { // new solution found
				if ( ISNAN(m_upperBound) ) {
					m_upperBound = m_solutionCost;
				} else {
					m_upperBound = std::min(m_upperBound, m_solutionCost);
				}

				double timestamp = m_timer.elapsed();
				totalAndNodesMAP += m_space->getAndNodes();
				totalOrNodesMAP += m_space->getOrNodes();
				totalAndNodes += ao_space->nodesAND + m_space->getAndNodes();
				totalOrNodes += ao_space->nodesOR + m_space->getOrNodes();

				// output solution
				std::cout << "["
					<< std::setw(9) << timestamp << "] w "
					<< std::setw(8) << m_epsilon << " "
					<< std::setw(12) << totalOrNodesMAP << " "
					<< std::setw(12) << totalAndNodesMAP << " "
					<< std::setw(12) << m_upperBound << " ("
					<< std::setw(12) << 1/ELEM_DECODE(m_upperBound) << ")"
					<< std::endl;

				if (m_epsilon == 1.0) { // proved optimality
					m_solved = true;
					break;
				}
			} else if (res == SEARCH_TIMEOUT) {
				break;
			}

			// update the weight and continue
			m_epsilon = m_schedule->getNextWeight();
		}

	} catch (std::bad_alloc& ba) {
		delete[] EmergencyMemory;
		EmergencyMemory = NULL;
		std::cerr << "OUT OF MEMORY";
		res = SEARCH_OUT_OF_MEMORY;
	} catch (int& e) {
		std::cout << "TIMEOUT";
		res = (e == SEARCH_TIMEOUT) ? SEARCH_TIMEOUT : SEARCH_FAILURE;
	} catch (...) {
		std::cerr << "unexpected error";
		res = SEARCH_FAILURE;
	}

	// stop timer
	m_timer.stop();
	m_tmSolve = m_timer.elapsed();


	// output solution (if found)
	std::cout << std::endl << std::endl;
	std::cout << "--- Search done ---" << std::endl;
	std::cout << "Problem name:        " << m_problem->getName() << std::endl;
	std::cout << "Status:              " << solver_status[res] << std::endl;
	std::cout << "OR nodes:            " << ao_space->nodesOR + m_space->getOrNodes() << std::endl;
	std::cout << "AND nodes:           " << ao_space->nodesAND + m_space->getAndNodes() << std::endl;
	std::cout << "OR nodes (MAP):      " << m_space->getOrNodes() << std::endl;
	std::cout << "AND nodes (MAP):     " << m_space->getAndNodes() << std::endl;
	std::cout << "Time elapsed:        " << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
	std::cout << "Preprocessing:       " << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
//	std::cout << "Solution:            " << m_solutionCost << " (" << ELEM_ENCODE(m_solutionCost) << ")" << std::endl;
	std::cout << "Solution:            " << 1/ELEM_DECODE(m_solutionCost) << " (" << -m_solutionCost  << ")" << std::endl;
	std::cout << "-------------------------------" << std::endl;

	// clean up
	if (EmergencyMemory) delete[] EmergencyMemory;

	return res;
}


