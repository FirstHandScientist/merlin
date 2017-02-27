/*
 * any_arastar.cpp
 *
 *  Created on: Oct 29, 2014
 *      Author: radu
 */


#include "any_arastar.h"

// repair the OPEN list
void AnyARAstar::repair() {

	double upbo = m_upperBound;
	std::vector<AstarSearchNode>& open = m_open.elements();
	std::vector<AstarSearchNode>::iterator it = open.begin();
	for (; it != open.end(); ) {

		if ( (it->g() + it->h()) >= upbo ) {
			it = open.erase(it);
		} else {
			it->f() = it->g() + m_epsilon*it->h();
			++it;
		}
	}

	m_open.heapify();
}

void AnyARAstar::astar() {

	// global constant
	double gconst = -ELEM_ENCODE(m_problem->getGlobalConstant());

	// initialize the weight
	m_epsilon = m_schedule->getInitialWeight();

	// initialize the search
	AstarSearchNode root(UNKNOWN, UNKNOWN);
	root.f() = 0;
	m_open.push(root);

	while (m_epsilon >= 1) {

		while ( !m_open.empty() ) {

			// check for timeout
			if (m_timer.timeout()) {
				throw SEARCH_TIMEOUT;
			}

			// output intermediate results if verbose mode is enabled
			if (m_options->verbose) {
				// logging
#ifdef SEARCH_CHECKPOINT_NODE
				if (m_expansions % SEARCH_CHECKPOINT_NODE == 0 &&
					m_expansions > m_prevNodeCheckpoint) {

					double timestamp = m_timer.elapsed();
					m_prevNodeCheckpoint = m_expansions;
					std::cout << "[*"
						<< std::setw(8) << timestamp << "] w "
						<< std::setw(8) << m_epsilon << " "
						<< std::setw(12) << m_expansions << " "
						<< std::setw(12) << m_upperBound << " (" << 1/ELEM_DECODE(m_upperBound) << ")"
						<< std::endl;
				}
#else
				double elapsed = m_timer.elapsed();
				if (elapsed - m_prevTimeCheckpoint >= SEARCH_CHECKPOINT_TIME) {

					double timestamp = elapsed;
					m_prevTimeCheckpoint = elapsed;
					std::cout << "[*"
						<< std::setw(8) << timestamp << "] w "
						<< std::setw(8) << m_epsilon << " "
						<< std::setw(12) << m_expansions << " "
						<< std::setw(12) << m_upperBound << " (" << 1/ELEM_DECODE(m_upperBound) << ")"
						<< std::endl;
				}
#endif
			}

			AstarSearchNode n = m_open.top();
			m_open.pop();
			if (n.g() + n.h() >= m_upperBound) {
				continue; // this should clear the OPEN list
			}

			if (n.terminal()) {

				// recover the full MAP assignment
				std::vector<int> assignment;
				getMapAssignment(n, assignment);

				// get the cost of the solution path
				double costMap = n.g(); //getMapAssignmentValue(assignment);
				//costMap = -ELEM_ENCODE(costMap);

				// solve the SUM subproblem
				double costSum = _ao(assignment);
				costSum = -ELEM_ENCODE(costSum);
				++m_sumEvals;

				// cost of the MAP solution:
				//double cost = costMap * costSum;
				double cost = costMap + costSum;
				if ( ISNAN(m_upperBound) ) {
					m_upperBound = cost;
				} else {
					m_upperBound = std::min(m_upperBound, cost);
				}

				// solution found
				double timestamp = m_timer.elapsed();
				std::cout << "["
					<< std::setw(9) << timestamp << "] w "
					<< std::setw(8) << m_epsilon << " "
					<< std::setw(12) << m_expansions << " "
					<< std::setw(12) << m_upperBound+gconst << " (" << 1/ELEM_DECODE(m_upperBound+gconst) << ")"
					<< std::endl;
				std::cout.flush();

				if (m_epsilon != 1.0) {
					break; // suboptimal solution found
				}

			} else {
				expand(n);
				++m_expansions;
			}
		} // while (!open.empty())

		// update the weight and repair the OPEN list
		if (m_epsilon == 1.0) { // proved optimality
			m_solved = true;
			m_solutionCost = m_upperBound;
			break;
		}

		// update the weight
		m_epsilon = m_schedule->getNextWeight();

		// repair the search space
		repair();

	} // end while (epsilon)

}

