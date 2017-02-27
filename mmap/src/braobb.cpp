/*
 * braobb.cpp
 *
 *  Created on: Jun 7, 2013
 *      Author: radu
 */


#include "braobb.h"
#include "bound_propagator.h"

AobbSearchNode* BRAOBB::nextNode() {
	AobbSearchNode* n = NULL;

	// check first if we're in a SUM problem
	if (!n && m_stack.size()) {
		n = m_stack.top();
		m_stack.pop();
		return n;
	}

	while (m_stacks.size()) {
		MyStack* st = m_stacks.front();
		if (st->getChildren()) {
			m_stacks.push(st);
			m_stackCount = 0;
		} else if (st->empty()) {
			if (st->getParent())
				st->getParent()->delChild();
			delete st;
			m_stackCount = 0;
		} else if (m_stackLimit && m_stackCount++ == m_stackLimit) {
			m_stacks.push(st);
			m_stackCount = 0;
		} else {
			n = st->top();
			st->pop();
			break; // while loop
		}
		m_stacks.pop();
	}
	return n;
}

bool BRAOBB::doExpand(AobbSearchNode* n) {
	assert(n);
	m_expand.clear();

	MyStack* stack = m_stacks.front();
	if (n->getType() == NODE_AND) {  // AND node
		if (generateChildrenAND(n, m_expand))
			return true; // no children

		if (m_problem->isSum(m_expand.at(0)->getVar())) {
			for (vector<AobbSearchNode*>::reverse_iterator it = m_expand.rbegin();
					it != m_expand.rend(); ++it) {
				m_stack.push(*it);
			}
		} else { // MAP variables
//			std::cout << "[DEBUG] expanded " << n->toString() << endl;
			if (m_expand.size() == 1) {  // no decomposition
				stack->push(m_expand.at(0));
			} else {  // decomposition, split stacks
				// reverse iterator needed since new stacks are put in queue (and not depth-first stack)
				for (vector<AobbSearchNode*>::reverse_iterator it = m_expand.rbegin();
						it != m_expand.rend(); ++it) {
					MyStack* s = new MyStack(stack);
					stack->addChild();
					s->push(*it);
					m_stacks.push(s);
				}
			}
		}
	} else {  // OR node
		if (generateChildrenOR(n, m_expand))
			return true; // no children

		if (m_problem->isSum(m_expand.at(0)->getVar())) {
			for (vector<AobbSearchNode*>::reverse_iterator it = m_expand.rbegin();
					it != m_expand.rend(); ++it) {
				m_stack.push(*it);
			}  // for loop
		} else { // MAP variables
//			std::cout << "[DEBUG] expanded " << n->toString() << endl;
			for (vector<AobbSearchNode*>::reverse_iterator it = m_expand.rbegin();
					it != m_expand.rend(); ++it) {
				stack->push(*it);
			}  // for loop
		}
	}
	// if over node type

	return false; // default false
}

// solve the problem
int BRAOBB::solve() {

#ifndef UAI_COMPETITION
	// check if solved during initialization
	if (m_solved) {
		std::cout << "--------- Solved during initialization ---------" << std::endl;
		std::cout << "Problem name:\t" << m_problem->getName() << std::endl;
		std::cout << "Status:\t\t" << solver_status[0] << std::endl;
		std::cout << "OR nodes:\t" << 0 << std::endl;
		std::cout << "AND nodes:\t" << 0 << std::endl;
		std::cout << "Cache hits:\t" << 0 << std::endl;
		std::cout << "Time elapsed:\t" << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
		std::cout << "Preprocessing:\t" << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
		std::cout << "Solution:\t" << std::setiosflags(std::ios::fixed) << std::setprecision(4) << m_solutionCost << std::endl;
		std::cout << "-------------------------------" << std::endl;

		return SEARCH_SUCCESS;
	}
#endif

	// prologue
	std::cout << "--- Starting search ---" << std::endl;
	int res = SEARCH_SUCCESS;
	char* EmergencyMemory = new char[16384]; // in order to handle bad alloc

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
		AobbSearchNode* first = initSearchSpace();
		if (first) {
			m_rootStack = new MyStack(NULL);
			m_rootStack->push(first);
			m_stacks.push(m_rootStack);
			m_stackCount = 0;
		}

		//std::cout << "[DEBUG] first node is " << first->toString() << std::endl;

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
		std::cout << "OUT OF MEMORY";
		res = SEARCH_OUT_OF_MEMORY;
	} catch (int& e) {
		std::cout << "TIMEOUT";
		res = (e == SEARCH_TIMEOUT) ? SEARCH_TIMEOUT : SEARCH_FAILURE;
	} catch (...) {
		std::cout << "unexpected error: ";
		res = SEARCH_FAILURE;
	}

	// stop timer
	m_timer.stop();
	m_tmSolve = m_timer.elapsed();
	m_solutionCost = m_propagator->getLowerBound();

	// output solution (if found)
	std::cout << std::endl << std::endl;
	std::cout << "----------------- Search done ------------------" << std::endl;
	std::cout << "Problem name:        " << m_problem->getName() << std::endl;
	std::cout << "Status:              " << solver_status[res] << std::endl;
	std::cout << "OR nodes:            " << m_space->nodesOR << std::endl;
	std::cout << "AND nodes:           " << m_space->nodesAND << std::endl;
	std::cout << "OR nodes (MAP):      " << m_space->nodesORmap << std::endl;
	std::cout << "AND nodes (MAP):     " << m_space->nodesANDmap << std::endl;
	std::cout << "Cache hits:          " << m_cacheHits << std::endl;
	std::cout << "Deadends (CP):       " << m_deadendsCP << std::endl;
	std::cout << "Deadends (UB):       " << m_deadendsUB << std::endl;
	std::cout << "Time elapsed:        " << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
	std::cout << "Preprocessing:       " << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
	std::cout << "Solution:            " << m_solutionCost << " (" << ELEM_ENCODE(m_solutionCost) << ")" << std::endl;
	std::cout << "------------------------------------------------" << std::endl;

	// clean up
	if (EmergencyMemory) delete[] EmergencyMemory;

	return res;
}

