/*
 * bb_park.cpp
 *
 *  Created on: Sep 26, 2013
 *      Author: radu
 */


#include "bb_park.h"
#include "join_tree_heur.h"

//#define DEBUG

void BBPark::init() {

	// init problem instance
	Search::init();

	// add dummy functions, just to be consistent (used only by the bucket tree)
	PseudotreeNode* r = m_pseudotree->getRoot();
	m_problem->addDummyFuns();
	std::cout << "Adding dummy functions: ";
	for (std::vector<PseudotreeNode*>::const_iterator it = r->getChildren().begin();
			it != r->getChildren().end(); ++it) {
		PseudotreeNode* c = (*it);
		int rVar = r->getVar();
		int cVar = c->getVar();
		int tabSize = m_problem->getDomainSize(cVar);
		double* tab = new double[tabSize];
		for (int i = 0; i < tabSize; ++i) tab[i] = 1.0;
		std::set<int> scope;
		scope.insert(rVar);
		scope.insert(cVar);
		int fid = m_problem->getC();
		Function* fn = new FunctionBayes(fid++, m_problem.get(), scope, tab, tabSize);
		m_problem->addFunction(fn);

		std::cout << *fn << " ";
	}
	std::cout << std::endl;

	// start a timer
	Timer tm;

	// create the mini-bucket heuristic
	assert(m_options->heuristic == HEUR_JOIN_TREE);

	tm.reset();
	tm.start();
	std::cout << "Computing the join-tree heuristic..." << std::endl;
	m_heuristic.reset(new JoinTreeHeur(m_problem.get(),
			NULL, m_options));
	size_t sz = m_heuristic->build(&m_assignment, true);

	tm.stop();
	m_tmHeuristic = tm.elapsed();

	std::cout << "\tJoin tree finished in " << m_tmHeuristic
			<< " seconds" << std::endl;
	std::cout << "\tUsing " << (sz / (1024*1024.0)) * sizeof(double)
			<< " MBytes of RAM" << std::endl;

	// Record time after initialization
	cout << "Initialization complete: " << (m_tmLoad + m_tmHeuristic)
			<< " seconds" << std::endl;

	// check for upper bound computation only
	if (m_options->upperBoundOnly) {
		std::cout << "Done." << std::endl;
		exit(EXIT_SUCCESS);
	}

}

// select next unassigned MAP variable
int BBPark::selectNextVar(const set<int>& M,
		std::vector<std::pair<int, double> >& values,
		const std::vector<int>& assignment) {

	values.clear();
	int nextVar = *(M.begin());
	double ratio = -INFINITY;
	for (std::set<int>::iterator it = M.begin(); it != M.end(); ++it) {
		int var = (*it);
		std::vector<double> bounds;
		m_heuristic->getHeur(var, bounds, assignment);

		std::vector<std::pair<int, double> > temp;
		double Tv = 0, Mv = -INFINITY;
		for (int val = 0; val < m_problem->getDomainSize(var); ++val) {
			temp.push_back(std::make_pair(val, bounds[val]));
			Mv = std::max(Mv, bounds[val]);
			if (bounds[val] >= m_solutionCost) {
				Tv += bounds[val];
			}
		}

		double r = Mv/Tv;
		if (r > ratio) {
			ratio = r;
			nextVar = var;
			values = temp;
		}
	}

	// safety checks
	assert(nextVar != NONE);
	return nextVar;
}

// recursive depth first branch and bound over the MAP variables
void BBPark::bb(std::set<int> M, std::vector<int> assignment) {

	// check for time limit violation
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
				<< std::setw(8) << timestamp << "] u "
				<< std::setw(8) << -1 << " "
				<< std::setw(12) << m_expansions << " "
				<< std::setw(12) << m_solutionCost
				<< std::endl;
		}
#else
		double elapsed = m_timer.elapsed();
		if (elapsed - m_prevTimeCheckpoint >= SEARCH_CHECKPOINT_TIME) {

			double timestamp = elapsed;
			m_prevTimeCheckpoint = elapsed;
			std::cout << "[*"
				<< std::setw(8) << timestamp << "] u "
				<< std::setw(8) << -1 << " "
				<< std::setw(12) << m_expansions << " "
				<< std::setw(12) << m_solutionCost
				<< std::endl;
		}
#endif
	}

	if (M.empty()) { // new solution found (all MAP variables assigned)
		m_solutionCost = m_heuristic->getGlobalUB();
		double timestamp = m_timer.elapsed();
		std::cout << "["
			<< std::setw(9) << timestamp << "] u "
			<< std::setw(8) << -1 << " "
			<< std::setw(12) << m_expansions << " "
			<< std::setw(12) << m_solutionCost
			<< std::endl;
		std::cout.flush();
		return;
	}

	// variable ordering heuristic
	std::vector<std::pair<int, double> > values;
	int nextVar = selectNextVar(M, values, assignment);
	M.erase(nextVar);
	assert(assignment[nextVar] == NONE);

	// sort values in descending order of their heuristics
	std::sort(values.begin(), values.end(), descHeurComp());

	// node expansion
	for (std::vector<std::pair<int, double> >::iterator it = values.begin();
			it != values.end(); ++it) {

		int val = it->first;
		double ub = it->second;
		if (ub > m_solutionCost) { // also prune infeasible values

			++m_expansions;
			assignment[nextVar] = val;
			m_heuristic->update(assignment); // propagate the new assignment
			bb(M, assignment); // recursive call
		}
	}
}

int BBPark::solve() {

	// prologue
	std::cout << "--- Starting search ---" << std::endl;
	int res = SEARCH_SUCCESS;
	char* EmergencyMemory = new char[16384]; // in order to handle bad alloc

	// init search (timer, search space, assignment)
	m_tmSolve = 0;
	m_timer.reset(m_options->timeLimit);
	m_timer.start();
	m_solutionCost = -INFINITY; // holds the current best MAP
	m_assignment.resize(m_problem->getN(), NONE);
	m_expansions = 0;
	m_assignment[m_problem->getDummyVar()] = 0; // assign the dummy variable

	try {

		// select MAP variables
		std::set<int> M;
		for (int var = 0; var < m_problem->getN(); ++var) {
			if (m_problem->isMap(var) && var != m_problem->getDummyVar()) {
			//if (m_problem->isMap(var)) {
				M.insert(var);
			}
		}

		// depth first recursive branch and bound over the MAP variables
		bb(M, m_assignment);

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
	m_solutionCost *= m_problem->getGlobalConstant();

	// output solution (if found)
	std::cout << std::endl << std::endl;
	std::cout << "----------------- Search done ------------------" << std::endl;
	std::cout << "Problem name:\t" << m_problem->getName() << std::endl;
	std::cout << "Status:\t\t" << solver_status[res] << std::endl;
	std::cout << "Nodes:\t\t" << m_expansions << std::endl;
	std::cout << "Time elapsed:\t" << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
	std::cout << "Preprocessing:\t" << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
	std::cout << "Solution:\t" << m_solutionCost << " (" << ELEM_ENCODE(m_solutionCost) << ")" << std::endl;
	std::cout << "------------------------------------------------" << std::endl;

	// clean up
	if (EmergencyMemory) delete[] EmergencyMemory;

	return res;
}
