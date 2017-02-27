/*
 * any_afse.cpp
 *
 *  Created on: 14 Aug 2016
 *      Author: radu
 */


#include "any_afse.h"
#include "graph.h"

// Copy DaoOpt Function class into mex::Factor class structures
mex::vector<mex::Factor> AnyAFSE::copyFactors(void) {
	Problem* _p = m_problem.get();
	mex::vector<mex::Factor> fs(_p->getC());
	for (int i = 0; i < _p->getC(); ++i)
		fs[i] = _p->getFunctions()[i]->asFactor();
	return fs;
}

void AnyAFSE::init() {

	// safety checks
	assert(m_options != NULL);
	assert(m_options->inputFile.empty() == false);

	// start a timer
	Timer tm;
	tm.start();

	m_solved = false;

	// load problem file
	m_problem.reset(new Problem());
	if (m_options->format == FORMAT_UAI) {
		m_problem->parseUAI(m_options->inputFile, m_options->evidenceFile, m_options->positive);
		cout << "Created problem with " << m_problem->getN() << " variables and "
				<< m_problem->getC() << " functions." << endl;
	} else if (m_options->format == FORMAT_ERGO) {
		m_problem->parseERG(m_options->inputFile, m_options->evidenceFile, m_options->positive);
		cout << "Created problem with " << m_problem->getN() << " variables and "
				<< m_problem->getC() << " functions." << endl;
	}

	// remove evidence
	m_problem->removeEvidence();
	cout << "Removed evidence, now " << m_problem->getN() << " variables and "
			<< m_problem->getC() << " functions." << endl;

	// some statistics
	cout << "Global constant:\t" << m_problem->getGlobalConstant() << endl;
	cout << "Max. domain size:\t" << (int) m_problem->getK() << endl;
	cout << "Max. function arity:\t" << m_problem->getR() << endl;

	// create primal graph
	Graph g(m_problem->getN());
	const std::vector<Function*>& fns = m_problem->getFunctions();
	for (std::vector<Function*>::const_iterator it = fns.begin();
			it != fns.end(); ++it) {
		g.addClique((*it)->getScope());
	}
	std::cout << "Graph with " << g.getNumNodes() << " nodes and "
			<< g.getNumEdges() << " edges created." << std::endl;

	// read the MAP (query) variables
	std::vector<int> mapVars;
	m_problem->parseMapVars(m_options->mapFile, mapVars);
	m_problem->updateVartypes(mapVars);

	cout << "Read MAP variables from file "
		<< m_options->mapFile << " (" << mapVars.size() << ")." << endl;
	cout << "MAP: ";
	std::copy(mapVars.begin(), mapVars.end(), std::ostream_iterator<int>(cout, " "));
	cout << endl;

	// read variable ordering (constrained)
	int w = g.eliminate(m_elim); // generate the min-fill ordering

	std::cout << "Generated elimination ordering (" << w << "): ";
	std::copy(m_elim.begin(), m_elim.end(), std::ostream_iterator<int>(std::cout, " "));
	std::cout << std::endl;

	// build pseudo tree (wrt to constrained ordering)
	m_pseudotree.reset(new Pseudotree(m_problem.get(), m_options->subprobOrder));
	m_pseudotree->build(g, m_elim, m_options->cbound);
	int h = m_pseudotree->getHeight();

	// pseudo tree has dummy node after build(), add to problem
	m_problem->addDummy(); // add dummy variable to problem, to be in sync with pseudo tree
	m_pseudotree->resetFunctionInfo(m_problem->getFunctions());

	cout << "Induced width:      " << w << std::endl;
	cout << "Pseudo tree height: " << h << std::endl;
	cout << "Problem variables:  " << m_problem->getN() << std::endl;

	// stop the timer and record the preprocessing time
	tm.stop();
	m_tmLoad = tm.elapsed();

	// Create pseudo tree and connect disconnected components with dummy functions
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

}

int AnyAFSE::solve() {

	// prologue
	std::cout << "--- Starting AFSE ---" << std::endl;
	int res = SEARCH_SUCCESS;
	char* EmergencyMemory = new char[16384]; // in order to handle bad alloc

	// init search (timer, search space, assignment)
	m_tmSolve = 0;
	m_timer.reset(m_options->timeLimit);
	m_timer.start();
	m_solutionCost = ELEM_NAN;
	m_assignment.resize(m_problem->getN(), UNKNOWN);

	try {

		// Mark MAP variables
		mex::vector<bool> vtypes(m_problem->getN(), false);
		for (int var = 0; var < m_problem->getN(); ++var) {
			if (m_problem->isMap(var))
				vtypes[var] = true;
		}

		// init the afse solver
		m_afse = mex::afse(copyFactors());
		m_afse.setProperties();
		m_afse.setVarTypes(vtypes);
		m_afse.init(); // re-compute the elim order
		m_iterations = m_afse.solve(m_options->timeLimit, m_options->afseStep); // throw timelimit, out-of-memory

		m_solutionCost = m_afse.logZ(); // optimal solution found

		m_solved = true;
	} catch (std::bad_alloc& ba) {
		delete[] EmergencyMemory;
		EmergencyMemory = NULL;
		std::cout << "OUT OF MEMORY";
		res = SEARCH_OUT_OF_MEMORY;
	}

	// stop timer
	m_timer.stop();
	m_tmSolve = m_timer.elapsed();
	m_solutionCost += -ELEM_ENCODE(m_problem->getGlobalConstant());

	// output solution (if found)
	std::cout << std::endl << std::endl;
	std::cout << "--- AFSE done ---" << std::endl;
	std::cout << "Problem name:        " << m_problem->getName() << std::endl;
	std::cout << "Status:              " << solver_status[res] << std::endl;
	std::cout << "Iterations:          " << m_iterations << std::endl;
	std::cout << "Time elapsed:        " << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
	std::cout << "Preprocessing:       " << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
	std::cout << "Lower bound:         " << ELEM_DECODE(m_afse.logZlb()) << " (" << -m_afse.logZlb() << ")" << std::endl;
	std::cout << "Upper bound:         " << ELEM_DECODE(m_afse.logZub()) << " (" << -m_afse.logZub() << ")" << std::endl;
	std::cout << "Solution:            " << ELEM_DECODE(m_solutionCost) << " (" << -m_solutionCost  << ")" << std::endl;
	std::cout << "-------------------------------" << std::endl;

	// clean up
	if (EmergencyMemory) delete[] EmergencyMemory;

	return res;
}
