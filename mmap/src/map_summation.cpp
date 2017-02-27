/*
 * evaluator.cpp
 *
 *  Created on: Oct 22, 2014
 *      Author: radu
 */

#include "map_summation.h"
#include "timer.h"
#include "graph.h"

// initialize the solver
void MapSummation::init() {

	// safety checks
	assert(m_options != NULL);
	assert(m_options->inputFile.empty() == false);

	// start a timer
	Timer tm;
	tm.start();

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
	//m_problem->removeEvidence();
	//cout << "Removed evidence, now " << m_problem->getN() << " variables and "
	//		<< m_problem->getC() << " functions." << endl;

	// some statistics
	//cout << "Global constant:\t" << m_problem->getGlobalConstant() << endl;
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

	//g.outputToFile("graph.dot");
	//std::cout << "Written to file." << std::endl;

	// check if GG is disconnected
//	std::set<int> vars;
//	Graph GGG(g);
//	for (int i = 0; i < m_problem->getN(); ++i) vars.insert(i);
//	std::map<int, std::set<int> > comps1 = GGG.connectedComponents(vars);
//	std::cout << "Input graph connected components: " << comps1.size();
//	std::cout << std::endl;


	// read the MAP (query) variables
	std::vector<int> mapVars;
	m_problem->loadMapVars(m_options->mapFile, mapVars);
	m_problem->updateVartypes(mapVars);

	std::cout << "Read MAP variables from file "
		<< m_options->mapFile << " (" << mapVars.size() << ")." << endl;
	std::cout << "MAP: ";
	std::copy(mapVars.begin(), mapVars.end(), std::ostream_iterator<int>(cout, " "));
	std::cout << endl;

	// read variable ordering (constrained)
	std::vector<int> elim;
//	int w = numeric_limits<int>::max();
//	m_problem->parseOrdering(m_options->orderingFile, elim);

	// generate the constrained min-fill ordering
	std::set<int> temp(mapVars.begin(), mapVars.end());
	int w = g.eliminate(temp, elim);
	std::cout << "Constrained elimination ordering (" << w << "): ";
	std::copy(elim.begin(), elim.end(), std::ostream_iterator<int>(std::cout, " "));
	std::cout << std::endl;
	std::cout << "Induced width:          " << w << std::endl;
	std::cout << "Problem variables:      " << m_problem->getN() << std::endl;

	// save the ordered sum variables (original indeces)
	std::vector<int> orderedSUM;
	for (size_t i = 0; i < elim.size(); ++i) {
		if (m_problem->isSum(elim[i]))
			orderedSUM.push_back(elim[i]);
	}

    // Remove MAP variables
    int N = m_problem->getN();
    std::vector<int> assignment;
    assignment.resize(N, UNKNOWN);
    for (size_t i = 0; i < mapVars.size(); ++i) {
    	assignment[ mapVars[i] ] = 0;
    }
    m_problem->setEvidence(assignment);
    m_problem->removeEvidence();
    Graph gg(m_problem->getN());
	const std::vector<Function*>& funs = m_problem->getFunctions();
	for (std::vector<Function*>::const_iterator it = funs.begin();
			it != funs.end(); ++it) {
		gg.addClique((*it)->getScope());
	}
	std::cout << "SUM graph with " << gg.getNumNodes() << " nodes and "
			<< gg.getNumEdges() << " edges created." << std::endl;

	//gg.outputToFile("graph2.dot");
	//std::cout << "Written to file." << std::endl;

	// Get the elimination order of the SUM variables from the original order
	// This will be consistent with the search order along the conditioned
	// pseudo tree (used by AOBB, AOBF, RBFAOO, and mostly by the AO search).
	const std::map<int, int>& old2new = m_problem->getOld2New();
	std::vector<int> elimS, elimCS;
	for (size_t i = 0; i < orderedSUM.size(); ++i) {
		int var = orderedSUM[i];
		std::map<int, int>::const_iterator mi = old2new.find(var);
		int nvar = mi->second;
		elimCS.push_back(nvar);
	}
	int ws = gg.eliminate(elimS);
	int wcs = gg.triangulate(elimCS);
	std::cout << "SUM induced width:                   " << ws << std::endl;
	std::cout << "SUM consistent induced width:        " << wcs << std::endl;
	std::cout << "SUM problem variables:               " << m_problem->getN() << std::endl;
	std::cout << "Done." << std::endl;

	// stop the timer and record the preprocessing time
	tm.stop();
	m_tmLoad = tm.elapsed();

}


int MapSummation::solve() {

	return SEARCH_SUCCESS;
}
