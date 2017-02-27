/*
 * map_generator.cpp
 *
 *  Created on: Sep 27, 2013
 *      Author: radu
 */


#include "map_generator.h"
#include "graph.h"
#include "pseudo_tree.h"

void MapGenerator::init() {

}

int MapGenerator::solve() {

	// safety checks
	assert(m_options != NULL);
	assert(m_options->inputFile.empty() == false);

	// load problem file
	m_problem.reset(new Problem());
	if (m_options->format == FORMAT_UAI) {
		m_problem->parseUAI(m_options->inputFile, "", m_options->positive);
		cout << "Created problem with " << m_problem->getN() << " variables and "
				<< m_problem->getC() << " functions." << endl;
	} else if (m_options->format == FORMAT_ERGO) {
		m_problem->parseERG(m_options->inputFile, "", m_options->positive);
		cout << "Created problem with " << m_problem->getN() << " variables and "
				<< m_problem->getC() << " functions." << endl;
		std::string uaiFileName = m_problem->getName();
		uaiFileName += ".uai";
		m_problem->writeUAI(uaiFileName);
	}

	// some statistics
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

	if (m_options->orderingFile.empty()) {

		// instance 0: find unconstrained minfill ordering and select last M
		std::cout << "Generating MAP instance " << 0 << " ..." << std::endl;

		// select the random MAP variables
		int N = m_problem->getN();
		double perc = m_options->percMapVars;
		int M = (int) std::floor((double)N * perc);
		if (M == 0) M = 10;
		if (M > N) M = N;

		std::vector<int> elim0;
		int w = g.eliminate(elim0); // unconstrained minfill

		// create the pseudo tree corresponding to this ordering
		std::vector<int> temp;
		Pseudotree pt(m_problem.get(), m_options->subprobOrder);
		pt.build(g, elim0, m_options->cbound);
		cout << "Read constrained elimination ordering from file "
				<< m_options->orderingFile << " (" << w << '/'
				<< pt.getHeight() << ")." << endl;
		pt.bfs(temp);
		std::vector<int> mapVars0;
		for (std::vector<int>::iterator it = temp.begin();
				it != temp.end(); ++it) {
			int var = *it;
			mapVars0.push_back(var);
			if ((int)mapVars0.size() == M) {
				break;
			}
		}

//		std::set<int> mapVars0;
//		for (std::vector<int>::reverse_iterator it = elim0.rbegin();
//				it != elim0.rend(); ++it) {
//			int var = *it;
//			mapVars0.insert(var);
//			if ((int)mapVars0.size() == M) {
//				break;
//			}
//		}


		std::cout << " - selected " << mapVars0.size() << " random MAP variables." << std::endl;
		std::cout << " - generated MAP elimination order (" << w << "/" << N << ")." << std::endl;
		std::cout << " - elimination order size: " << elim0.size() << std::endl;
		std::cout << " - nodes in the graph: " << g.getNumNodes() << std::endl;
		std::cout.flush();


		// elimination order file
		char version[128];
		std::string elimOrderFile;
		elimOrderFile = m_problem->getName();
		sprintf(version, ".vo.%d", 0);
		elimOrderFile += std::string(version);
		std::cout << " - saved MAP elimination order to " << elimOrderFile << "." << std::endl;
		std::ofstream out1(elimOrderFile.c_str());
		out1 << "# constrained elimination order for " << m_problem->getName()
				<< "; width = " << w << std::endl;
		out1 << elim0.size() << std::endl;
		std::copy(elim0.begin(), elim0.end(), std::ostream_iterator<int>(out1, "\n"));
		out1.close();

		std::string mapFile;
		mapFile = m_problem->getName();
		sprintf(version, ".map.%d", 0);
		mapFile += std::string(version);
		std::cout << " - saved MAP variables to " << mapFile << "." << std::endl;
		std::ofstream out2(mapFile.c_str());
		out2 << "# MAP variables for " << m_problem->getName()
				<< "; count = " << M << std::endl;
		out2 << mapVars0.size() << std::endl;
		std::copy(mapVars0.begin(), mapVars0.end(), std::ostream_iterator<int>(out2, "\n"));
		out2.close();

		// generate random MAP problem instances (instances 1, 2, ...)
		for (int inst = 1; inst <= m_options->numInstances; ++inst) {
			std::cout << "Generating MAP instance " << inst << " ..." << std::endl;

			// select the M variables randomly
			std::set<int> mapVars;
			for (int i = 0; i < M;) {
				int var = rand::next(N);
				if (mapVars.find(var) == mapVars.end()) {
					mapVars.insert(var);
					++i;
				}
			}

			std::cout << " - selected " << M << " random MAP variables." << std::endl;

			std::vector<int> elim;
			int w = g.eliminate(mapVars, elim); // constrained minfill
			std::cout << " - generated MAP elimination order (" << w << "/" << N << ")." << std::endl;
			std::cout << " - elimination order size: " << elim.size() << std::endl;
			std::cout << " - nodes in the graph: " << g.getNumNodes() << std::endl;
			std::cout.flush();

			// elimination order file
			elimOrderFile = m_problem->getName();
			if (m_options->numInstances == 1) sprintf(version, ".vo");
			else sprintf(version, ".vo.%d", inst);
			elimOrderFile += std::string(version);
			std::cout << " - saved MAP elimination order to " << elimOrderFile << "." << std::endl;
			std::ofstream outElim(elimOrderFile.c_str());
			outElim << "# constrained elimination order for " << m_problem->getName()
					<< "; width = " << w << std::endl;
			outElim << elim.size() << std::endl;
			std::copy(elim.begin(), elim.end(), std::ostream_iterator<int>(outElim, "\n"));
			outElim.close();

			// MAP file
			mapFile = m_problem->getName();
			if (m_options->numInstances == 1) sprintf(version, ".map");
			else sprintf(version, ".map.%d", inst);
			mapFile += std::string(version);
			std::cout << " - saved MAP variables to " << mapFile << "." << std::endl;
			std::ofstream outMAP(mapFile.c_str());
			outMAP << "# MAP variables for " << m_problem->getName()
					<< "; count = " << M << std::endl;
			outMAP << mapVars.size() << std::endl;
			std::copy(mapVars.begin(), mapVars.end(), std::ostream_iterator<int>(outMAP, "\n"));
			outMAP.close();
		}

		// instance r: MAP variables are the root variables
		std::cout << "Generating MAP instance " << "R" << " ..." << std::endl;

		// select the root MAP variables
		std::set<int> mapVarsR;
		N = m_problem->getN();
		for (std::vector<Function*>::const_iterator it = m_problem->getFunctions().begin();
				it != m_problem->getFunctions().end(); ++it) {
			Function* f = (*it);
			const std::set<int>& scope = f->getScope();
			if (scope.size() == 1) { // root variable
				int var = *(scope.begin());
				mapVarsR.insert(var);
			}
		}

		M = (int) mapVarsR.size();

		if (M == 0) {
			std::cout << "No root MAP variables selected!" << std::endl;
		} else {

			std::vector<int> elimR;
			w = g.eliminate(mapVarsR, elimR); // constrained minfill
			std::cout << " - selected " << M << " root MAP variables." << std::endl;
			std::cout << " - generated MAP elimination order (" << w << "/" << N << ")." << std::endl;
			std::cout << " - elimination order size: " << elimR.size() << std::endl;
			std::cout << " - nodes in the graph: " << g.getNumNodes() << std::endl;
			std::cout.flush();

			// elimination order file
			elimOrderFile = m_problem->getName();
			sprintf(version, ".vo.R");
			elimOrderFile += std::string(version);
			std::cout << " - saved MAP elimination order to " << elimOrderFile << "." << std::endl;
			std::ofstream outElimR(elimOrderFile.c_str());
			outElimR << "# constrained elimination order for " << m_problem->getName()
					<< "; width = " << w << std::endl;
			outElimR << elimR.size() << std::endl;
			std::copy(elimR.begin(), elimR.end(), std::ostream_iterator<int>(outElimR, "\n"));
			outElimR.close();

			// MAP file
			mapFile = m_problem->getName();
			sprintf(version, ".map.R");
			mapFile += std::string(version);
			std::cout << " - saved MAP variables to " << mapFile << "." << std::endl;
			std::ofstream outMAPR(mapFile.c_str());
			outMAPR << "# MAP variables for " << m_problem->getName()
					<< "; count = " << M << std::endl;
			outMAPR << mapVarsR.size() << std::endl;
			std::copy(mapVarsR.begin(), mapVarsR.end(), std::ostream_iterator<int>(outMAPR, "\n"));
			outMAPR.close();
		}

		// instance h: MAP variables form a start hypergraph pseudo tree
		std::cout << "Generating MAP instance " << "H" << " ..." << std::endl;

		// select the root MAP variables
		if (m_options->mapFile.empty()) {
			std::cout << " - nothing to generate." << std::endl;
		} else {

			// read the MAP vars from the file
			std::vector<int> temp;
			m_problem->loadMapVars(m_options->mapFile, temp);
			std::set<int> mapVarsH(temp.begin(), temp.end());

			std::vector<int> elim;
			int w = g.eliminate(mapVarsH, elim); // constrained minfill
			std::cout << " - number of MAP varibles: " << mapVarsH.size() << std::endl;
			std::cout << " - generated MAP elimination order (" << w << "/" << N << ")." << std::endl;
			std::cout << " - elimination order size: " << elim.size() << std::endl;
			std::cout << " - nodes in the graph: " << g.getNumNodes() << std::endl;
			std::cout.flush();

			// elimination order file
			elimOrderFile = m_problem->getName();
			sprintf(version, ".vo.H");
			elimOrderFile += std::string(version);
			std::cout << " - saved MAP elimination order to " << elimOrderFile << "." << std::endl;
			std::ofstream outElim(elimOrderFile.c_str());
			outElim << "# constrained elimination order for " << m_problem->getName()
					<< "; width = " << w << std::endl;
			outElim << elim.size() << std::endl;
			std::copy(elim.begin(), elim.end(), std::ostream_iterator<int>(outElim, "\n"));
			outElim.close();

		}

	} else { // this is a standard minfill or a hypergraph ordering
		std::vector<int> elim;
		m_problem->parseOrdering(m_options->orderingFile, elim);

		// select the MAP variables to be the last in the ordering
		int N = m_problem->getN();
		double perc = m_options->percMapVars;
		int M = std::floor((double)N * perc);
		if (M == 0) M = 10;
		if (M > N) M = N;
		std::vector<int> temp;

		// create the pseudo tree corresponding to this ordering
		Pseudotree pt(m_problem.get(), m_options->subprobOrder);
		pt.build(g, elim, m_options->cbound);
		int w = pt.getWidth();
		cout << "Read constrained elimination ordering from file "
				<< m_options->orderingFile << " (" << w << '/'
				<< pt.getHeight() << ")." << endl;
		pt.bfs(temp);

		std::vector<int> mapVars;
		for (std::vector<int>::iterator it = temp.begin();
				it != temp.end(); ++it) {
			int var = *it;
			mapVars.push_back(var);
			if ((int)mapVars.size() == M) {
				break;
			}
		}

		char version[128];
		std::string mapFile;
		mapFile = m_problem->getName();
		sprintf(version, ".map");
		mapFile += std::string(version);
		std::cout << " - saved MAP variables to " << mapFile << "." << std::endl;
		std::ofstream out2(mapFile.c_str());
		out2 << "# MAP variables for " << m_problem->getName()
				<< "; count = " << M << std::endl;
		out2 << mapVars.size() << std::endl;
		std::copy(mapVars.begin(), mapVars.end(), std::ostream_iterator<int>(out2, "\n"));
		out2.close();
	}

	return 0;
}

