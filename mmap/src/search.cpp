/*
 * search.cpp
 *
 *  Created on: Apr 17, 2013
 *      Author: radu
 */


#include "search.h"
#include "timer.h"
#include "hash_murmur.h"
#include "function.h"

size_t strHasher(const std::string& a) {
#ifdef HASH_MURMUR
	size_t hash = 0;
	size_t seed = 0x9e3779b9;
	int len = (int)a.size();
	register unsigned char *key = (unsigned char*)a.c_str();
	MurmurHash3_x64_64(key, len, seed, &hash);
	return hash;
#elif defined(HASH_DEFAULT)
	return boost::hash<std::string>()(a.c_str());
#endif
}

using namespace minisat;

// encode the zero probability tuples into clauses
void Search::encodeDeterminism() {

	if (m_options->propagation == CP_FULL) { // full SAT (minisat)

		vec<Lit> lits;

		// create the SAT variables V = (var, val)
		int N = m_problem->getN();
		m_var2sat.resize(N);
		for (int var = 0; var < N; ++var) {
			m_var2sat[var].resize(m_problem->getDomainSize(var));
			for (int val = 0; val < m_problem->getDomainSize(var); ++val) {
				Var vid = m_minisat.newVar();
				m_var2sat[var][val] = vid;
			}
		}

		std::vector<int> assignment;
		assignment.resize(m_problem->getN(), UNKNOWN);

		// add nogood clauses
		std::vector<Function*>::const_iterator fi = m_problem->getFunctions().begin();
		for (; fi != m_problem->getFunctions().end(); ++fi) {
			Function* fn = (*fi);

			std::vector<int> argv, vals;
			std::set<int>::const_iterator si = fn->getScope().begin();
			for (; si != fn->getScope().end(); ++si) {
				argv.push_back(*si);
			}
			vals.resize(argv.size(), 0);
			vals[vals.size()-1] = -1;

			int i, argc = (int)argv.size();
			while (true) {
				for (i = argc-1; i >= 0; --i) {
					int arg = argv[i];
					int last = m_problem->getDomainSize(arg) - 1;
					if (vals[i] < last) break;
					vals[i] = 0;
				}

				if (i < 0) break;	// done.
				++vals[i];

				// Full scope instantiation.
				for (int a = 0; a < argc; ++a) assignment[argv[a]] = vals[a];

				// Add a nogood clause.
				double val = fn->getValue(assignment);
				if (0.0 == val) {
					lits.clear();
					for (int a = 0; a < argc; ++a) {
						int vid = m_var2sat[argv[a]][vals[a]];
						lits.push(~mkLit(vid));
					}

					m_minisat.addClause_(lits);
				}
			}

			assignment.resize(m_problem->getN(), UNKNOWN);
		}

		// add at-most-one clauses
		lits.clear();
		for (int var = 0; var < m_problem->getN(); ++var) {
			int domSize = m_problem->getDomainSize(var);
			for (int val1 = 0; val1 < domSize-1; ++val1) {
				for (int val2 = val1 + 1; val2 < domSize; ++val2) {
					int v1 = m_var2sat[var][val1];
					int v2 = m_var2sat[var][val2];

					lits.clear();
					lits.push(~mkLit(v1));
					lits.push(~mkLit(v2));
					m_minisat.addClause_(lits);
				}
			}
		}

		// add at-least-one clauses
		lits.clear();
		for (int var = 0; var < m_problem->getN(); ++var) {
			lits.clear();
			for (int val = 0; val < m_problem->getDomainSize(var); ++val) {
				int vid = m_var2sat[var][val];
				lits.push(mkLit(vid));
			}
			m_minisat.addClause_(lits);
		}

		m_minisat.simplify();
		std::cout << "Created SAT (minisat) instance with " << m_minisat.nVars()
				<< " variables and " << m_minisat.nClauses() << " clauses" << std::endl;

		bool ok = m_minisat.solve();
		if (ok == true) {
			std::cout << "Initial SAT instance is satisfiable!" << std::endl;
		} else {
			std::cout << "Initial SAT instance is not satisfiable!" << std::endl;
		}

		assert(m_currentDomains.empty());

	} else if (m_options->propagation == CP_UNIT) { // unit prop (zchaff)

		int vid = 1; // SAT vars are indexed from 1
		int N = m_problem->getN(), nvars = 0;
		for (int i = 0; i < N; ++i) {
			nvars += m_problem->getDomainSize(i);
		}

		// Create variables from the belief network. For each variable
		// value assignment (e.g. X=x) create a sat variable X_x.
		m_var2sat.resize(N);
		m_sat2var.resize(nvars+1, std::make_pair(-1, -1));
		m_zchaff.set_variable_number(nvars);
		for (int i = 0; i < N; ++i)	{
			int dom = m_problem->getDomainSize(i);
			m_var2sat[i].resize(dom, 0);
			for (int k = 0; k < dom; ++k)	{
				m_var2sat[i][k] = vid;
				m_sat2var[vid] = std::make_pair(i, k);
				++vid;
			}
		}

		// Set the sat vars in the solver
		m_zchaff.set_var_maps(m_var2sat, m_sat2var);

		// From clauses from the 0 probability tuples (conflict clauses).
		int num_cl = 0;
		std::vector<Function*>::const_iterator fi = m_problem->getFunctions().begin();
		for (; fi != m_problem->getFunctions().end(); ++fi) {

			Function *fn = (*fi);

			// Get the scope.
			int argc = fn->getArity();
			const std::set<int>& scope = fn->getScope();
			int *argv = new int[argc];
			int l = 0;
			for (std::set<int>::const_iterator si = scope.begin();
					si != scope.end(); ++si) {
				argv[l++] = *si;
			}

			int *vals = new int[argc];
			memset(vals, 0, argc*sizeof(int));
			vals[argc-1] = -1;

			// Enumerate scope instantiations.
			std::vector<int> assignment(N, -1);
			int i;
			while (true) {
				for (i = argc-1; i >= 0; --i) {
					int arg = argv[i];
					int last = m_problem->getDomainSize(arg) - 1;
					if (vals[i] < last) break;
					vals[i] = 0;
				}

				if (i < 0) break;	// done.
				++vals[i];

				// Full scope instantiation.
				for (int a = 0; a < argc; ++a) {
					assignment[argv[a]] = vals[a];
				}

				// Add a nogood clause.
				if (fn->getValue(assignment) == 0.0) {
					m_zchaff.add_nogood_clause(argc, argv, vals);
					++num_cl;
				}
			}

			// Restore assignment.
			assignment.resize(N, -1);
			delete[] vals;
			delete[] argv;
		}

		// Add at-most-one clauses for the domain values (to ensure that
		// no BN variable is assigned more than one values at once).
		for (int i = 0; i < N; ++i) {
			int dom = m_problem->getDomainSize(i);
			for (int v1 = 0; v1 < dom-1; ++v1) {
				for (int v2 = v1+1; v2 < dom; ++v2) {
					m_zchaff.add_atmostone_clause(i, v1, v2);
					++num_cl;
				}
			}
		}

		// Add at-least-one clauses for the domain values (to ensure that
		// each BN variable is assigned at least one value from its domain).
		for (int i = 0; i < N; ++i) {
			int dom = m_problem->getDomainSize(i);
			m_zchaff.add_atleastone_clause(i, dom);
			++num_cl;
		}

		// Init the SAT solver and preprocess.
		m_zchaff.init_solve();

		// Preprocess.
		m_zchaff.preprocess();

		std::cout << "Created SAT (zchaff) instance with " << m_zchaff.num_variables()
				<< " variables and " << m_zchaff.num_clauses() << " clauses" << std::endl;

		// Init current domains
		m_currentDomains.resize(N);
		for (int i = 0; i < N; ++i) {
			m_currentDomains[i].resize(m_problem->getDomainSize(i), true);
		}

		// Update domains
		std::list<std::pair<int, int> > changes;
		m_zchaff.update_domains(m_currentDomains, changes);
	} else {
		assert(m_currentDomains.empty());
	}
}

// returns 'true' if conflict, otherwise returns 'false'
bool Search::propagate(const std::vector<int>& vars) {

	// safety checks
	assert( m_options->propagation == CP_FULL );
	vec<Lit> lits;
	for (std::vector<int>::const_iterator it = vars.begin();
			it != vars.end(); ++it) {
		int var = (*it);
		lits.push(mkLit(var));
	}

	// propagate the assignment (full SAT solving)
	bool ok = m_minisat.solve(lits);
	if ( ok == true ) return false; // no conflict
	else return true; // conflict
}

// returns 'false' if conflict, otherwise returns 'true'
bool Search::lookahead(int var, int val, list<pair<int, int> >& changes,
		const std::set<int>& subtree) {

	// safety checks
	assert( m_options->propagation == CP_UNIT );

	// Unit resolution.
	int res = m_zchaff.propagate(var, val);
	m_zchaff.update_domains(m_currentDomains, changes);
	if (zchaff::CONFLICT == res) {
		return false;
	}

	// Look for a future empty domain in the subproblem.
	std::set<int>::const_iterator si = subtree.begin();
	for (; si != subtree.end(); ++si) {
		if (*si == var) {
			continue;
		}

		//assert(!isAssigned(subtree[i]));
		int desc = *si;;
		std::vector<bool> &domain = m_currentDomains[desc];
		size_t sz = domain.size(), pruned = 0;
		for (size_t j = 0; j < sz; ++j) {
			if (!domain[j]) {
				++pruned;
			}
		}

		if (pruned == sz) {
			return false;	// future empty domain detected.
		}
	}

	return true;
}


// initialize the solver
void Search::init() {

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

	// check if GG is disconnected
//	std::set<int> vars;
//	Graph GGG(g);
//	for (int i = 0; i < m_problem->getN(); ++i) vars.insert(i);
//	std::map<int, std::set<int> > comps1 = GGG.connectedComponents(vars);
//	std::cout << "Input graph connected components: " << comps1.size();
//	std::cout << std::endl;


	// read the MAP (query) variables
	std::vector<int> mapVars;
	if (m_options->algorithm != ALGO_EVALUATOR)
		m_problem->parseMapVars(m_options->mapFile, mapVars);
	m_problem->updateVartypes(mapVars);

	cout << "Read MAP variables from file "
		<< m_options->mapFile << " (" << mapVars.size() << ")." << endl;
	cout << "MAP: ";
	std::copy(mapVars.begin(), mapVars.end(), std::ostream_iterator<int>(cout, " "));
	cout << endl;

	// read variable ordering (constrained)
#ifndef UAI_COMPETITION
	std::vector<int> elim;
	int w; // = numeric_limits<int>::max();
	if (m_options->orderingFile.empty() == false) {
		m_problem->parseOrdering(m_options->orderingFile, elim);
		w = g.triangulate(elim);
	} else {
		// generate the constrained min-fill ordering
		std::set<int> temp(mapVars.begin(), mapVars.end());
		w = g.eliminate(temp, elim);
	}

	std::cout << "Generated constrained elimination ordering (" << w << "): ";
	std::copy(elim.begin(), elim.end(), std::ostream_iterator<int>(std::cout, " "));
	std::cout << std::endl;

#else
	std::vector<int> elim, best;
	int w = numeric_limits<int>::max();

	// generate the constrained min-fill ordering
	Timer tmVO(m_options->orderingTime);
	tmVO.start();
	std::cout << "Randomized minfill ordering ..." << std::endl;
	std::cout << "  iter: " << m_options->orderingIter << "; time: " << m_options->orderingTime << std::endl;
	std::set<int> temp(mapVars.begin(), mapVars.end());
	int stopIter = m_options->orderingIter;
	for (int iter = 1; iter <= stopIter; ++iter) {
		int ww = g.eliminate(temp, elim);
		if (ww < w) {
			w = ww;
			best = elim;
			std::cout << "  found best constrained induced width: " << w << std::endl;
		}

		if (tmVO.timeout()) {
			break;
		}
	}

	elim = best;
	std::cout << "Generated constrained elimination ordering (" << w << "): ";
	std::copy(elim.begin(), elim.end(), std::ostream_iterator<int>(std::cout, " "));
	std::cout << std::endl;

	// save the elimination order
	m_elimOrder = elim;
#endif

	// build pseudo tree (wrt to constrained ordering)
	m_pseudotree.reset(new Pseudotree(m_problem.get(), m_options->subprobOrder));
	m_pseudotree->build(g, elim, m_options->cbound);
	w = m_pseudotree->getWidth();
	if (m_options->orderingFile.empty() == false) {
		std::cout << "Read constrained elimination ordering from file "
				<< m_options->orderingFile << " (" << w << '/'
				<< m_pseudotree->getHeight() << ")." << endl;
	} else {
		std::cout << "Generated constrained elimination (" << w << '/'
				<< m_pseudotree->getHeight() << ")." << endl;
	}

	// pseudo tree has dummy node after build(), add to problem
	m_problem->addDummy(); // add dummy variable to problem, to be in sync with pseudo tree
	m_pseudotree->resetFunctionInfo(m_problem->getFunctions());

	// get the start (MAP) pseudo tree depth
	int mapDepth = 0;
	for (size_t j = 0; j < m_pseudotree->getN(); ++j) {
		PseudotreeNode* node = m_pseudotree->getNode(j);
		if (m_problem->isMap(j)) {
			mapDepth = std::max(mapDepth, node->getDepth());
		}
	}

	cout << "Induced width:\t\t" << m_pseudotree->getWidth() << std::endl;
	cout << "Pseudotree depth:\t" << m_pseudotree->getHeight() << std::endl;
	cout << "Start pseudotree:\t" << mapDepth << std::endl;
	cout << "Problem variables:\t" << m_pseudotree->getSize() << std::endl;
	cout << "Disconn. components:\t" << m_pseudotree->getComponents() << std::endl;
	if (m_problem->getDeterminism() < DETERMINISM_THRESHOLD) {
		m_options->propagation = CP_NONE; // disable constraint propagation
		cout << "Determinism below threshold (" << DETERMINISM_THRESHOLD
				<< "), so unit propagation is disabled." << std::endl;
	}

	// output the pseudo tree
	if (m_options->pseudotreeFile.empty() == false) {
		m_pseudotree->outputToFile(m_options->pseudotreeFile);
	}

	// stop the timer and record the preprocessing time
	tm.stop();
	m_tmLoad = tm.elapsed();

}

