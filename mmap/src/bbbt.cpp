/*
 * bbbt.cpp
 *
 *  Created on: Nov 21, 2013
 *      Author: radu
 */


#include "bbbt.h"
#include "mini_cluster_tree_heur.h"
#include "mini_cluster_tree_incr_heur.h"

#include "mini_bucket_heur.h"
#include "mini_bucket_mm_heur.h"

// initialize
void BBBT::init() {

	// initialize
	Search::init();

#ifdef UAI_COMPETITION
	// preprocessing by mixed BP
	MiniBucketMMHeur* pre = new MiniBucketMMHeur(m_problem.get(),
			m_pseudotree.get(), m_options, m_options->ibound);
	delete pre;
#endif

	// safety checks
	assert(m_options->heuristic == HEUR_MINI_CLUSTER_TREE);

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
	tm.start();

	// encode determinism
	encodeDeterminism();

#ifdef UAI_COMPETITION
	// adjust the i-bound to fit the memory limit
	// determine the optimal i-bound given 4GB of RAM
	size_t memlimit = m_options->memoryLimit*1024;
	MiniBucketHeur *temp = new MiniBucketHeur(m_problem.get(),
			m_pseudotree.get(), m_options, m_options->ibound);
	temp->limitSize(memlimit, &m_assignment);
	delete temp;
#endif


	// init heuristic generator
	if (m_options->incremental == false) {
		m_useStaticVo = false;
		std::cout << "Computing the mini-cluster-tree heuristic..." << std::endl;
		m_heuristic.reset(new MiniClusterTreeHeur(m_problem.get(),
				m_pseudotree.get(), m_options));
	}  else {
		m_useStaticVo = true;
		std::cout << "Computing the incremental mini-cluster-tree heuristic..." << std::endl;
		m_heuristic.reset(new MiniClusterTreeIncrHeur(m_problem.get(),
				m_pseudotree.get(), m_options));
	}

	size_t sz = m_heuristic->build(&m_assignment, true);
	tm.stop();
	m_tmHeuristic = tm.elapsed();

	std::cout << "\tCluster tree finished in " << m_tmHeuristic
			<< " seconds" << std::endl;
	std::cout << "\tUsing " << (sz / (1024*1024.0)) * sizeof(double)
			<< " MBytes of RAM" << std::endl;

	// Now, we need to rearrange the pseudo tree such that the MAP variables
	// form a chain at the top of the pseudo tree (and they're searched in the
	// order suggested by the bucket tree heuristic ...
	std::list<int> mapSearchOrder;
	m_heuristic->getStaticSearchOrder(mapSearchOrder);
	m_pseudotree->forceMapChain(mapSearchOrder, m_options->cbound);
	m_lastMapVar = mapSearchOrder.back();
	m_mapVars = mapSearchOrder;
	assert( m_problem->isMap(m_lastMapVar) );

	// heuristic might have changed problem functions, pseudotree needs remapping
	m_pseudotree->resetFunctionInfo(m_problem->getFunctions());

	// Record time after initialization
	cout << "Last MAP variable: " << m_lastMapVar << std::endl;
	cout << "Initialization complete: " << (m_tmLoad + m_tmHeuristic)
			<< " seconds" << std::endl;

	// check for upper bound computation only
	if (m_options->upperBoundOnly) {
		std::cout << "Done." << std::endl;
		exit(EXIT_SUCCESS);
	}

}

// select next unassigned MAP variable
int BBBT::selectNextVar(const set<int>& M,
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
void BBBT::bbDynamic(std::set<int> M, std::vector<int> assignment, std::set<int> F) {

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
				<< std::setw(12) << m_lowerBound << " (" << -ELEM_ENCODE(m_lowerBound) << ")"
				<< std::endl;
		}
#endif
	}

	if (M.empty()) { // new solution found (all MAP variables assigned)

		// get the cost of the solution path
		double costMap = getMapAssignmentValue(assignment);

		// solve the SUM subproblem
		double costSum = _ao(assignment);
		++m_numSumEvals;

		// cost of the MAP solution:
		double cost = costMap * costSum;

		if (m_solutionCost < cost) {
			m_solutionCost = cost;
			m_lowerBound = m_solutionCost * m_problem->getGlobalConstant();
			double timestamp = m_timer.elapsed();
			std::cout << "["
				<< std::setw(9) << timestamp << "] u "
				<< std::setw(8) << -1 << " "
				<< std::setw(12) << m_expansions << " "
				<< std::setw(12) << m_lowerBound << " (" << -ELEM_ENCODE(m_lowerBound) << ")"
				//<< " (sum=" << costSum << ", map=" << costMap << ")"
				<< std::endl;
			std::cout.flush();
		}

		return;
	}

	// variable ordering heuristic
	std::vector<std::pair<int, double> > values;
	int nextVar = selectNextVar(M, values, assignment);
	M.erase(nextVar);
	F.erase(nextVar);
	assert(assignment[nextVar] == NONE);

	// sort values in descending order of their heuristics
	std::sort(values.begin(), values.end(), heurGreater());

	// node expansion
	for (std::vector<std::pair<int, double> >::iterator it = values.begin();
			it != values.end(); ++it) {

		int val = it->first;
		list<pair<int, int> > changes;

		// check if consistent value
		if (m_currentDomains.empty() == false &&
			m_currentDomains[nextVar][val] == false) {
			continue;
		}

		double ub = it->second;
		if (ub > m_solutionCost) { // also prune infeasible values

			++m_expansions;
			assignment[nextVar] = val;

			// propagate
			if (m_options->propagation == CP_FULL) {
				// propagate the current partial MAP assignemt
				std::vector<int> vars;
				for (std::list<int>::iterator it = m_mapVars.begin();
						it != m_mapVars.end(); ++it) {
					int x = (*it);
					int y = assignment[x];
					if (y >= 0) {
						int vid = m_var2sat[x][y];
						vars.push_back(vid);
					}
				}

				bool conflict = propagate(vars);
				if (conflict == false) {
					m_heuristic->update(assignment); // propagate the new assignment
					bbDynamic(M, assignment, F); // recursive call
				} else {
					++m_deadendsCP;
				}
			} else if (m_options->propagation == CP_UNIT) {
				// propagate the assignment
				m_zchaff.dlevel()++;
				if (lookahead(nextVar, val, changes, F)) {
					m_heuristic->update(assignment); // propagate the new assignment
					bbDynamic(M, assignment, F); // recursive call
				} else {
					++m_deadendsCP;
				}
			} else {
				m_heuristic->update(assignment); // propagate the new assignment
				bbDynamic(M, assignment, F); // recursive call
			}

			// backtrack
			if (m_options->propagation == CP_UNIT) {
				std::list<std::pair<int,int> >::iterator ci = changes.begin();
				for (; ci != changes.end(); ++ci) {
					std::pair<int,int> &p = (*ci);
					m_currentDomains[p.first][p.second] = true;
				}
				m_zchaff.release_assignments();
				m_zchaff.dlevel()--;
			}
		} else {
			++m_deadendsUB;
		}
	}
}


// recursive depth first branch and bound over the MAP variables
void BBBT::bbStatic(std::list<int> M, std::vector<int> assignment) {

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
				<< std::setw(12) << m_lowerBound << " (" << -ELEM_ENCODE(m_lowerBound) << ")"
				<< std::endl;
		}
#endif
	}

	if (M.empty()) { // new solution found

#ifdef UAI_COMPETITION
		// output the solution
		writeSolution(assignment);
#endif

		// get the cost of the solution path
		double costMap = getMapAssignmentValue(assignment);

		// solve the SUM subproblem
		double costSum = _ao(assignment);
		++m_numSumEvals;

		// cost of the MAP solution:
		double cost = costMap * costSum;

		if (m_solutionCost < cost) {
			m_solutionCost = cost;
			m_lowerBound = m_solutionCost*m_problem->getGlobalConstant();
			double timestamp = m_timer.elapsed();
			std::cout << "["
				<< std::setw(9) << timestamp << "] u "
				<< std::setw(8) << -1 << " "
				<< std::setw(12) << m_expansions << " "
				<< std::setw(12) << m_lowerBound << " (" << -ELEM_ENCODE(m_lowerBound) << ")"
				<< std::endl;
			std::cout.flush();
		}

		return;
	}

	// select next MAP variable (static ordering)
	int currVar = M.front();
	M.pop_front();
	assert(assignment[currVar] == NONE);

	// value ordering heuristic
	std::vector<double> bounds;
	m_heuristic->getHeur(currVar, bounds, assignment);
	std::vector<std::pair<int, double> > values;
	for (int val = 0; val < m_problem->getDomainSize(currVar); ++val) {
		values.push_back(std::make_pair(val, bounds[val]));
	}

	// sort values in descending order of their heuristics - largest UB first
	std::sort(values.begin(), values.end(), heurGreater());

	// node expansion
	for (std::vector<std::pair<int, double> >::iterator it = values.begin();
			it != values.end(); ++it) {

		int val = it->first;
		list<pair<int, int> > changes;

		// check if consistent value
		if (m_currentDomains.empty() == false &&
			m_currentDomains[currVar][val] == false) {
			continue;
		}


		double ub = it->second;
		if (ub > m_solutionCost) { // also prune infeasible values

			++m_expansions;
			assignment[currVar] = val;

			// propagate assignment
			if (m_options->propagation == CP_FULL) {

				// propagate the current partial assignemt
				std::vector<int> vars;
				for (std::list<int>::iterator it = m_mapVars.begin();
						it != m_mapVars.end(); ++it) {
					int x = (*it);
					int y = assignment[x];
					if (y >= 0) {
						int vid = m_var2sat[x][y];
						vars.push_back(vid);
					}
				}

				bool conflict = propagate(vars);
				if (conflict == false) {
					int nextVar = ( M.empty() ? NONE : M.front() );
					m_heuristic->updateIncr(currVar, nextVar, assignment);
					bbStatic(M, assignment); // recursive call
					m_heuristic->rollback();
				} else {
					++m_deadendsCP;
				}
			} else if (m_options->propagation == CP_UNIT) {

				// propagate the assignment
				m_zchaff.dlevel()++;
				if (lookahead(currVar, val, changes,
						m_pseudotree->getNode(currVar)->getSubprobVars())) {
					int nextVar = ( M.empty() ? NONE : M.front() );
					m_heuristic->updateIncr(currVar, nextVar, assignment);
					bbStatic(M, assignment); // recursive call
					m_heuristic->rollback();
				} else {
					++m_deadendsCP;
				}
			} else { // no propagation
				int nextVar = ( M.empty() ? NONE : M.front() );
				m_heuristic->updateIncr(currVar, nextVar, assignment);
				bbStatic(M, assignment); // recursive call
				m_heuristic->rollback();
			}

			// backtrack
			if (m_options->propagation == CP_UNIT) {

				std::list<std::pair<int,int> >::iterator ci = changes.begin();
				for (; ci != changes.end(); ++ci) {
					std::pair<int,int> &p = (*ci);
					m_currentDomains[p.first][p.second] = true;
				}
				m_zchaff.release_assignments();
				m_zchaff.dlevel()--;
			}
		} else {
			++m_deadendsUB;
		}
	}
}



// Main solver routine
int BBBT::solve() {

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
	m_cacheHits = 0;

	try {

		// init the AO search space
		ao_space.reset(new AobbSearchSpace());
		ao_propagator.reset(new BoundPropagator(m_problem.get(), ao_space.get(), m_pseudotree.get(), m_options));
		if (m_options->propagation == CP_UNIT) {
			ao_propagator->setSatSolver(&m_zchaff);
		} else {
			ao_propagator->setSatSolver(NULL);
		}
		ao_space->cache = new AobbCacheTable(m_problem->getN());

		if (m_useStaticVo) {
			// get the static order of the MAP variables and corresponding clusters
			std::list<int> M;
			m_heuristic->getStaticSearchOrder(M);

			// depth first recursive branch and bound over the MAP variables
			bbStatic(M, m_assignment);
		} else {
			// select MAP variables
			std::list<int> L;
			m_heuristic->getStaticSearchOrder(L);
			std::set<int> M(L.begin(), L.end()), F;
			for (int v = 0; v < m_problem->getN(); ++v) F.insert(v);

			// depth first recursive branch and bound over the MAP variables
			bbDynamic(M, m_assignment, F);
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
	m_solutionCost *= m_problem->getGlobalConstant();

	// output solution (if found)
	std::cout << std::endl << std::endl;
	std::cout << "----------------- Search done ------------------" << std::endl;
	std::cout << "Problem name:        " << m_problem->getName() << std::endl;
	std::cout << "Status:              " << solver_status[res] << std::endl;
	std::cout << "Nodes:               " << m_expansions << std::endl;
	std::cout << "OR nodes:            " << ao_space->nodesOR << std::endl;
	std::cout << "AND nodes:           " << ao_space->nodesAND << std::endl;
	std::cout << "Cache hits:          " << m_cacheHits << std::endl;
	std::cout << "Deadends (CP):       " << m_deadendsCP << std::endl;
	std::cout << "Deadends (UB):       " << m_deadendsUB << std::endl;
	std::cout << "SUM evaluations:     " << m_numSumEvals << std::endl;
	std::cout << "Time elapsed:        " << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
	std::cout << "Preprocessing:       " << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
	std::cout << "Solution:            " << m_solutionCost << " (" << ELEM_ENCODE(m_solutionCost) << ")" << std::endl;
	std::cout << "------------------------------------------------" << std::endl;

	// clean up
	if (EmergencyMemory) delete[] EmergencyMemory;

	return res;
}


// DONE: radu
double BBBT::getMapAssignmentValue(const std::vector<int>& assignment) {

	// Safety checks
	double cost = 1.0;
	for (std::list<int>::iterator it = m_mapVars.begin();
			it != m_mapVars.end(); ++it) {

		int v = (*it);
		const list<Function*>& funs = m_pseudotree->getFunctions(v);
		for (std::list<Function*>::const_iterator it = funs.begin();
				it != funs.end(); ++it) {
			cost *= (*it)->getValue(assignment);
		}
	}

	return cost;
}

////////////////////////////////////////////////////////////////////////////////

// DONE: radu
void BBBT::_initAoSearch() {

	// create root OR node (dummy variable)
	AobbSearchNode* node = new AobbSearchNodeOR(NULL, m_lastMapVar, 0);
	ao_space->root = node;

	// create dummy AND node (domain size 1) with global constant as label
	AobbSearchNode* next = new AobbSearchNodeAND(node, m_assignment[m_lastMapVar], 1.0);
	ao_space->root->setChild(next);

	assert(ao_stack.empty());
	ao_stack.push(next);
}

// DONE: radu
void BBBT::_heuristicOR(AobbSearchNode* n) {

	// Safety checks
	assert(n->getType() == NODE_OR);

	int v = n->getVar();
	double d;
	double* dv = new double[m_problem->getDomainSize(v) * 2];
	const list<Function*>& funs = m_pseudotree->getFunctions(v);

	for (val_t i = 0; i < m_problem->getDomainSize(v); ++i) {
		m_assignment[v] = i;

		// compute heuristic value
		dv[2 * i] = ELEM_NAN;

		// precompute label value
		d = ELEM_ONE;
		for (std::list<Function*>::const_iterator it = funs.begin();
				it != funs.end(); ++it) {
			d *= (*it)->getValue(m_assignment);
		}

		// store label and heuristic into cache table
		dv[2 * i + 1] = d; // label
		dv[2 * i] *= d; // heuristic (includes label)

	}

	n->setHeurCache(dv);

} // AOBB::heuristicOR

// DONE: radu
bool BBBT::_doProcess(AobbSearchNode* node) {
	assert(node);
	if (node->getType() == NODE_AND) {
		int var = node->getVar();
		int val = node->getVal();
		m_assignment[var] = val; // record assignment

	} else { // NODE_OR
		// do nothing
	}

	return false; // default
}

// DONE: radu
void BBBT::_addCacheContext(AobbSearchNode* node, const set<int>& ctxt) const {

	context_t sig;
	sig.reserve(ctxt.size());
	for (set<int>::const_iterator itC = ctxt.begin(); itC != ctxt.end(); ++itC) {
		sig.push_back(m_assignment[*itC]);
	}

	node->setCacheContext(sig);
}

// DONE: radu
bool BBBT::_doCaching(AobbSearchNode* node) {

	assert(node);
	int var = node->getVar();
	PseudotreeNode* ptnode = m_pseudotree->getNode(var);

	if (node->getType() == NODE_AND) { // AND node -> reset associated adaptive cache tables

		// do nothing at AND nodes
	    const list<int>& resetList = ptnode->getCacheReset();
	    for (list<int>::const_iterator it=resetList.begin(); it!=resetList.end(); ++it)
	    	ao_space->cache->reset(*it);

	} else { // OR node, try actual caching

		if (!ptnode->getParent())
			return false;

		if (ptnode->getFullContext().size() <= ptnode->getParent()->getFullContext().size()) {
			// add cache context information
			//addCacheContext(node, ptnode->getFullContext());
			_addCacheContext(node, ptnode->getCacheContext());

			// try to get value from cache
			try {
				// will throw int(UNKNOWN) if not found
#ifndef NO_ASSIGNMENT
				pair<double,vector<val_t> > entry = ao_space->cache->read(var, node->getCacheContext());
				node->setValue( entry.first ); // set value
				node->setOptAssig( entry.second ); // set assignment
#else
				double entry = ao_space->cache->read(var, node->getCacheContext());
				node->setValue(entry); // set value
#endif
				node->setLeaf(); // mark as leaf
				++m_cacheHits;

				return true;
			} catch (...) { // cache lookup failed
				node->setCachable(); // mark for caching later
			}
		}
	} // if on node type

	return false; // default, no caching applied

} // BBBT::doCaching

// DONE: radu
bool BBBT::_doPruning(AobbSearchNode* node) {

	assert(node);

	if (_canBePruned(node)) {
		node->setLeaf();
		node->setPruned();
		node->setValue(ELEM_ZERO); // dead end
		++m_deadendsCP;

		return true;
	}

	return false; // default false
} // AOBB::doPruning

bool BBBT::_canBePruned(AobbSearchNode* n) {

	// propagate current assignment (AND nodes only)
	if (n->getType() == NODE_AND) {

		if (m_options->propagation == CP_FULL) {

			// propagate the assignment (as unit propagation)
			std::vector<int> vars = m_assumptions;
			getPathAssignment(n, vars);

			return propagate(vars); // 'true' if conflict

		} else if (m_options->propagation == CP_UNIT) {

			// Propagate the assignment
			m_zchaff.dlevel()++;
			if (!lookahead(n->getVar(), n->getVal(), n->changes(),
					m_pseudotree->getNode(n->getVar())->getSubprobVars())) {
				// found inconsistency
				++m_deadendsCP;
				return true; // 'true' if conflict
			}
		}
	}

	// default, no pruning possible
	return false;

} // BBBT::canBePruned

void BBBT::getPathAssignment(AobbSearchNode* node, std::vector<int>& vars) {

	assert(node->getType() == NODE_AND);
	AobbSearchNode *crt = node, *par = node->getParent();
	while (true) {
		int vvar = crt->getVar();
		int vval = crt->getVal();
		int vid = m_var2sat[vvar][vval];
		vars.push_back(vid);

		crt = par->getParent(); // AND
		if (crt == NULL) break;
		par = crt->getParent(); // OR
	}

}

// DONE: radu
AobbSearchNode* BBBT::_nextLeaf() {

	AobbSearchNode* node = this->_nextNode();
	while (node) {

		// check for time limit violation
		if (m_timer.timeout()) {
			throw SEARCH_TIMEOUT;
		}

		if (_doProcess(node)) { // initial processing
			return node;
		}
		if (_doCaching(node)) { // caching?
			return node;
		}
		if (_doPruning(node)) { // pruning?
			return node;
		}
		if (_doExpand(node)) { // node expansion
			return node;
		}
		node = this->_nextNode();
	}
	return NULL;
}

// DONE: radu
bool BBBT::_generateChildrenAND(AobbSearchNode* n, vector<AobbSearchNode*>& chi) {

	assert(n && n->getType() == NODE_AND);
	//assert(n->getChildren() == NULL);

	if (n->getChildren()) {  // node previously expanded
		if (n->getChildCountAct() == 0) {
			n->clearChildren(); // previously expanded, but all children deleted
		} else {
			for (size_t i = 0; i < n->getChildCountFull(); ++i) {
				if (n->getChildren()[i])
					chi.push_back(n->getChildren()[i]);
			}
			return false;
		}
	}

	int var = n->getVar();
	PseudotreeNode* ptnode = m_pseudotree->getNode(var);
	int depth = ptnode->getDepth();

	// increase AND node expansions
	ao_space->nodesAND += 1;

	// create new OR children (going in reverse due to reversal on stack)
	for (vector<PseudotreeNode*>::const_iterator it =
			ptnode->getChildren().begin(); it != ptnode->getChildren().end();
			++it) {
		int vChild = (*it)->getVar();
		AobbSearchNodeOR* c = new AobbSearchNodeOR(n, vChild, depth + 1);
		// Compute and set heuristic estimate, includes child labels
		_heuristicOR(c);
		chi.push_back(c);

	} // for loop over new OR children

	if (chi.empty()) {
		n->setLeaf(); // terminal node
		n->setValue(ELEM_ONE);
		return true; // no children
	}

	// (use reverse iterator due to stack reversal)
	n->addChildren(chi);

	return false; // default

} // AO::generateChildrenAND

// DONE: radu
bool BBBT::_generateChildrenOR(AobbSearchNode* n, vector<AobbSearchNode*>& chi) {

	assert(n && n->getType() == NODE_OR);
	//assert(n->getChildren() == NULL);

	if (n->getChildren()) {  // node previously expanded
		if (n->getChildCountAct() == 0) {
			n->clearChildren(); // previously expanded, but all children deleted
		} else {
			for (size_t i = 0; i < n->getChildCountFull(); ++i) {
				if (n->getChildren()[i])
					chi.push_back(n->getChildren()[i]);
			}
			return false;
		}
	}


	int var = n->getVar();

	// increase OR node expansions
	ao_space->nodesOR += 1;

	// retrieve precomputed labels
	double* heur = n->getHeurCache();
	for (val_t i = 0 ; i < m_problem->getDomainSize(var); ++i) {

		// consistent value
		if (m_currentDomains.empty() == false &&
			m_currentDomains[var][i] == false) {
			continue;
		}

		double label = heur[2 * i + 1];
		AobbSearchNodeAND* c = new AobbSearchNodeAND(n, i, label); // uses cached label
		// set cached heur. value (includes the weight)
		c->setHeur(heur[2 * i]);
		chi.push_back(c);
	}

	if (chi.empty()) { // deadend
		n->setLeaf();
		n->setValue(ELEM_ZERO);
		return true; // no children
	}

	// (use reverse iterator due to stack reversal)
	n->addChildren(chi);

	return false; // default

} // AO::generateChildrenOR

// DONE: radu
AobbSearchNode* BBBT::_nextNode() {
	AobbSearchNode* n = NULL;
	if (!n && ao_stack.size()) {
		n = ao_stack.top();
		ao_stack.pop();
	}
	return n;
}

// DONE: radu
bool BBBT::_doExpand(AobbSearchNode* n) {
	assert(n);
	ao_expand.clear();

	if (n->getType() == NODE_AND) {  // AND node

		if (_generateChildrenAND(n, ao_expand))
			return true; // no children
		for (vector<AobbSearchNode*>::reverse_iterator it = ao_expand.rbegin();
				it != ao_expand.rend(); ++it)
			ao_stack.push(*it);

	} else {  // OR node

		if (_generateChildrenOR(n, ao_expand))
			return true; // no children
		for (vector<AobbSearchNode*>::reverse_iterator it = ao_expand.rbegin();
				it != ao_expand.rend(); ++it) {
			ao_stack.push(*it);
		} // for loop

	} // if over node type

	return false; // default false
}

// solve the problem
double BBBT::_ao(const std::vector<int>& assignment) {

	// a temporary buffer
	double result = ELEM_NAN;

	// reinit the search space
	ao_space->root = NULL;
	m_assignment = assignment;

	// initialize the assumptions corresponding to the MAP assignment
	if (m_options->propagation == CP_FULL) {
		m_assumptions.clear();
		for (std::list<int>::iterator it = m_mapVars.begin();
				it != m_mapVars.end(); ++it) {
			int var = (*it);
			int val = assignment[var];
			int vid = m_var2sat[var][val];
			m_assumptions.push_back(vid);
		}
	}

	// init root of the search space
	_initAoSearch();

	AobbSearchNode* n = _nextLeaf();
	while (n) { // throws timeout
		ao_propagator->propagate(n, m_timer, m_currentDomains, false); // true = report solutions
		n = _nextLeaf();
	}

	// Conditioned SUM problem solved
	result = ao_space->root->getValue();

	return result;
}

#ifdef UAI_COMPETITION
void BBBT::writeSolution(std::vector<int>& sol) {

	// reconstruct full solution (original variables)
	int nOrg = m_problem->getNOrg();
	vector<val_t> curSolution;
	curSolution.resize(nOrg, UNKNOWN);

	const map<int, int>& old2new = m_problem->getOld2New();
	const map<int, int>& evidence = m_problem->getEvidence();
	for (int i = 0; i < nOrg; ++i) {
		map<int, int>::const_iterator itRen = old2new.find(i);
		if (itRen != old2new.end()) { // var part of solution
			curSolution.at(i) = sol.at(itRen->second);
		} else {
			map<int, val_t>::const_iterator itEvid = evidence.find(i);
			if (itEvid != evidence.end())  // var part of evidence
				curSolution.at(i) = itEvid->second;
			else
				// var had unary domain
				curSolution.at(i) = 0;
		}
	}

	// save the assignment to the output file
	std::string output = m_options->outputFile;
	if (output.empty() == false) {
		std::ofstream out;
		out.open(output.c_str(), std::ofstream::out | std::ofstream::app);
		if (out.fail()) {
			std::cerr << "Cannot append to output file " << output << std::endl;
			exit(EXIT_FAILURE);
		}

		out << "-BEGIN-" << std::endl;
		const std::set<int>& query = m_problem->getQuery();
		std::set<int>::const_iterator li = query.begin();
		out << query.size();
		for (; li != query.end(); ++li) {
			out << " " << (*li) << " " << curSolution[*li];
		}
		out << std::endl;
		out.close();
	}
}
#endif

