/*
 * astar.cpp
 *
 *  Created on: Mar 31, 2014
 *      Author: radu
 */


#include "astar.h"
#include "mini_cluster_tree_heur.h"
#include "mini_cluster_tree_incr_heur.h"



// initialize
void Astar::init() {

	// initialize
	Search::init();

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

	// init heuristic generator
	std::cout << "Computing the mini-cluster-tree heuristic..." << std::endl;
	m_heuristic.reset(new MiniClusterTreeHeur(m_problem.get(),
				m_pseudotree.get(), m_options));

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
	std::copy(mapSearchOrder.begin(), mapSearchOrder.end(),
			std::inserter(m_mapVars, m_mapVars.begin()));
	m_lastMapVar = mapSearchOrder.back();
	m_firstMapVar = mapSearchOrder.front();
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

void Astar::getMapAssignment(const AstarSearchNode& n,
		std::vector<int>& assignment) {
	assignment.resize(m_problem->getN(), UNKNOWN);
	for (list<pair<int, int> >::const_iterator li = n.assignment().begin();
			li != n.assignment().end(); ++li) {
		int var = li->first;
		int val = li->second;
		assignment[var] = val;
	}
}

// A* search over the MAP variables
void Astar::astar() {

	// initialize the search
	AstarSearchNode root(UNKNOWN, UNKNOWN);
	root.f() = 0;
	m_open.push(root);

	while ( !m_open.empty() ) {

		// check for timeout
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
					<< std::setw(12) << m_expansions << " (opened " << m_open.size() << ")"
					<< std::endl;
			}
#endif
		}

		AstarSearchNode n = m_open.top();
		m_open.pop();

		if (n.f() >= m_solutionCost) {
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

			// cost of the MAP solution:
			double cost = costMap + costSum;

			//if (cost < m_solutionCost)
			m_solutionCost = cost;
			double timestamp = m_timer.elapsed();
			std::cout << "["
				<< std::setw(9) << timestamp << "] u "
				<< std::setw(8) << -1 << " "
				<< std::setw(12) << m_expansions << " "
				<< std::setw(12) << m_solutionCost
				<< " (sum=" << costSum << ", map=" << costMap << ")"
				<< std::endl;
			std::cout.flush();
			++m_sumEvals;

		} else {
			expand(n);
			++m_expansions;
		}
	}
}

/*
// select next unassigned MAP variable
int Astar::selectNextVar(const set<int>& M,
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
*/

// select next unassigned MAP variable (A* uses lexicographic order)
int Astar::selectNextVar(int currVar,
		std::vector<std::pair<int, double> >& values,
		const std::vector<int>& assignment) {

	values.clear();

	// select next uninstantiated MAP variable
	int nextVar = UNKNOWN;
	if (currVar == UNKNOWN) {
		nextVar = m_firstMapVar;
	} else {
		PseudotreeNode* ptNode = m_pseudotree->getNode(currVar);
		const std::vector<PseudotreeNode*>& children = ptNode->getChildren();
		assert(children.size() == 1);
		nextVar = children[0]->getVar();
	}

	// safety checks
	assert(nextVar != NONE);

	// compute the lower bounds associated with the domain values
	std::vector<double> bounds;
	m_heuristic->getHeur(nextVar, bounds, assignment);
	for (int val = 0; val < m_problem->getDomainSize(nextVar); ++val) {
		values.push_back(std::make_pair(val, bounds[val]));
	}

	return nextVar;
}

int Astar::selectNextVar(const set<int>& M,
		std::vector<std::pair<int, double> >& values,
		const std::vector<int>& assignment) {

	values.clear();
	int nextVar = *(M.begin());
	std::vector<double> bounds;
	m_heuristic->getHeur(nextVar, bounds, assignment);

	for (int val = 0; val < m_problem->getDomainSize(nextVar); ++val) {
		values.push_back(std::make_pair(val, bounds[val]));
	}

	// safety checks
	assert(nextVar != NONE);
	return nextVar;
}

// expand the search node
void Astar::expand(const AstarSearchNode& n) {

	// get the current partial MAP assignment
	std::vector<int> assignment;
	getMapAssignment(n, assignment);
	m_heuristic->update(assignment); // propagate the new assignment

	// get the current instantiated MAP variables
	std::set<int> currMapVars;
	for (list<pair<int, int> >::const_iterator li = n.assignment().begin();
			li != n.assignment().end(); ++li) {
		currMapVars.insert(li->first);
	}

	// get the future uninstantiated MAP variables
	std::set<int> M;
	std::set_difference(m_mapVars.begin(), m_mapVars.end(),
			currMapVars.begin(), currMapVars.end(),
			std::inserter(M, M.begin()));

	// select the most promissing unassigned MAP variable
	std::vector<std::pair<int, double> > values;
//	int nextVar = selectNextVar(n.var(), values, assignment);
//	bool term = (nextVar == m_lastMapVar);
	int nextVar = selectNextVar(M, values, assignment);
	bool term = (M.size() == 1);
	for (std::vector<std::pair<int, double> >::iterator vi = values.begin();
			vi != values.end(); ++vi) {
		int val = vi->first;
		double upbo = vi->second;

		AstarSearchNode succ(nextVar, val);
//		succ.f() = -ELEM_ENCODE(upbo);
		assignment[nextVar] = val;
		double lb = -ELEM_ENCODE( upbo );
		double co = -ELEM_ENCODE( getPartialMapAssignmentValue(assignment) );
		succ.g() = co;
		succ.h() = lb - succ.g();
		succ.f() = succ.g() + m_epsilon*succ.h();

		succ.terminal() = term;
		succ.depth() = n.depth() + 1;
		succ.assignment() = n.assignment();
		succ.assignment().push_back(std::make_pair(nextVar, val));

		m_open.push(succ);
	}
}

// Main solver routine
int Astar::solve() {

	// prologue
	std::cout << "--- Starting search ---" << std::endl;
	int res = SEARCH_SUCCESS;
	char* EmergencyMemory = new char[16384]; // in order to handle bad alloc

	// init search (timer, search space, assignment)
	m_tmSolve = 0;
	m_timer.reset(m_options->timeLimit);
	m_timer.start();
	m_solutionCost = INFINITY; // holds the current best MAP
	m_assignment.resize(m_problem->getN(), NONE);
	m_expansions = 0;

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

		// A* search over the MAP variables
		astar();

		// Proved optimality
		m_solved = true;

	} catch (std::bad_alloc& ba) {
		delete[] EmergencyMemory;
		EmergencyMemory = NULL;
		std::cerr << "OUT OF MEMORY";
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
	m_solutionCost += -ELEM_ENCODE( m_problem->getGlobalConstant() );

	// output solution (if found)
	std::cout << std::endl << std::endl;
	std::cout << "----------------- Search done ------------------" << std::endl;
	std::cout << "Problem name:        " << m_problem->getName() << std::endl;
	std::cout << "Status:              " << solver_status[res] << std::endl;
	std::cout << "Nodes:               " << m_expansions << std::endl;
	std::cout << "OR nodes:            " << ao_space->nodesOR << std::endl;
	std::cout << "AND nodes:           " << ao_space->nodesAND << std::endl;
	std::cout << "SUM evaluations:     " << m_sumEvals << std::endl;
	std::cout << "Time elapsed:        " << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
	std::cout << "Preprocessing:       " << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
//	std::cout << "Solution:            " << m_solutionCost << " (" << 1/ELEM_DECODE(m_solutionCost) << ")" << std::endl;
	std::cout << "Solution:            " << 1/ELEM_DECODE(m_solutionCost) << " (" << -m_solutionCost  << ")" << std::endl;
	std::cout << "------------------------------------------------" << std::endl;

	// clean up
	if (EmergencyMemory) delete[] EmergencyMemory;

	return res;
}

////////////////////////////////////////////////////////////////////////////////

// DONE: radu
void Astar::_initAoSearch() {

	assert( m_lastMapVar != UNKNOWN );
	assert( m_assignment[m_lastMapVar] != UNKNOWN );

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
void Astar::_heuristicOR(AobbSearchNode* n) {

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
double Astar::getMapAssignmentValue(const std::vector<int>& assignment) {

	// Safety checks
	double cost = 1.0;
	for (std::set<int>::iterator it = m_mapVars.begin();
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

// DONE: radu
double Astar::getPartialMapAssignmentValue(const std::vector<int>& assignment) {

	// Safety checks
	double cost = 1.0;
	for (std::set<int>::iterator it = m_mapVars.begin();
			it != m_mapVars.end(); ++it) {

		int v = (*it);
		const list<Function*>& funs = m_pseudotree->getFunctions(v);
		for (std::list<Function*>::const_iterator it = funs.begin();
				it != funs.end(); ++it) {
			Function* f = (*it);
			if (f->isInstantiated(assignment)) {
				cost *= (*it)->getValue(assignment);
			}
		}
	}

	return cost;
}


// DONE: radu
bool Astar::_doProcess(AobbSearchNode* node) {
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
void Astar::_addCacheContext(AobbSearchNode* node, const set<int>& ctxt) const {

	context_t sig;
	sig.reserve(ctxt.size());
	for (set<int>::const_iterator itC = ctxt.begin(); itC != ctxt.end(); ++itC) {
		sig.push_back(m_assignment[*itC]);
	}

	node->setCacheContext(sig);
}

// DONE: radu
bool Astar::_doCaching(AobbSearchNode* node) {

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

#ifndef NO_ASSIGNMENT
				pair<double,vector<val_t> > entry = ao_space->cache->read(var, node->getCacheContext());
				node->setValue( entry.first ); // set value
				node->setOptAssig( entry.second ); // set assignment
#else
				double entry = ao_space->cache->read(var, node->getCacheContext());
				node->setValue(entry); // set value
#endif

				// will throw int(UNKNOWN) if not found
//				double entry = ao_space->cache->read(var, node->getCacheContext());
//				node->setValue(entry); // set value
				node->setLeaf(); // mark as leaf

				return true;
			} catch (...) { // cache lookup failed
				node->setCachable(); // mark for caching later
			}
		}
	} // if on node type

	return false; // default, no caching applied

} // BBBT::doCaching

// DONE: radu
bool Astar::_doPruning(AobbSearchNode* node) {

	assert(node);

	if (_canBePruned(node)) {
		node->setLeaf();
		node->setPruned();
		node->setValue(ELEM_ZERO); // dead end

		return true;
	}

	return false; // default false
} // AOBB::doPruning

bool Astar::_canBePruned(AobbSearchNode* n) {

	// propagate current assignment (AND nodes only)
	if (n->getType() == NODE_AND) {

//		if (m_options->propagation == CP_FULL) {
//
//			// propagate the assignment (as unit propagation)
//			std::vector<int> vars = m_assumptions;
//			getPathAssignment(n, vars);
//
//			return propagate(vars); // 'true' if conflict
//
//		} else if (m_options->propagation == CP_UNIT) {
//
//			// Propagate the assignment
//			m_zchaff.dlevel()++;
//			if (!lookahead(n->getVar(), n->getVal(), n->changes(),
//					m_pseudotree->getNode(n->getVar())->getSubprobVars())) {
//				// found inconsistency
//				return true; // 'true' if conflict
//			}
//		}
	}

	// default, no pruning possible
	return false;

} // BBBT::canBePruned

// DONE: radu
AobbSearchNode* Astar::_nextLeaf() {

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
bool Astar::_generateChildrenAND(AobbSearchNode* n, vector<AobbSearchNode*>& chi) {

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
bool Astar::_generateChildrenOR(AobbSearchNode* n, vector<AobbSearchNode*>& chi) {

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
AobbSearchNode* Astar::_nextNode() {
	AobbSearchNode* n = NULL;
	if (!n && ao_stack.size()) {
		n = ao_stack.top();
		ao_stack.pop();
	}
	return n;
}

// DONE: radu
bool Astar::_doExpand(AobbSearchNode* n) {
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
double Astar::_ao(const std::vector<int>& assignment) {

	// a temporary buffer
	double result = ELEM_NAN;

	// reinit the search space
	ao_space->root = NULL;
	m_assignment = assignment;

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


