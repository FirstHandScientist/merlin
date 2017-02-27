/*
 * ao.cpp
 *
 *  Created on: Sep 6, 2013
 *      Author: radu
 */

#include "ao.h"

//#define DEBUG

// initialize
void AO::init() {

	// init problem instance
	Search::init();

	Timer tm;
	tm.start();

	// Encode determinism
	encodeDeterminism();

	tm.stop();

	// heuristic might have changed problem functions, pseudotree needs remapping
	m_pseudotree->resetFunctionInfo(m_problem->getFunctions());

	// Record time after initialization
	cout << "Initialization complete: " << (m_tmLoad + m_tmHeuristic)
			<< " seconds" << std::endl;
}


void AO::getPathAssignment(AobbSearchNode* node, std::vector<int>& vars) {

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
AobbSearchNode* AO::initSearch() {
	assert(m_space->root == NULL);

	// Add initial set of dummy nodes.

	// create root OR node (dummy variable)
	PseudotreeNode* ptroot = m_pseudotree->getRoot();
	AobbSearchNode* node = new AobbSearchNodeOR(NULL, ptroot->getVar(), -1);
	m_space->root = node;
	double h = heuristicOR(node);
	node->setHeur(h);

	// create dummy AND node (domain size 1) with global constant as label
	AobbSearchNode* next = new AobbSearchNodeAND(node, 0, m_problem->getGlobalConstant());
	m_space->root->setChild(next);
	next->setHeur(node->getHeurCache()[0]);

	return next;
}

// DONE: radu
double AO::heuristicOR(AobbSearchNode* n) {

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

	return ELEM_ZERO;

} // AOBB::heuristicOR


// DONE: radu
bool AO::doProcess(AobbSearchNode* node) {
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
void AO::addCacheContext(AobbSearchNode* node, const set<int>& ctxt) const {

	context_t sig;
	sig.reserve(ctxt.size());
	for (set<int>::const_iterator itC = ctxt.begin(); itC != ctxt.end(); ++itC) {
		sig.push_back(m_assignment[*itC]);
	}

	node->setCacheContext(sig);
}

// DONE: radu
bool AO::doCaching(AobbSearchNode* node) {

	assert(node);
	int var = node->getVar();
	PseudotreeNode* ptnode = m_pseudotree->getNode(var);

	if (node->getType() == NODE_AND) { // AND node -> reset associated adaptive cache tables

		// do nothing at AND nodes
	    const list<int>& resetList = ptnode->getCacheReset();
	    for (list<int>::const_iterator it=resetList.begin(); it!=resetList.end(); ++it)
	    	m_space->cache->reset(*it);

	} else { // OR node, try actual caching

		if (!ptnode->getParent())
			return false;

		if (ptnode->getFullContext().size() <= ptnode->getParent()->getFullContext().size()) {
			// add cache context information
			//addCacheContext(node, ptnode->getFullContext());
			addCacheContext(node, ptnode->getCacheContext());

			// try to get value from cache
			try {
				// will throw int(UNKNOWN) if not found
#ifndef NO_ASSIGNMENT
				pair<double,vector<val_t> > entry = m_space->cache->read(var, node->getCacheContext());
				node->setValue( entry.first ); // set value
				node->setOptAssig( entry.second ); // set assignment
#else
				double entry = m_space->cache->read(var, node->getCacheContext());
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

} // AOBB::doCaching

// DONE: radu
AobbSearchNode* AO::nextLeaf() {

	AobbSearchNode* node = this->nextNode();
	while (node) {

		// check for time limit violation
		if (m_timer.timeout()) {
			throw SEARCH_TIMEOUT;
		}

		// output intermediate results if verbose mode is enabled
		if (m_options->verbose) {
			// logging
#ifdef SEARCH_CHECKPOINT_NODE
			size_t expansions = m_space->nodesAND;
			if (expansions % SEARCH_CHECKPOINT_NODE == 0 &&
				expansions > m_prevNodeCheckpoint) {

				double timestamp = m_timer.elapsed();
				m_prevNodeCheckpoint = expansions;
				std::cout << "[*"
					<< std::setw(8) << timestamp << "] u "
					<< std::setw(8) << -1 << " "
					<< std::setw(12) << m_space->nodesOR << " "
					<< std::setw(12) << m_space->nodesAND << " "
					<< std::setw(12) << m_lowerBound
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
					<< std::setw(12) << m_space->nodesOR << " "
					<< std::setw(12) << m_space->nodesAND << " "
					<< std::setw(12) << m_lowerBound
					<< std::endl;
			}
#endif
		}


		if (doProcess(node)) { // initial processing
			return node;
		}
		if (doCaching(node)) { // caching?
			return node;
		}
		if (doPruning(node)) { // pruning?
			return node;
		}
		if (doExpand(node)) { // node expansion
			return node;
		}
		node = this->nextNode();
	}
	return NULL;
}

// DONE: radu
bool AO::doPruning(AobbSearchNode* node) {

	assert(node);

	if (canBePruned(node)) {
		node->setLeaf();
		node->setPruned();
		node->setValue(ELEM_ZERO); // dead end
		++m_deadendsCP;

		return true;
	}

	return false; // default false
} // AO::doPruning

// DONE: radu
bool AO::canBePruned(AobbSearchNode* n) {

	// propagate current assignment (AND nodes only)
	if (n->getType() == NODE_AND) {

		if (m_options->propagation == CP_FULL) {

			// propagate the assignment (as unit propagation)
			std::vector<int> vars;
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

} // AO::canBePruned


// DONE: radu
bool AO::generateChildrenAND(AobbSearchNode* n, vector<AobbSearchNode*>& chi) {

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
	m_space->nodesAND += 1;
	if (m_problem->isMap(var))
		m_space->nodesANDmap += 1;

	// create new OR children (going in reverse due to reversal on stack)
	for (vector<PseudotreeNode*>::const_iterator it =
			ptnode->getChildren().begin(); it != ptnode->getChildren().end();
			++it) {
		int vChild = (*it)->getVar();
		AobbSearchNodeOR* c = new AobbSearchNodeOR(n, vChild, depth + 1);
		// Compute and set heuristic estimate, includes child labels
		heuristicOR(c);
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
bool AO::generateChildrenOR(AobbSearchNode* n, vector<AobbSearchNode*>& chi) {

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
	m_space->nodesOR += 1;
	if (m_problem->isMap(var))
		m_space->nodesORmap += 1;

	// retrieve precomputed labels
	double* heur = n->getHeurCache();
	for (val_t i = 0 ; i < m_problem->getDomainSize(var); ++i) {
		if (m_currentDomains.empty() == false &&
			m_currentDomains[var][i] == false) {
			continue; // skip inconsistent values
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
AobbSearchNode* AO::nextNode() {
	AobbSearchNode* n = NULL;
	if (!n && m_stack.size()) {
		n = m_stack.top();
		m_stack.pop();
	}
	return n;
}

// DONE: radu
bool AO::doExpand(AobbSearchNode* n) {
	assert(n);
	m_expand.clear();

	if (n->getType() == NODE_AND) {  // AND node

		if (generateChildrenAND(n, m_expand))
			return true; // no children
		for (vector<AobbSearchNode*>::reverse_iterator it = m_expand.rbegin();
				it != m_expand.rend(); ++it)
			m_stack.push(*it);

	} else {  // OR node

		if (generateChildrenOR(n, m_expand))
			return true; // no children
		for (vector<AobbSearchNode*>::reverse_iterator it = m_expand.rbegin();
				it != m_expand.rend(); ++it) {
			m_stack.push(*it);
		} // for loop

	} // if over node type

	return false; // default false
}

// solve the problem
int AO::solve() {

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
		AobbSearchNode* first = initSearch();
		if (first) {
			m_stack.push(first);
		}

		AobbSearchNode* n = nextLeaf();
		while (n) { // throws timeout
#ifdef DEBUG
			cout << "Propagate leaf: " << n->toString() << endl;
#endif
			m_propagator->propagate(n, m_timer, m_currentDomains, true); // true = report solutions
			m_lowerBound = m_propagator->getLowerBound();
			n = nextLeaf();
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
	m_solutionCost = m_propagator->getLowerBound();

	// output solution (if found)
	std::cout <<  std::endl << std::endl;
	std::cout << "----------------- Search done ------------------" << std::endl;
	std::cout << "Problem name:\t\t" << m_problem->getName() << std::endl;
	std::cout << "Status:\t\t\t" << solver_status[res] << std::endl;
	std::cout << "OR nodes:\t\t" << m_space->nodesOR << std::endl;
	std::cout << "AND nodes:\t\t" << m_space->nodesAND << std::endl;
	std::cout << "OR nodes (MAP):\t\t" << m_space->nodesORmap << std::endl;
	std::cout << "AND nodes (MAP):\t" << m_space->nodesANDmap << std::endl;
	std::cout << "Cache hits:\t\t" << m_cacheHits << std::endl;
	std::cout << "Time elapsed:\t\t" << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
	std::cout << "Preprocessing:\t\t" << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
	std::cout << "Solution:\t\t" << m_solutionCost << " (" << ELEM_ENCODE(m_solutionCost) << ")" << std::endl;
	std::cout << "------------------------------------------------" << std::endl;

	// clean up
	if (EmergencyMemory) delete[] EmergencyMemory;

	return res;
}


