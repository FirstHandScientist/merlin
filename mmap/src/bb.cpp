/*
 * bb.cpp
 *
 *  Created on: Jun 7, 2013
 *      Author: radu
 */


#include "bb.h"
#include "mini_cluster_tree_heur.h"
#include "mini_cluster_tree_incr_heur.h"

// initialize
void BB::init() {

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
	if (m_options->incremental == false) {
		std::cout << "Computing the mini-cluster-tree heuristic..." << std::endl;
		m_heuristic.reset(new MiniClusterTreeHeur(m_problem.get(),
				m_pseudotree.get(), m_options));
	}  else {
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

	// heuristic might have changed problem functions, pseudotree needs remapping
	m_pseudotree->resetFunctionInfo(m_problem->getFunctions());

	// Record time after initialization
	cout << "Initialization complete: " << (m_tmLoad + m_tmHeuristic)
			<< " seconds" << std::endl;

	// MAP variables
	m_nextMapVar.resize(m_problem->getN(), NONE);
	std::vector<int> vars(mapSearchOrder.begin(), mapSearchOrder.end());
	for (int i = 0; i < (int)vars.size() - 1; ++i) {
		int currVar = vars[i];
		int nextVar = vars[i+1];
		m_nextMapVar[currVar] = nextVar;
	}

	// check for upper bound computation only
	if (m_options->upperBoundOnly) {
		std::cout << "Done." << std::endl;
		exit(EXIT_SUCCESS);
	}
}


// DONE: radu
AobbSearchNode* BB::initSearchSpace() {
	assert(m_space->root == NULL);

	// Create root OR node (root of the pseudo tree)
	PseudotreeNode* ptroot = m_pseudotree->getRoot();
	AobbSearchNode* root = new AobbSearchNodeOR(NULL, ptroot->getVar(), -1);
	m_space->root = root;
	double h = heuristicOR(root);
	root->setHeur(h);
	root->setDepth(1);

	return root;
}

// DONE: radu
double BB::heuristicOR(AobbSearchNode* n) {

	// Safety checks
	assert(n->getType() == NODE_OR);

	// check if MAP variable
	bool mapVar = m_problem->isMap(n->getVar());

	int var = n->getVar();
	double* dv = new double[m_problem->getDomainSize(var) * 2];
	const list<Function*>& funs = m_pseudotree->getFunctions(var);

	// compute the bounds for each domain value
	std::vector<double> bounds;
	std::vector<int> assignment;
	assignment.resize(m_assignment.size(), NONE);
	n->getPathAssignment(assignment);
	double h = (mapVar) ? -INFINITY : ELEM_ZERO; // the new OR nodes h value
	if (mapVar) { // MAP variable
		m_heuristic->getHeur(var, bounds, assignment);
	} else { // SUM variable
		int pvar = n->getParent()->getVar();
		if (m_problem->isMap(pvar)) {
			h = n->getParent()->getHeur();
		}
	}

	for (val_t i = 0; i < m_problem->getDomainSize(var); ++i) {

		assignment[var] = i;

		// precompute label value
		double d = ELEM_ONE;
		for (std::list<Function*>::const_iterator it = funs.begin();
				it != funs.end(); ++it) {
			d *= (*it)->getValue(assignment);
		}

		if (mapVar) {
			//double g = pathCost * d;
			//dv[2 * i] = ( g == 0.0) ? ELEM_ZERO : (bounds[i] / g );
			dv[2 * i] = bounds[i];
		} else {
			dv[2 * i] = INFINITY; // not used
		}

		// store label and heuristic into cache table
		dv[2 * i + 1] = d; // label
		//dv[2 * i] *= d; // heuristic (includes label)

		if (mapVar) {
			if (dv[2 * i] > h) {
				h = dv[2 * i]; // keep max. for OR node heuristic (MAP var)
			}
		}

	}

	n->setHeur(h); // OR node
	n->setHeurCache(dv);

	return h;

} // AOBB::heuristicOR


// DONE: radu
bool BB::doProcess(AobbSearchNode* node) {
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
void BB::addCacheContext(AobbSearchNode* node, const set<int>& ctxt) const {

	context_t sig;
	sig.reserve(ctxt.size());
	for (set<int>::const_iterator itC = ctxt.begin(); itC != ctxt.end(); ++itC) {
		sig.push_back(m_assignment[*itC]);
	}

	node->setCacheContext(sig);
}

// DONE: radu
bool BB::doCaching(AobbSearchNode* node) {

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

#ifndef NO_ASSIGNMENT
				pair<double,vector<val_t> > entry = m_space->cache->read(var, node->getCacheContext());
				node->setValue( entry.first ); // set value
				node->setOptAssig( entry.second ); // set assignment
#else
				double entry = m_space->cache->read(var, node->getCacheContext());
				node->setValue(entry); // set value
#endif

				// will throw int(UNKNOWN) if not found
//				double entry = m_space->cache->read(var, node->getCacheContext());
//				node->setValue(entry); // set value
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
bool BB::doPruning(AobbSearchNode* node) {

	assert(node);

	if (canBePruned(node)) {
		node->setLeaf();
		node->setPruned();
		if (node->getType() == NODE_OR) {
			if (ISNAN(node->getValue())) // value could be set by LDS
				node->setValue(ELEM_ZERO);
		}
		return true;
	}

	return false; // default false
} // AOBB::doPruning

// DONE: radu
AobbSearchNode* BB::nextLeaf() {

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
					<< std::setw(12) << m_space->nodesORmap << " "
					<< std::setw(12) << m_space->nodesANDmap << " "
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
					<< std::setw(12) << m_space->nodesORmap << " "
					<< std::setw(12) << m_space->nodesANDmap << " "
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

bool BB::canBePruned(AobbSearchNode* n) const {

	// do not prune SUM nodes
	if (m_problem->isMap(n->getVar()) == false ||
		n->getType() == NODE_OR) {
		return false;
	}

	if (n->getHeur() <= m_lowerBound) {
		return true;
	} else {
		return false;
	}
} // AOBB::canBePruned

// DONE: radu
bool BB::generateChildrenAND(AobbSearchNode* n, vector<AobbSearchNode*>& chi) {

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
	for (vector<PseudotreeNode*>::const_reverse_iterator it =
			ptnode->getChildren().rbegin(); it != ptnode->getChildren().rend();
			++it) {
		int vChild = (*it)->getVar();
		AobbSearchNodeOR* c = new AobbSearchNodeOR(n, vChild, depth + 1);
		// Compute and set heuristic estimate, includes child labels
		heuristicOR(c);
		c->setDepth(n->getDepth() + 1);
		chi.push_back(c);

	} // for loop over new OR children

	if (chi.empty()) {
		n->setLeaf(); // terminal node
		n->setValue(ELEM_ONE);
		return true; // no children
	}

	// order subproblems in decreasing order of their heuristic - largest UB first
	// (use reverse iterator due to stack reversal)
	if (m_problem->isMap(n->getVar())) {
		std::sort(chi.begin(), chi.end(), AobbSearchNode::heurGreater);
	}

	n->addChildren(chi);

	return false; // default

} // AOBB::generateChildrenAND

// DONE: radu
bool BB::generateChildrenOR(AobbSearchNode* n, vector<AobbSearchNode*>& chi) {

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

	// retrieve precomputed labels and heuristic values
	double* heur = n->getHeurCache();
	for (val_t i = m_problem->getDomainSize(var) - 1; i >= 0; --i) {
		// early pruning if heuristic is zero (since it's an upper bound)
		if (m_problem->isMap(var) && heur[2 * i] == ELEM_ZERO) { // 2*i=heuristic, 2*i+1=label
			continue;
		}

		AobbSearchNodeAND* c = new AobbSearchNodeAND(n, i, heur[2 * i + 1]); // uses cached label
		// set cached heur. value (includes the weight)
		c->setHeur(heur[2 * i]);
		c->setDepth(n->getDepth() + 1);
		chi.push_back(c);

		if (m_problem->isMap(n->getVar())) {
			assert(!ISNAN(heur[2*i]));
		}
	}

	if (chi.empty()) { // deadend
		n->setLeaf();
		n->setValue(ELEM_ZERO);
		return true; // no children
	}

	// sort new nodes by decreasing heuristic value - largest UB first
	// (use reverse iterator due to stack reversal)
	if (m_problem->isMap(n->getVar())) {
		sort(chi.begin(), chi.end(), AobbSearchNode::heurGreater);
	}

	n->addChildren(chi);

	return false; // default

} // AOBB::generateChildrenOR

// DONE: radu
AobbSearchNode* BB::nextNode() {
	AobbSearchNode* n = NULL;
	if (!n && m_stack.size()) {
		n = m_stack.top();
		m_stack.pop();
	}
	return n;
}

// DONE: radu
bool BB::doExpand(AobbSearchNode* n) {
	assert(n);
	m_expand.clear();

	if (n->getType() == NODE_AND) {  // AND node

//		if (m_problem->isMap(n->getVar()))
//			std::cout << "Expand " << n->toString() << std::endl;

		// update the heuristic
		std::vector<int> assignment;
		assignment.resize(m_assignment.size(), NONE);
		n->getPathAssignment(assignment);
		int currVar = n->getVar();
		if (m_problem->isMap(currVar)) {
			if (m_options->incremental == false) {
				m_heuristic->update(assignment);
			} else {
				n->setExpanded(); // mark the AND node as expanded
				m_heuristic->updateIncr(currVar, m_nextMapVar[currVar], assignment);
			}
		}

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
int BB::solve() {

	// check if solved during initialization
	if (m_solved) {
		std::cout << "--------- Solved during initialization ---------" << std::endl;
		std::cout << "Problem name:\t\t" << m_problem->getName() << std::endl;
		std::cout << "Status:\t\t\t" << solver_status[0] << std::endl;
		std::cout << "OR nodes:\t\t" << 0 << std::endl;
		std::cout << "AND nodes:\t\t" << 0 << std::endl;
		std::cout << "Cache hits:\t\t" << 0 << std::endl;
		std::cout << "Time elapsed:\t\t" << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
		std::cout << "Preprocessing:\t\t" << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
		std::cout << "Solution:\t\t" << m_solutionCost << " (" << ELEM_ENCODE(m_solutionCost) << ")" << std::endl;
		std::cout << "------------------------------------------------" << std::endl;

		return SEARCH_SUCCESS;
	}

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
	m_assignment.resize(m_problem->getN(), NONE);
	m_propagator.reset(new BoundPropagator(m_problem.get(),
			m_space.get(), m_pseudotree.get(), m_options, true, m_heuristic.get()));

	try {

		m_lowerBound = 0;

		// init cache table
		if (!m_space->cache) {
			m_space->cache = new AobbCacheTable(m_problem->getN());
		}

		// init root of the search space
		AobbSearchNode* first = initSearchSpace();
		if (first) {
			m_stack.push(first);
		}

		AobbSearchNode* n = nextLeaf();
		while (n) { // throws timeout
			m_propagator->propagate(n, m_timer, true); // true = report solutions
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
	std::cout << std::endl << std::endl;
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

