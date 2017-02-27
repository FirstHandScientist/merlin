/*
 * any_hbfs.cpp
 *
 *  Created on: 23 Sep 2015
 *      Author: radu
 */

#include "any_wrbfs.h"
#include "timer.h"
#include "utils.h"

#include "mini_bucket_heur.h"
#include "mini_bucket_mm_heur.h"
#include "mini_bucket_jglp_heur.h"


const float epsilon = 1e-10;

// Radu: DONE
// initialize
void AnyWRBFS::init() {

	// init problem instance
	Search::init();

	Timer tm;
	tm.start();

	// init heuristic generator
	if (m_options->heuristic == HEUR_MBE) {
		std::cout << "Computing MBE heuristic (i=" << m_options->ibound << ") ..." << std::endl;
		m_heuristic.reset(new MiniBucketHeur(m_problem.get(), m_pseudotree.get(),
				m_options, m_options->ibound));
	} else 	if (m_options->heuristic == HEUR_WMB_MM) {
		std::cout << "Computing WMB-MM heuristic (i=" << m_options->ibound << ") ..." << std::endl;
		m_heuristic.reset(new MiniBucketMMHeur(m_problem.get(), m_pseudotree.get(),
				m_options, m_options->ibound));
	} else 	if (m_options->heuristic == HEUR_WMB_JG) {
		std::cout << "Computing WMB-JG heuristic (i=" << m_options->ibound << ") ..." << std::endl;
		m_heuristic.reset(new MiniBucketJglpHeur(m_problem.get(), m_pseudotree.get(),
				m_options, m_options->ibound));
	}

	size_t sz = m_heuristic->build(&m_assignment, true);
	tm.stop();
	m_tmHeuristic = tm.elapsed();

	std::cout << "\tMini bucket finished in " << m_tmHeuristic
			<< " seconds" << std::endl;
	std::cout << "\tUsing " << (sz / (1024*1024.0)) * sizeof(double)
			<< " MBytes of RAM" << std::endl;

	// Now, we need to rearrange the pseudo tree such that the MAP variables
	// form a chain at the top of the pseudo tree (and they're searched in the
	// order suggested by the bucket tree heuristic ...
	std::list<int> order;
	m_pseudotree->dfs(order, true);
	m_pseudotree->forceMapChain(order, m_options->cbound);
	m_lastMapVar = order.back();
	assert( m_problem->isMap(m_lastMapVar) );

	// heuristic might have changed problem functions, pseudotree needs remapping
	m_pseudotree->resetFunctionInfo(m_problem->getFunctions());

	// check if heuristic is accurate
	if (m_heuristic->isAccurate()) {
		std::cout << "Heuristic is exact!" << std::endl;
		m_solutionCost = (m_heuristic->getGlobalUB());
	    m_solved = true;
	}

	// Record time after initialization
	cout << "Last MAP variable: " << m_lastMapVar << std::endl;
	cout << "Initialization complete: " << (m_tmLoad + m_tmHeuristic)
			<< " seconds" << std::endl;

	// output the pseudo tree
	m_pseudotree->outputToDot("pseudotree.dot");

	// check for upper bound computation only
	if (m_options->upperBoundOnly) {
		std::cout << "Done." << std::endl;
		exit(EXIT_SUCCESS);
	}
}

// Radu: DONE
AnyWRBFS::AnyWRBFS(ProgramOptions* opt) : Search(opt) {

	// Preallocate space for expansion vector.
	m_upperBound = ELEM_INF;
	m_lowerBound = ELEM_NAN;
	m_numExpanded = 0;
	m_numBacktracks = 0;
	m_prevTimeCheckpoint = 0;
	m_prevNodeCheckpoint = 0;
	m_overestimation = opt->overestimation;
	m_weight = opt->weight;
	m_numSumEvals = 0;
	m_lastMapVar = UNKNOWN;
}

// Radu: DONE
int AnyWRBFS::solve() {

	// check if solved during initialization
	if (m_solved) {
		std::cout << "--------- Solved during initialization ---------" << std::endl;
		std::cout << "Problem name:        " << m_problem->getName() << std::endl;
		std::cout << "Status:              " << solver_status[0] << std::endl;
		std::cout << "OR nodes:            " << 0 << std::endl;
		std::cout << "AND nodes:           " << 0 << std::endl;
		std::cout << "OR nodes (MAP):      " << 0 << std::endl;
		std::cout << "AND nodes (MAP):     " << 0 << std::endl;
		std::cout << "SUM evaluations:     " << 0 << std::endl;
		std::cout << "Time elapsed:        " << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
		std::cout << "Preprocessing:       " << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
		std::cout << "Solution:            " << m_solutionCost << " (" << ELEM_ENCODE(m_solutionCost) << ")" << std::endl;
		std::cout << "-------------------------------" << std::endl;

		return SEARCH_SUCCESS;
	}

	// init search space
	Timer tm;
	tm.start();
//	m_space.reset(new RbfaooSearchSpace(m_pseudotree.get(), m_options));
//	assert( m_options->cacheSize > 0 ); // in kilo bytes
//	size_t cache_kilobytes = (m_options->cacheSize);
//	if (!m_space->dfpncache) {
//		m_space->dfpncache = new RbfaooCacheTable;
//		m_space->dfpncache->init(cache_kilobytes);
//	}
	tm.stop();
	std::cout << "Cache allocation complete: " << tm.elapsed() << " seconds" << std::endl;

	// prologue
	std::cout << "--- Starting search ---" << std::endl;
	int res = SEARCH_SUCCESS;
	char* EmergencyMemory = new char[16384]; // in order to handle bad alloc

	m_tmSolve = 0;
	m_timer.reset(m_options->timeLimit);
	m_timer.start();
	m_solutionCost = ELEM_NAN;
	m_assignment.resize(m_problem->getN(), UNKNOWN);
	m_lowerBound = std::fabs( -ELEM_ENCODE(m_heuristic->getGlobalUB()) );

	try {

		// init the AO search space
		ao_space.reset(new AobbSearchSpace());
		ao_propagator.reset(new BoundPropagator(m_problem.get(), ao_space.get(), m_pseudotree.get(), m_options));
		ao_propagator->setSatSolver(NULL);
		ao_space->cache = new AobbCacheTable(m_problem->getN());


		// get the root of the pseudo tree
		int var = m_pseudotree->getRoot()->getVar();
		assert( m_problem->getDomainSize(var) == 1 ); // dummy variable
		int val = 0;

		// create the root node
		SearchNode root(var, val);
		double g = std::fabs( -ELEM_ENCODE(m_problem->getGlobalConstant()) );
		double h = std::fabs( -ELEM_ENCODE(m_heuristic->getGlobalUB()) );
		root.cost() = g;
		root.goal() = false;
		root.g() = g;
		root.h() = h;
		root.f() = g + h;
		root.F() = root.f();

		// do the RBFS search
		rbfs(root, root.F(), ELEM_INF, 0);
		m_solutionCost = m_upperBound;

	} catch (std::bad_alloc& ba) {
		delete[] EmergencyMemory;
		EmergencyMemory = NULL;
		std::cerr << "out of memory";
		res = SEARCH_FAILURE;
	} catch (...) {
		std::cout << "unexpected error: ";
		res = SEARCH_FAILURE;
	}

	// stop timer
	m_timer.stop();
	m_tmSolve = m_timer.elapsed();

	// output solution (if found)
	std::cout << std::endl << std::endl;
	std::cout << "--- Search done ---" << std::endl;
	std::cout << "Problem name:        " << m_problem->getName() << std::endl;
	std::cout << "Status:              " << solver_status[res] << std::endl;
//	std::cout << "OR nodes:            " << m_num_expanded_or << std::endl;
//	std::cout << "AND nodes:           " << m_num_expanded_and << std::endl;
	std::cout << "Nodes (MAP):         " << m_numExpanded << std::endl;
	std::cout << "SUM evaluations:     " << m_numSumEvals << std::endl;
	std::cout << "Time elapsed:        " << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
	std::cout << "Preprocessing:       " << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
//	std::cout << "Solution:            " << m_solutionCost << " (" << ELEM_ENCODE(m_solutionCost) << ")" << std::endl;
	std::cout << "Solution:            " << 1/ELEM_DECODE(m_solutionCost) << " (" << -m_solutionCost  << ")" << std::endl;
	std::cout << "-------------------------------" << std::endl;

	// clean up
	if (EmergencyMemory) delete[] EmergencyMemory;
	return res;
}

// Evaluate a goal node (ie, full MAP assignment)
double AnyWRBFS::eval(SearchNode& n) {
	assert(n.goal() && n.var() == m_lastMapVar);

	std::vector<int> backup = m_assignment;
	backup[n.var()] = n.val();
	double cost = ELEM_ZERO;
	const std::vector<PseudotreeNode*>& children = m_pseudotree->getNode(n.var())->getChildren();
	for (size_t i = 0; i < children.size(); ++i) {
		int ch = children[i]->getVar();
		double s = _ao(ch, backup);
		cost += std::fabs( -ELEM_ENCODE(s) );
	}
	m_assignment = backup;
	m_numSumEvals++;

	return cost;
}

double AnyWRBFS::rbfs(SearchNode& n, double f, double threshold, int depth) {

	double elapsed = m_timer.elapsed();
	if (elapsed - m_prevTimeCheckpoint >= SEARCH_CHECKPOINT_TIME) {

		double timestamp = elapsed;
		m_prevTimeCheckpoint = elapsed;
		std::cout << "[*"
			<< std::setw(8) << timestamp << "] w "
			<< std::setw(8) << m_weight << " "
			<< std::setw(12) << m_numExpanded << " "
			<< std::setw(12) << m_numBacktracks << " "
			<< std::setw(12) << m_numSumEvals << " "
			<< std::setw(12) << m_upperBound << " (" << 1/ELEM_DECODE(m_upperBound) << ") "
			<< std::setw(12) << m_lowerBound << " (" << 1/ELEM_DECODE(m_lowerBound) << ")"
			<< std::endl;
	}

	// check for goal
	if (n.goal()) { // complete MAP assignment
		return n.f();
	}

	double w = m_weight;

	// expand the node
	std::vector<SearchNode> succ = expand(n);
	assert(succ.empty() == false);

	for (size_t i = 0; i < succ.size(); ++i) {

		// check for solution
		if (succ[i].goal() && succ[i].f() < m_upperBound) {
			m_upperBound = succ[i].f();
			std::cout << "["
				<< std::setw(9) << m_timer.elapsed() << "] w "
				<< std::setw(8) << m_weight << " "
				<< std::setw(12) << m_numExpanded << " "
				<< std::setw(12) << m_numBacktracks << " "
				<< std::setw(12) << m_numSumEvals << " "
				<< std::setw(12) << m_upperBound << " (" << 1/ELEM_DECODE(m_upperBound) << ") "
				<< std::setw(12) << m_lowerBound << " (" << 1/ELEM_DECODE(m_lowerBound) << ")"
				<< std::endl;
		}

		if (succ[i].f() >= m_upperBound) {
			succ[i].F() = ELEM_INF; // prune
		} else if (n.f() < n.F()) {
			succ[i].F() = std::max(n.F(), succ[i].f());
		} else {
			succ[i].F() = succ[i].f();
		}

		//std::cout << " " << succ[i] << std::endl;
	}

	while (true) {
		int best = getBestFValueIndex(succ); // get the best F-value
		double bestCost = succ[best].g() + w*(succ[best].F() - succ[best].g()); // and the corresp. inflated value
		if (succ[best].F() >= m_upperBound || bestCost > threshold) {
			++m_numBacktracks;
			return succ[best].F();
		}

		int alt = getNextBestFValueIndex(succ, best); // get the second best F-value
		double altCost = succ[alt].g() + w*(succ[alt].F() - succ[alt].g()); // and the corresp. inflated value
		double new_threshold = std::min(threshold, altCost + m_overestimation); // new threshold
		//std::cout << "RBFS call on " << succ[best] << " with threshold " << new_threshold << std::endl;
		succ[best].F() = rbfs(succ[best], succ[best].F(), new_threshold, depth+1); // RBFS call
		//std::cout << "Backed up " << succ[best] << std::endl;
	}
}

size_t AnyWRBFS::getBestFValueIndex(std::vector<SearchNode>& nodes) {
	size_t lidx = 0;
	double lowest = ELEM_INF;
	for (size_t i = 0; i < nodes.size(); ++i) {
		if (nodes[i].F() < lowest) {
			lowest = nodes[i].F();
			lidx = i;
		}
	}

	return lidx;
}

size_t AnyWRBFS::getNextBestFValueIndex(std::vector<SearchNode>& nodes, size_t bestIndex) {
	size_t lidx = bestIndex;
	double lowest = ELEM_INF;
	for (size_t i = 0; i < nodes.size(); ++i) {
		if (i != bestIndex && nodes[i].F() < lowest) {
			lowest = nodes[i].F();
			lidx = i;
		}
	}

	return lidx;
}

// Expand a search node by generating its successors
std::vector<AnyWRBFS::SearchNode> AnyWRBFS::expand(SearchNode& n) {
	std::vector<SearchNode> succ;

	assert(n.goal() == false);
	int var = n.var();
	int val = n.val();
	m_assignment[var] = val;

	//std::cout << "Expanding " << n << std::endl;

	PseudotreeNode* ptn = m_pseudotree->getNode(var);
	const std::vector<PseudotreeNode*>& children = ptn->getChildren();
	assert(children.size() == 1); // OR chain over the MAP variables
	int child = children.at(0)->getVar();
	double* dv = heuristic(child);
	assert(dv != NULL);
	for (int k = 0; k < m_problem->getDomainSize(child); ++k) {
		SearchNode ch(child, k);
		ch.cost() = dv[2*k+1];
		ch.g() = n.g() + ch.cost();
		ch.h() = dv[2*k];
		ch.f() = ch.g() + ch.h();
		ch.goal() = (child == m_lastMapVar);
		if (ch.goal()) {
			double sval = eval(ch);
			ch.f() = ch.g() + sval;
			m_assignment[ch.var()] = UNKNOWN;
		}

		succ.push_back(ch);
	}

	m_numExpanded++;
	delete[] dv;
	return succ;
}

double* AnyWRBFS::heuristic(int var) {
	assert(m_problem->isMap(var));
	int oldValue = m_assignment[var];
	double* dv = new double[m_problem->getDomainSize(var) * 2];
	const list<Function*>& funs = m_pseudotree->getFunctions(var);
	for (val_t i = 0; i < m_problem->getDomainSize(var); ++i) {
		m_assignment[var] = i;

		// compute heuristic value
		dv[2 * i] = std::fabs( -ELEM_ENCODE(m_heuristic->getHeur(var, m_assignment)) );

		// precompute label value
		double d = ELEM_ONE;
		for (std::list<Function*>::const_iterator it = funs.begin();
				it != funs.end(); ++it) {
			d *= (*it)->getValue(m_assignment);
		}

		// store label and heuristic into cache table
		dv[2 * i + 1] = std::fabs( -ELEM_ENCODE(d) ); // label
		//dv[2 * i] += std::fabs( -ELEM_ENCODE(d) ); // heuristic
	}

	m_assignment[var] = oldValue;

	return dv;
}

////////////////////////////////////////////////////////////////////////////////
// solve the conditioned summation problem (SUM) via AOC search
double AnyWRBFS::_ao(int var, const std::vector<int>& assignment) {

	// a temporary buffer
	double result = ELEM_NAN;

	// reinit the search space
	ao_space->root = NULL;
	m_assignment = assignment;

	// init root of the search space
	_initAoSearch(var);

	AobbSearchNode* n = _nextLeaf();
	while (n) {
		ao_propagator->propagate(n, m_timer, m_currentDomains, false); // true = report solutions
		n = _nextLeaf();
	}

	// Conditioned SUM problem solved
	result = ao_space->root->getValue();
	delete ao_space->root;
	ao_space->root = NULL;

	return result;
}

// DONE: radu
void AnyWRBFS::_initAoSearch(int var) {

	// create root OR node
	AobbSearchNode* node = new AobbSearchNodeOR(NULL, var, 0);
	ao_space->root = node;
	_heuristicOR(node);

	assert(ao_stack.empty());
	ao_stack.push(node);
}

// DONE: radu
void AnyWRBFS::_heuristicOR(AobbSearchNode* n) {

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
bool AnyWRBFS::_doProcess(AobbSearchNode* node) {
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
void AnyWRBFS::_addCacheContext(AobbSearchNode* node, const set<int>& ctxt) const {

	context_t sig;
	sig.reserve(ctxt.size());
	for (set<int>::const_iterator itC = ctxt.begin(); itC != ctxt.end(); ++itC) {
		sig.push_back(m_assignment[*itC]);
	}

	node->setCacheContext(sig);
}

// DONE: radu
bool AnyWRBFS::_doCaching(AobbSearchNode* node) {

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
AobbSearchNode* AnyWRBFS::_nextLeaf() {

	AobbSearchNode* node = _nextNode();
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
		if (_doExpand(node)) { // node expansion
			return node;
		}
		node = _nextNode();
	}
	return NULL;
}

// DONE: radu
bool AnyWRBFS::_generateChildrenAND(AobbSearchNode* n, vector<AobbSearchNode*>& chi) {

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
bool AnyWRBFS::_generateChildrenOR(AobbSearchNode* n, vector<AobbSearchNode*>& chi) {

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
AobbSearchNode* AnyWRBFS::_nextNode() {
	AobbSearchNode* n = NULL;
	if (!n && ao_stack.size()) {
		n = ao_stack.top();
		ao_stack.pop();
	}
	return n;
}

// DONE: radu
bool AnyWRBFS::_doExpand(AobbSearchNode* n) {
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


