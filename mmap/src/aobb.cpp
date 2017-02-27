/*
 * aobb.cpp
 *
 *  Created on: Jun 7, 2013
 *      Author: radu
 */


#include "aobb.h"
#include "mini_bucket_heur.h"
#include "mini_bucket_mm_heur.h"
#include "mini_bucket_jglp_heur.h"
#include "mini_bucket_wmbjg_heur.h"

// initialize
void AOBB::init() {

	// init problem instance
	Search::init();

	Timer tm;
	tm.start();

	// Encode determinism
	encodeDeterminism();

#ifdef UAI_COMPETITION
	// determine the optimal i-bound given 4GB of RAM
	size_t memlimit = m_options->memoryLimit*1024;
	MiniBucketHeur *temp = new MiniBucketHeur(m_problem.get(),
			m_pseudotree.get(), m_options, m_options->ibound);
	temp->limitSize(memlimit, &m_assignment);
	delete temp;

	// slightly smaller ibound for WMB-JG
	if (m_options->heuristic == HEUR_WMB_JG || m_options->heuristic == HEUR_WMB_JG2) {
		int ib = m_options->ibound;
		int delta = (int)floor((ib*0.3));
		ib -= delta;
		if (ib < 2) ib = 2; // default
		m_options->ibound = ib;
	}

	// get the static order of the MAP variables and re-arrange the pseudo tree
	std::vector<int> elim = m_elimOrder;
	std::list<int> mapSearchOrder;
	mapSearchOrder.push_back(m_pseudotree->getRoot()->getVar()); // the dummy var
	for (std::vector<int>::reverse_iterator ri = elim.rbegin(); ri != elim.rend(); ++ri) {
		if (m_problem->isMap(*ri)) mapSearchOrder.push_back(*ri);
	}
	std::cout << "MAP search order: ";
	std::copy(mapSearchOrder.begin(), mapSearchOrder.end(), std::ostream_iterator<int>(std::cout, " "));
	std::cout << std::endl;
//	m_pseudotree->forceMapChain(mapSearchOrder, m_options->cbound);
//	m_pseudotree->resetFunctionInfo(m_problem->getFunctions());
//	std::cout << "Forcing chain pseudo tree for MAP variables done." << std::endl;

#endif

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
	} else 	if (m_options->heuristic == HEUR_WMB_JG2) {
		std::cout << "Computing WMB-JG heuristic (i=" << m_options->ibound << ") ..." << std::endl;
		m_heuristic.reset(new MiniBucketJGHeur(m_problem.get(), m_pseudotree.get(),
				m_options, m_options->ibound));
		m_options->upperBoundOnly = true;
	}

	size_t sz = m_heuristic->build(&m_assignment, true);
	tm.stop();
	m_tmHeuristic = tm.elapsed();

	std::cout << "\tMini bucket finished in " << m_tmHeuristic
			<< " seconds" << std::endl;
	std::cout << "\tUsing " << (sz / (1024*1024.0)) * sizeof(double)
			<< " MBytes of RAM" << std::endl;

	// heuristic might have changed problem functions, pseudotree needs remapping
	m_pseudotree->resetFunctionInfo(m_problem->getFunctions());

	// check if heuristic is accurate
	if (m_heuristic->isAccurate()) {
		std::cout << "Heuristic is exact!" << std::endl;
		m_solutionCost = (m_heuristic->getGlobalUB());
	    m_solved = true;
	}

	// Record time after initialization
	cout << "Initialization complete: " << (m_tmLoad + m_tmHeuristic)
			<< " seconds" << std::endl;

	// check for upper bound computation only
	if (m_options->upperBoundOnly) {
		std::cout << "Done." << std::endl;
		exit(EXIT_SUCCESS);
	}
}

// DONE: radu
AobbSearchNode* AOBB::initSearchSpace() {
	assert(m_space->root == NULL);

	// Add initial set of dummy nodes.

	// create root OR node (dummy variable)
	PseudotreeNode* ptroot = m_pseudotree->getRoot();
	AobbSearchNode* node = new AobbSearchNodeOR(NULL, ptroot->getVar(), -1);
	m_space->root = node;
//	double h = heuristicOR(node);
	double h = m_heuristic->getGlobalUB();
	node->setHeur(h);

	// create dummy AND node (domain size 1) with global constant as label
	AobbSearchNode* next = new AobbSearchNodeAND(node, 0, m_problem->getGlobalConstant());
	m_space->root->setChild(next);
//	next->setHeur(node->getHeurCache()[0]);
	next->setHeur( h/next->getLabel() );

	return next;
}

// DONE: radu
double AOBB::heuristicOR(AobbSearchNode* n) {

	// Safety checks
	assert(n->getType() == NODE_OR);

	// check if MAP variable
	bool mapVar = m_problem->isMap(n->getVar());

	int var = n->getVar();
	double* dv = new double[m_problem->getDomainSize(var) * 2];
	const list<Function*>& funs = m_pseudotree->getFunctions(var);
	double h = (mapVar) ? -INFINITY : ELEM_ZERO; // the new OR nodes h value
	for (val_t i = 0; i < m_problem->getDomainSize(var); ++i) {
		m_assignment[var] = i;

		// compute heuristic value
		dv[2 * i] = m_heuristic->getHeur(var, m_assignment);

		// precompute label value
		double d = ELEM_ONE;
		for (std::list<Function*>::const_iterator it = funs.begin();
				it != funs.end(); ++it) {
			d *= (*it)->getValue(m_assignment);
		}

		// store label and heuristic into cache table
		dv[2 * i + 1] = d; // label
		dv[2 * i] *= d; // heuristic (includes label)

		if (mapVar) {
			if (dv[2 * i] > h) {
				h = dv[2 * i]; // keep max. for OR node heuristic (MAP var)
			}
		} else {
			h += dv[2 * i]; // keep sum. for OR node heuristic (SUM var)
		}
	}

	n->setHeur(h);
	n->setHeurCache(dv);

	return h;

} // AOBB::heuristicOR


// DONE: radu
bool AOBB::doProcess(AobbSearchNode* node) {
	assert(node);
	if (node->getType() == NODE_AND) {
		int var = node->getVar();
		int val = node->getVal();
		m_assignment[var] = val; // record assignment

/*
		// check if we have a full MAP assignment (new)
		bool full = true;
		for (int i = 0; i < m_problem->getN(); ++i) {
			if (m_problem->isMap(i) && m_assignment[i] == UNKNOWN) {
				full = false;
				break;
			}
		}

		if (full) {

			bool equal = true;
			if (m_prevAssignment.empty()) {
				equal = false;
			} else {
				for (int i = 0; i < m_problem->getN(); ++i) {
					if (m_problem->isMap(i)) {
						if (m_assignment[i] != m_prevAssignment[i]) {
							equal = false;
							break;
						}
					}
				}
			}

			if (equal == false) {
				m_prevAssignment = m_assignment;
				std::vector<val_t> sol = m_assignment;

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
					std::ofstream out(output.c_str(), std::ios::app);
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
		}
*/
	} else { // NODE_OR
		// do nothing
	}

	return false; // default
}

// DONE: radu
void AOBB::addCacheContext(AobbSearchNode* node, const set<int>& ctxt) const {

	context_t sig;
	sig.reserve(ctxt.size());
	for (set<int>::const_iterator itC = ctxt.begin(); itC != ctxt.end(); ++itC) {
		sig.push_back(m_assignment[*itC]);
	}

	node->setCacheContext(sig);
}

// DONE: radu
bool AOBB::doCaching(AobbSearchNode* node) {

	assert(node);
	int var = node->getVar();
	PseudotreeNode* ptnode = m_pseudotree->getNode(var);

	if (node->getType() == NODE_AND) { // AND node -> reset associated adaptive cache tables

		// if adaptive caching then purge corresponding caches
	    const list<int>& resetList = ptnode->getCacheReset();
	    for (list<int>::const_iterator it = resetList.begin();
	    		it != resetList.end(); ++it) {
	    	m_space->cache->reset(*it);
	    }

	} else { // OR node, try actual caching

		if (!ptnode->getParent())
			return false;

		if (ptnode->getFullContext().size() <=
			ptnode->getParent()->getFullContext().size()) {

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
bool AOBB::doPruning(AobbSearchNode* node) {

	assert(node);

	if (canBePruned(node)) {
		node->setLeaf();
		node->setPruned();
		if (node->getType() == NODE_OR) {
			if (ISNAN(node->getValue())) // value could be set by LDS
				node->setValue(ELEM_ZERO);
		} else if (node->getType() == NODE_AND) {
			node->setValue(ELEM_ZERO); // dead end
		}

		return true;
	}

	return false; // default false
} // AOBB::doPruning

// DONE: radu
AobbSearchNode* AOBB::nextLeaf() {

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
					<< std::setw(12) << m_lowerBound << " (" << -(ELEM_ENCODE(m_lowerBound)) << ")"
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

void AOBB::getPathAssignment(AobbSearchNode* node, std::vector<int>& vars) {

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

bool AOBB::canBePruned(AobbSearchNode* n) {

	// propagate current assignment (AND nodes only)
	if (n->getType() == NODE_AND) {
		if (m_options->propagation == CP_FULL) {

			// propagate the assignment
			std::vector<int> vars;
			getPathAssignment(n, vars);

			bool conflict = propagate(vars);
			if (conflict) {
				++m_deadendsCP;
				return true;
			}
		} else if (m_options->propagation == CP_UNIT) {

			// Propagate the assignment
			m_zchaff.dlevel()++;
			if (!lookahead(n->getVar(), n->getVal(), n->changes(),
					m_pseudotree->getNode(n->getVar())->getSubprobVars())) {
				// found inconsistency
				++m_deadendsCP;
				return true;
			}
		}
	}

	// do not prune SUM nodes
	if (m_problem->isMap(n->getVar()) == false) {
		return false;
	}

	// heuristic is an upper bound, hence can use to prune if value=0
	if (n->getHeur() == ELEM_ZERO) {
		++m_deadendsUB;
		return true;
	}

	AobbSearchNode* curAND;
	AobbSearchNode* curOR;
	double curPSTVal;

	if (n->getType() == NODE_AND) {
		curAND = n;
		curOR = n->getParent();
		curPSTVal = curAND->getHeur(); // includes label
	} else { // NODE_OR
		curAND = NULL;
		curOR = n;
		curPSTVal = curOR->getHeur(); // n->getHeur()
	}

	list<AobbSearchNode*> notOptOR; // marks nodes for tagging as possibly not optimal

	// up to root node, if we have to
	while (curOR->getParent()) {

		//if ( fpGt(curPSTVal, curOR->getValue()) ) {
		if ( curPSTVal <= curOR->getValue() ) {
			for (list<AobbSearchNode*>::iterator it=notOptOR.begin(); it!=notOptOR.end(); ++it)
				if (m_problem->isMap((*it)->getVar()))
					(*it)->setNotOpt(); // mark as possibly not optimal
			++m_deadendsUB;
			return true;// pruning is possible!
		}

		notOptOR.push_back(curOR);

		// climb up, update values
		curAND = curOR->getParent();

		// collect AND node label
		curPSTVal *= curAND->getLabel();
		// incorporate already solved sibling OR nodes
		curPSTVal *= curAND->getSubSolved();
		// incorporate new not-yet-solved sibling OR nodes through their heuristic
		NodeP* children = curAND->getChildren();
		for (size_t i = 0; i < curAND->getChildCountFull(); ++i) {
			if (!children[i] || children[i] == curOR) continue;
			else curPSTVal *= children[i]->getHeur();
		}
		curOR = curAND->getParent();
	}

	// default, no pruning possible
	return false;

} // AOBB::canBePruned

// DONE: radu
bool AOBB::generateChildrenAND(AobbSearchNode* n, vector<AobbSearchNode*>& chi) {

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
bool AOBB::generateChildrenOR(AobbSearchNode* n, vector<AobbSearchNode*>& chi) {

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

	// check if subproblem is solved exactly (by MBE)
	if ( m_options->exactSubprob ) { // check if switch is on
		PseudotreeNode* ptnode = m_pseudotree->getNode(var);
		int tw = ptnode->getSubWidth();
		if (tw <= m_options->ibound) {
			n->setValue(n->getHeur());
			n->setCachable();
			n->setLeaf();
			return true; // no children generated below it
		}
	}

	// increase OR node expansions
	m_space->nodesOR += 1;
	if (m_problem->isMap(var))
		m_space->nodesORmap += 1;

	// retrieve precomputed labels and heuristic values
	double* heur = n->getHeurCache();
	for (val_t i = m_problem->getDomainSize(var) - 1; i >= 0; --i) {

		// consistent values only
		if (m_currentDomains.empty() == false && m_currentDomains[var][i] == false) {
			continue;
		}

		// early pruning if heuristic is zero (since it's an upper bound)
		if (m_problem->isMap(var) && heur[2 * i] == ELEM_ZERO) { // 2*i=heuristic, 2*i+1=label
			continue;
		}

		AobbSearchNodeAND* c = new AobbSearchNodeAND(n, i, heur[2 * i + 1]); // uses cached label
		// set cached heur. value (includes the weight)
		c->setHeur(heur[2 * i]);
		chi.push_back(c);
	}

	if (chi.empty()) { // deadend
		n->setLeaf();
		n->setValue(ELEM_ZERO);
		return true; // no children
	}

	// sort new nodes by decreasing order of heuristic value - largest UB first
	// (use reverse iterator due to stack reversal)
	if (m_problem->isMap(n->getVar())) {
		sort(chi.begin(), chi.end(), AobbSearchNode::heurGreater);
	}

	n->addChildren(chi);

	return false; // default

} // AOBB::generateChildrenOR

// DONE: radu
AobbSearchNode* AOBB::nextNode() {
	AobbSearchNode* n = NULL;
	if (!n && m_stack.size()) {
		n = m_stack.top();
		m_stack.pop();
	}
	return n;
}

// DONE: radu
bool AOBB::doExpand(AobbSearchNode* n) {
	assert(n);
	m_expand.clear();

	if (n->getType() == NODE_AND) {  // AND node

		if (generateChildrenAND(n, m_expand))
			return true; // no children
		for (vector<AobbSearchNode*>::reverse_iterator it = m_expand.rbegin();
				it != m_expand.rend(); ++it) {
			m_stack.push(*it);
		}
		//std::cout << "[DEBUG] expanded " << n->toString() << endl;
	} else {  // OR node

		if (generateChildrenOR(n, m_expand))
			return true; // no children
		for (vector<AobbSearchNode*>::reverse_iterator it = m_expand.rbegin();
				it != m_expand.rend(); ++it) {
			m_stack.push(*it);
		} // for loop
		//std::cout << "[DEBUG] expanded " << n->toString() << endl;
	} // if over node type

	return false; // default false
}

// solve the problem
int AOBB::solve() {

#ifndef UAI_COMPETITION
	// check if solved during initialization
	if (m_solved) {
		std::cout << "--------- Solved during initialization ---------" << std::endl;
		std::cout << "Problem name:        " << m_problem->getName() << std::endl;
		std::cout << "Status:              " << solver_status[0] << std::endl;
		std::cout << "OR nodes:            " << 0 << std::endl;
		std::cout << "AND nodes:           " << 0 << std::endl;
		std::cout << "Cache hits:          " << 0 << std::endl;
		std::cout << "SUM evaluations:     " << 0 << std::endl;
		std::cout << "Time elapsed:        " << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
		std::cout << "Preprocessing:       " << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
		std::cout << "Solution:            " << m_solutionCost << " (" << ELEM_ENCODE(m_solutionCost) << ")" << std::endl;
		std::cout << "------------------------------------------------" << std::endl;

		return SEARCH_SUCCESS;
	}
#endif

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

	// search
	try {

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
			m_propagator->propagate(n, m_timer, m_currentDomains, true); // true = report solutions
			m_lowerBound = m_propagator->getLowerBound();
			n = nextLeaf();
		}

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
	m_solutionCost = m_propagator->getLowerBound();

	// output solution (if found)
	std::cout << std::endl << std::endl;
	std::cout << "----------------- Search done ------------------" << std::endl;
	std::cout << "Problem name:        " << m_problem->getName() << std::endl;
	std::cout << "Status:              " << solver_status[res] << std::endl;
	std::cout << "OR nodes:            " << m_space->nodesOR << std::endl;
	std::cout << "AND nodes:           " << m_space->nodesAND << std::endl;
	std::cout << "OR nodes (MAP):      " << m_space->nodesORmap << std::endl;
	std::cout << "AND nodes (MAP):     " << m_space->nodesANDmap << std::endl;
	std::cout << "Cache hits:          " << m_cacheHits << std::endl;
	std::cout << "Deadends (CP):       " << m_deadendsCP << std::endl;
	std::cout << "Deadends (UB):       " << m_deadendsUB << std::endl;
	std::cout << "SUM evaluations:     " << m_space->numSumEvals << std::endl;
	std::cout << "Time elapsed:        " << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
	std::cout << "Preprocessing:       " << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
	std::cout << "Solution:            " << m_solutionCost << " (" << ELEM_ENCODE(m_solutionCost) << ")" << std::endl;
	std::cout << "------------------------------------------------" << std::endl;

	// clean up
	if (EmergencyMemory) delete[] EmergencyMemory;

	return res;
}

