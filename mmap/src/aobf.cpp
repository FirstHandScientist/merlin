/*
 * aobf.cpp
 *
 *  Created on: Apr 17, 2013
 *      Author: radu
 */



#include "aobf.h"
#include "timer.h"
#include "utils.h"

#include "mini_bucket_heur.h"
#include "mini_bucket_mm_heur.h"
#include "mini_bucket_jglp_heur.h"

//#define DEBUG

// initialize
void AOBF::init() {

	// init problem instance
	Search::init();

	Timer tm;
	tm.start();

	// Encode determinism
	encodeDeterminism();

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

	// heuristic might have changed problem functions, pseudotree needs remapping
	m_pseudotree->resetFunctionInfo(m_problem->getFunctions());

	m_pseudotree->outputToDot("pseudotree.dot");

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

// solve the problem
int AOBF::solve() {

	// check if solved during initialization
	if (m_solved) {
		std::cout << "--------- Solved during initialization ---------" << std::endl;
		std::cout << "Problem name:        " << m_problem->getName() << std::endl;
		std::cout << "Status:              " << solver_status[0] << std::endl;
		std::cout << "OR nodes:            " << 0 << std::endl;
		std::cout << "AND nodes:           " << 0 << std::endl;
		std::cout << "OR nodes (MAP):      " << 0 << std::endl;
		std::cout << "AND nodes (MAP):     " << 0 << std::endl;
		std::cout << "Deadends (CP):       " << 0 << std::endl;
		std::cout << "SUM evaluations:     " << 0 << std::endl;
		std::cout << "Time elapsed:        " << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
		std::cout << "Preprocessing:       " << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
		std::cout << "Solution:            " << m_solutionCost << " (" << ELEM_ENCODE(m_solutionCost) << ")" << std::endl;
		std::cout << "-------------------------------" << std::endl;

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
	m_space.reset(new AobfSearchSpace());
	m_solutionCost = ELEM_NAN;
	m_assignment.resize(m_problem->getN(), UNKNOWN);

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

		// AO* search
		res = aostar(m_options->verbose);

		m_solved = true;
	} catch (std::bad_alloc& ba) {
		delete[] EmergencyMemory;
		EmergencyMemory = NULL;
		std::cerr << "OUT_OF_MEMORY";
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
	m_solutionCost += -ELEM_ENCODE(m_problem->getGlobalConstant());

	// output solution (if found)
	std::cout << std::endl << std::endl;
	std::cout << "--- Search done ---" << std::endl;
	std::cout << "Problem name:        " << m_problem->getName() << std::endl;
	std::cout << "Status:              " << solver_status[res] << std::endl;
	std::cout << "OR nodes:            " << ao_space->nodesOR + m_space->getOrNodes() << std::endl;
	std::cout << "AND nodes:           " << ao_space->nodesAND + m_space->getAndNodes() << std::endl;
	std::cout << "OR nodes (MAP):      " << m_space->getOrNodes() << std::endl;
	std::cout << "AND nodes (MAP):     " << m_space->getAndNodes() << std::endl;
	std::cout << "Cache hits (MAP):    " << m_numCacheHits << std::endl;
	std::cout << "Deadends (CP):       " << m_deadendsCP << std::endl;
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

void AOBF::getPathAssignment(AobfSearchNode* node, std::vector<int>& vars) {

	assert(node->getType() == NODE_OR);
	if (node->getCurrentParent() == NULL) return;
	AobfSearchNode *crt = node->getCurrentParent(), *par = crt->getCurrentParent();
	while (true) {
		int vvar = crt->getVar();
		int vval = crt->getVal();
		int vid = m_var2sat[vvar][vval];
		vars.push_back(vid);

		crt = par->getCurrentParent(); // AND
		if (crt == NULL) break;
		par = crt->getCurrentParent(); // OR
	}

}

bool AOBF::isDeadend(AobfSearchNode* n, int var, int val) {

	// propagate current assignment (AND nodes only)
	if (n->getType() == NODE_OR) {

		if (m_options->propagation == CP_FULL) {

			// propagate the assignment
			std::vector<int> vars;
			getPathAssignment(n, vars);
			vars.push_back(m_var2sat[var][val]);

			bool conflict = propagate(vars);
			if (conflict) {
				++m_deadendsCP;
				return true;
			}

		}
	}

	// default, no pruning possible
	return false;

} // BBBT::canBePruned

// Revise the value of node 'n' based on its children values; log arithmetic
bool AOBF::revise(AobfSearchNode *node) {

	// safety checks
	assert(node);
	assert(node->getType() == NODE_AND || node->getType() == NODE_OR);
	assert(m_problem->isMap(node->getVar()));

	bool change = true;
	if (node->getType() == NODE_AND) { // revise the cost of an AND node
		if (node->isTerminal()) {

			// terminal AND node
			node->setValue(ELEM_ZERO);
			node->setSolved(true);
			node->setFringe(false);

			change = true;
		} else {

			// non-terminal AND node
			double oldValue = node->getValue();
			bool solved = true;
			double qval = ELEM_ZERO;
			std::list<AobfSearchNode*>::const_iterator li =
					node->getChildren().begin();
			for (; li != node->getChildren().end(); ++li) {
				AobfSearchNode *ch = (*li);
				bool sol = ch->isSolved();
				solved &= sol;
				qval += ch->getValue();
			}

			// label the node 'solved' if all its OR children are 'solved'
			node->setValue(qval);
			if (solved) {
				node->setSolved(true);
				node->setFringe(false);
			}

			change = ( solved || (qval != oldValue) );
		}
	} else { // revise the cost of an OR node

		double oldValue = node->getValue();
		double qval = INFINITY;
		AobfSearchNode *best = NULL;

		std::list<AobfSearchNode*>::const_iterator li = node->getChildren().begin();
		for (; li != node->getChildren().end(); ++li) {
			AobfSearchNode *ch = (*li);
			int val = ch->getVal();
			double w = node->getWeight(val);
			double q = w + ch->getValue(); // could be infinity

			if (q + DOUBLE_PRECISION < qval) {
				qval = q;
				best = ch;
			} else if ( q == qval /*equals(q, qval)*/) {
				if (best == NULL) {
					best = ch;
				}
			}
		}

		// safety checks (OR nodes must always be marked)
		assert( best != NULL );

		// label the node 'solved' if the marked AND child is 'solved'
		node->setValue(qval);
		node->setBestChild(best);
		bool solved = best->isSolved();
		if (solved) {
			node->setSolved(true);
			node->setFringe(false);
		}

		change = ( solved || (qval != oldValue) );
	}

	return change;
}

// returns a string of the context instantiation
std::string AOBF::context(int var, int val, const std::set<int>& C) {

	std::ostringstream oss;
	if ( C.empty() || val == NONE) { // root and dummy node
		oss << "s" << m_globalSearchIndex++;
	} else {
		for (std::set<int>::const_iterator si = C.begin(); si != C.end(); ++si) {
			int cvar = (*si);
			if (cvar == var) {
				oss << "x" << var << "=" << val << ";";
			} else {
				assert(m_assignment[cvar] != UNKNOWN);
				oss << "x" << cvar << "=" << m_assignment[cvar] << ";";
			}
		}
	}

	return oss.str();
}

// expand a search node by generating its successors
bool AOBF::expand(AobfSearchNode *node) {

	// safety checks
	assert(node);
	assert(node->getType() == NODE_AND || node->getType() == NODE_OR);

	bool noChildren = true;
	if (node->getType() == NODE_AND) { // expand an AND node

		int var = node->getVar();
		int depth = node->getDepth();
		PseudotreeNode* ptNode = m_pseudotree->getNode(var);
		const std::vector<PseudotreeNode*>& children = ptNode->getChildren();
		std::vector<PseudotreeNode*>::const_iterator it = children.begin();
		for (; it != children.end(); ++it) {
			int vChild = (*it)->getVar();
			AobfSearchNodeOR* c = new AobfSearchNodeOR(vChild, depth + 1);
			c->addParent(node); // only 1 parent per OR node
			node->addChild(c);

			if (m_problem->isMap(vChild)) {
				double h = heuristicOR(c); // log space
				c->setHeur( h );
				c->setValue( m_epsilon*h );
			} else { // conditioned SUM problem rooted at an OR node
				assert(m_problem->isSum(vChild));
				std::vector<int> backup = m_assignment;
				double qval = _ao(vChild, m_assignment);
				c->setHeur(ELEM_NAN);
				c->setValue( -ELEM_ENCODE(qval) );
				c->setTerminal(true);
				c->setSolved(true);
				c->setFringe(false);
				c->setExpanded(true);
				m_assignment = backup;
				++m_numSumEvals;
			}

			// add the OR child to the search space
			std::string str = context(vChild, NONE, ptNode->getFullContext());
			AobfSearchState state(NODE_OR, str);
			m_space->add(state, c);
			noChildren = false;
		}

		// set the flags: expanded, fringe, terminal
		node->setExpanded(true);
		node->setFringe(false);
		node->setTerminal( children.empty() );
		m_space->incNodesExpanded(NODE_AND);
	} else { // expand an OR node
		int var = node->getVar();
		int depth = node->getDepth();
		PseudotreeNode* ptNode = m_pseudotree->getNode(var);
		assert( m_problem->isMap(var) == true );
		assert( node->getHeurCache() != NULL );
		for (int val = 0; val < m_problem->getDomainSize(var); ++val) {

			// check if deadend (due to constraint propagation)
			if (isDeadend(node, var, val)) {
				++m_deadendsCP;
				continue;
			}

			// create the state corresponding to the AND child
			std::string str = context(var, val, ptNode->getAndContext());
			AobfSearchState state(NODE_AND, str);
			if (m_space->find(state) == true) {
				AobfSearchNodeAND* c = (AobfSearchNodeAND*) m_space->get(state);
				node->addChild(c);
				c->addParent(node);
				++m_numCacheHits;
			} else {
				AobfSearchNodeAND* c = new AobfSearchNodeAND(var, val, depth + 1);
				node->addChild(c);
				c->addParent(node);
				double *heurCache = node->getHeurCache();
				double h = heurCache[2*val]; // includes the weight
				double w = heurCache[2*val+1];
				c->setHeur(h - w);
				c->setValue(m_epsilon*(h - w));
				m_space->add(state, c);
			}

			noChildren = false;
		}

		// set the flags: expanded, fringe, terminal
		node->setExpanded(true);
		node->setFringe(false);
		node->setTerminal(false);
		m_space->incNodesExpanded(NODE_OR);
	}

	return noChildren; // "true" if no children, "false" otherwise
}

// DONE: radu
double AOBF::heuristicOR(AobfSearchNode* n) {

	// Safety checks
	assert(n->getType() == NODE_OR);

	// check if MAP variable
	int var = n->getVar();
	bool mapVar = m_problem->isMap(var);
	int oldValue = m_assignment[var];
	double* dv = new double[m_problem->getDomainSize(var) * 2];
	const list<Function*>& funs = m_pseudotree->getFunctions(var);
	double h = (mapVar) ? INFINITY : ELEM_ZERO; // the new OR nodes h value
	for (val_t i = 0; i < m_problem->getDomainSize(var); ++i) {
		m_assignment[var] = i;

		// compute heuristic value
		dv[2 * i] = -ELEM_ENCODE(m_heuristic->getHeur(var, m_assignment));

		// precompute label value
		double d = ELEM_ONE;
		for (std::list<Function*>::const_iterator it = funs.begin();
				it != funs.end(); ++it) {
			d *= (*it)->getValue(m_assignment);
		}

		// store label and heuristic into cache table
		dv[2 * i + 1] = -ELEM_ENCODE(d); // label
		dv[2 * i] += -ELEM_ENCODE(d); // heuristic (includes label)

		if (mapVar) {
			if (dv[2 * i] < h) {
				h = dv[2 * i]; // keep min. for OR node heuristic (MAP var)
			}
		} else {
			h += dv[2 * i]; // keep sum. for OR node heuristic (SUM var)
		}
	}

	n->setHeur(h);
	n->setHeurCache(dv);
	m_assignment[var] = oldValue;

	return h;

} // AOBF::heuristicOR

// expand and update a search node
void AOBF::expandAndRevise(AobfSearchNode *node) {

	std::string sol;
	std::multiset<AobfSearchNode*, CompNodeIndexAsc> S;

	// expand node n: Steps 6 and 7 from Nilsson's
	assert( node->isFringe() && !node->isSolved() );
	expand(node);

#ifdef DEBUG
	std::cout << "Expanded "<< node->toString() << std::endl;
#endif

	S.insert(node);
	node->setVisited(true);

	// iterate over S: Step 10 from Nilsson's
	while ( !S.empty() ) {
		AobfSearchNode *e = *S.begin();
		S.erase(S.begin());
		e->setVisited(false);

		assert( e->getIndex() == 0 );
		bool change = revise(e);
#ifdef DEBUG
	std::cout << "Revised "<< e->toString() << std::endl;
#endif

		if (change) {

			// update parents: Step 12 from Nilsson's
			if (e->getType() == NODE_AND) {
				// AND nodes may have multiple parents in the context-minimal AND/OR graph.
				std::list<AobfSearchNode*>::iterator pi = e->getParents().begin();
				for (; pi != e->getParents().end(); ++pi) {
					AobfSearchNode *p = (*pi);
					assert(m_problem->isMap(p->getVar()));

					size_t index = 0;
					AobfSearchNode *best = p->getBestChild();
					bool found = (best == e);	// check if 'e' is the marked AND child of 'p'
					std::list<AobfSearchNode*>::iterator ci = p->getChildren().begin();
					for (; ci != p->getChildren().end(); ++ci) {
						AobfSearchNode *c = (*ci);
						if (c->isVisited()) {
							++index;
						}
					}

					if (p->isVisited()) {
						// if the parent already in S, decrease the index.
						std::multiset<AobfSearchNode*, CompNodeIndexAsc>::iterator si = S.begin();
						for (; si != S.end(); ++si) {
							if (p == *si) {
								/*si = */S.erase(si);
								break;
							}
						}
						p->decIndex();
						S.insert(p);
					} else if (found) {
						// new parent through marked connector
						p->setIndex(index);
						p->setVisited(true);
						S.insert(p);
					}
				}
			} else {
				// OR nodes have only one parent in the context-minimal AND/OR graph.
				std::list<AobfSearchNode*>::iterator pi = e->getParents().begin();
				for (; pi != e->getParents().end(); ++pi) {
					AobfSearchNode *p = (*pi);
					assert(m_problem->isMap(p->getVar()));

					size_t index = 0;
					std::list<AobfSearchNode*>::iterator ci = p->getChildren().begin();
					for (; ci != p->getChildren().end(); ++ci) {
						AobfSearchNode *c = (*ci);
						if (c->isVisited()) {
							++index;
						}
					}

					p->setIndex(index);
					p->setVisited(true);
					S.insert(p);
				}
			}
		} else {

			// update the index of potential visited parents
			std::list<AobfSearchNode*>::iterator pi = e->getParents().begin();
			for (; pi != e->getParents().end(); ++pi) {
				AobfSearchNode *p = (*pi);
				assert(m_problem->isMap(p->getVar()));

				if (p->isVisited()) {

					// if the parent already in S, decrease the index
					std::multiset<AobfSearchNode*, CompNodeIndexAsc>::iterator si = S.begin();
					for (; si != S.end(); ++si) {
						if (p == *si) {
							/*si = */S.erase(si);
							break;
						}
					}
					p->decIndex();
					S.insert(p);
				}
			}
		} // end if
	}

	assert( S.empty() );

	// compute best partial solution tree: Step 4 from Nilsson's
	m_tips.clear();
	findBestPartialTree();

	assert( m_space->getRoot()->isSolved() || (m_tips.size() > 0) );
}

// find the new best partial solution tree
bool AOBF::findBestPartialTree() {

	// get the root of the search space
	AobfSearchNode *root = m_space->getRoot();

	// clear the previous value assignment.
	std::fill(m_assignment.begin(), m_assignment.end(), UNKNOWN);

	// follow the markings
	std::stack<AobfSearchNode*> s;
	//root->setCurrentParent(NULL);
	s.push(root);
	while (!s.empty()) {

		AobfSearchNode *e = s.top();
		s.pop();

		if (e->getType() == NODE_AND) { // this is an AND node
			// add all its pseudo-tree children (if any).
			int var = e->getVar();
			int val = e->getVal();
			m_assignment[var] = val;

			const std::list<AobfSearchNode*>& succ = e->getChildren();
			if (!succ.empty()) { // expanded AND node
				std::list<AobfSearchNode*>::const_iterator ci = succ.begin();
				for (; ci != succ.end(); ++ci) {
					AobfSearchNode *c = (*ci);
					c->setCurrentParent(e);
					if (!c->isSolved()) {
						s.push(c);
					}
				}
			} else { // unexpanded AND node
				e->setFringe(true);
				m_tips.push_back(e);
			}
		} else { // this is an OR node
			assert( m_problem->isMap(e->getVar()) );
			const std::list<AobfSearchNode*>& succ = e->getChildren();
			if (!succ.empty()) { // expanded OR node
				AobfSearchNode* m = e->getBestChild();
				assert(m != NULL);
				m->setCurrentParent(e);
				s.push(m);
			} else { // unexpanded OR node
				e->setFringe(true);
				m_tips.push_back(e);
			}
		}
	} // end while.

	// problem solved if no tips found
	return ( !m_tips.empty() );
}

// sort the tip nodes (ascending, descending, random shuffle)
void AOBF::arrangeTipNodes() {

	std::sort(m_tips.begin(), m_tips.end(), CompNodeHeurAsc());

//	switch (sortType) {
//	case SORT_NODE_HEUR_ASCENDING:
//		std::sort(m_tips.begin(), m_tips.end(), compNodeHeurAsc());
//		break;
//	case SORT_NODE_HEUR_DESCENDING:
//		std::sort(m_tips.begin(), m_tips.end(), compNodeHeurDesc());
//		break;
//	case SORT_NODE_RANDOM:
//		std::random_shuffle(m_tips.begin(), m_tips.end());
//		break;
//	default:
//		// do nothing
//		break;
//	}
}

// select a tip node to expand
AobfSearchNode* AOBF::chooseTipNode() {

	if (m_tips.empty() == false) {
		return *m_tips.begin();
	} else {
		return NULL;
	}
}

// AO* search over the context minimal AND/OR graph restricted to MAP vars
int AOBF::aostar(bool verbose) {

	// init the search space
	int varRoot = m_pseudotree->getRoot()->getVar();
	double h = -ELEM_ENCODE( m_heuristic->getGlobalUB() );
	double w = ELEM_ZERO;
	double* dv = new double[2];
	dv[0] = h;
	dv[1] = w;

	// OR root node (corresponding to dummy variable)
	AobfSearchNode *root = new AobfSearchNodeOR(varRoot, 0);
	root->setHeur(h);
	root->setValue(m_epsilon*h);
	root->setTerminal(false);
	root->setFringe(true);
	root->setHeurCache(dv);

	AobfSearchState state(NODE_OR, "s-2");
	m_space->setRoot(root);
	m_space->add(state, root);
	m_tips.push_back(root);

	int result = SEARCH_SUCCESS;
	while ( !root->isSolved() ) {

		// check for time limit violation
		if (m_timer.timeout()) {
			result = SEARCH_TIMEOUT;
			break;
		}

		// output intermediate results if verbose mode is enabled
		if (verbose) {
			// logging
#ifdef SEARCH_CHECKPOINT_NODE
			size_t expansions = m_space->getAndNodes();
			if (expansions % SEARCH_CHECKPOINT_NODE == 0 &&
				expansions > m_prevNodeCheckpoint) {

				double timestamp = m_timer.elapsed();
				m_prevNodeCheckpoint = expansions;
				std::cout << "[*" << std::setiosflags(std::ios::fixed) << std::setprecision(2)
					<< std::setw(8) << timestamp << "] w " << std::setiosflags(std::ios::fixed) << std::setprecision(4)
					<< std::setw(8) << m_epsilon << " "
					<< std::setw(12) << m_space->getOrNodes() << " "
					<< std::setw(12) << m_space->getAndNodes()
					<< std::endl;
			}
#else
			double elapsed = m_timer.elapsed();
			if (elapsed - m_prevTimeCheckpoint >= SEARCH_CHECKPOINT_TIME) {

				double timestamp = elapsed;
				m_prevTimeCheckpoint = elapsed;
				std::cout << "[*"
					<< std::setw(8) << timestamp << "] w "
					<< std::setw(8) << m_epsilon << " "
					<< std::setw(12) << m_space->getOrNodes() << " "
					<< std::setw(12) << m_space->getAndNodes()
					<< std::endl;
			}
#endif
		}

		// AOstar graph search
		assert(m_tips.size() != 0);
		arrangeTipNodes();
		AobfSearchNode *n = chooseTipNode();
		expandAndRevise(n);
	}

	// get the optimal solution (if solved)
	if (root->isSolved()) {
		m_solutionCost = root->getValue();
	}

	return result;
}


////////////////////////////////////////////////////////////////////////////////
// solve the conditioned summation problem (SUM) via AOC search
double AOBF::_ao(int var, const std::vector<int>& assignment) {

	// a temporary buffer
	double result = ELEM_NAN;

	// reinit the search space
	ao_space->root = NULL;
	m_assignment = assignment;

	// init root of the search space
	_initAoSearch(var);

	AobbSearchNode* n = _nextLeaf();
	while (n) { // throws timeout
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
void AOBF::_initAoSearch(int var) {

	// create root OR node
	AobbSearchNode* node = new AobbSearchNodeOR(NULL, var, 0);
	ao_space->root = node;
	_heuristicOR(node);

	assert(ao_stack.empty());
	ao_stack.push(node);
}

// DONE: radu
void AOBF::_heuristicOR(AobbSearchNode* n) {

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
bool AOBF::_doProcess(AobbSearchNode* node) {
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
void AOBF::_addCacheContext(AobbSearchNode* node, const set<int>& ctxt) const {

	context_t sig;
	sig.reserve(ctxt.size());
	for (set<int>::const_iterator itC = ctxt.begin(); itC != ctxt.end(); ++itC) {
		sig.push_back(m_assignment[*itC]);
	}

	node->setCacheContext(sig);
}

// DONE: radu
bool AOBF::_doCaching(AobbSearchNode* node) {

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

//				// will throw int(UNKNOWN) if not found
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
bool AOBF::_doPruning(AobbSearchNode* node) {

	assert(node);

	if (_canBePruned(node)) {
		node->setLeaf();
		node->setPruned();
		node->setValue(ELEM_ZERO); // dead end

		return true;
	}

	return false; // default false
} // AOBB::doPruning

void AOBF::_getPathAssignment(AobbSearchNode* node, std::vector<int>& vars) {

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

bool AOBF::_canBePruned(AobbSearchNode* n) {

	// propagate current assignment (AND nodes only)
	if (n->getType() == NODE_AND) {

		if (m_options->propagation == CP_FULL) {

			// propagate the assignment
			std::vector<int> vars;
			_getPathAssignment(n, vars);

			bool conflict = propagate(vars);
			if (conflict) {
				++m_deadendsCP;
				return true;
			}

		} else if (m_options->propagation == CP_UNIT) {

//			// Propagate the assignment
//			m_zchaff.dlevel()++;
//			if (!lookahead(n->getVar(), n->getVal(), n->changes(),
//					m_pseudotree->getNode(n->getVar())->getSubprobVars())) {
//				// found inconsistency
//				return true; // 'true' if conflict
//			}
		}
	}

	// default, no pruning possible
	return false;

} // BBBT::canBePruned

// DONE: radu
AobbSearchNode* AOBF::_nextLeaf() {

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
		if (_doPruning(node)) { // pruning?
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
bool AOBF::_generateChildrenAND(AobbSearchNode* n, vector<AobbSearchNode*>& chi) {

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
bool AOBF::_generateChildrenOR(AobbSearchNode* n, vector<AobbSearchNode*>& chi) {

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
AobbSearchNode* AOBF::_nextNode() {
	AobbSearchNode* n = NULL;
	if (!n && ao_stack.size()) {
		n = ao_stack.top();
		ao_stack.pop();
	}
	return n;
}

// DONE: radu
bool AOBF::_doExpand(AobbSearchNode* n) {
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

