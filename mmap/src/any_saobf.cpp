/*
 * any_saobf.cpp
 *
 *  Created on: 29 Sep 2015
 *      Author: radu
 */

#include "any_saobf.h"
#include "mini_bucket_heur.h"
#include "mini_bucket_mm_heur.h"
#include "mini_bucket_jglp_heur.h"

// initialize
void AnySAOBF::init() {

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
int AnySAOBF::solve() {

	// check if solved during initialization
	if (m_solved) {
		std::cout << "--------- Solved during initialization ---------" << std::endl;
		std::cout << "Problem name:        " << m_problem->getName() << std::endl;
		std::cout << "Status:              " << solver_status[0] << std::endl;
		std::cout << "OR nodes:            " << 0 << std::endl;
		std::cout << "AND nodes:           " << 0 << std::endl;
		std::cout << "OR nodes (MAP):      " << 0 << std::endl;
		std::cout << "AND nodes (MAP):     " << 0 << std::endl;
		std::cout << "OUT expansions:      " << 0 << std::endl;
		std::cout << "IN expansions:       " << 0 << std::endl;
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
	m_spaceBFS.reset(new AobfSearchSpace2(m_problem->getN()));
	m_solutionCost = ELEM_NAN;
	m_assignment.resize(m_problem->getN(), UNKNOWN);

	try {

		// init the AO search space
		m_spaceDFS.reset(new AobbSearchSpace());
		m_propagator.reset(new BoundPropagator(m_problem.get(), m_spaceDFS.get(), m_pseudotree.get(), m_options));
		if (m_options->propagation == CP_UNIT) {
			m_propagator->setSatSolver(&m_zchaff);
		} else {
			m_propagator->setSatSolver(NULL);
		}
		m_spaceDFS->cache = new AobbCacheTable(m_problem->getN());

		// AO* search
		res = aostar(m_options->verbose);

		m_solved = true;
	} catch (std::bad_alloc& ba) {
		delete[] EmergencyMemory;
		EmergencyMemory = NULL;
		std::cout << "OUT OF MEMORY";
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
	std::cout << "OR nodes:            " << m_spaceDFS->nodesOR + m_spaceBFS->getOrNodes() << std::endl;
	std::cout << "AND nodes:           " << m_spaceDFS->nodesAND + m_spaceBFS->getAndNodes() << std::endl;
	std::cout << "OR nodes (MAP):      " << m_spaceBFS->getOrNodes() << std::endl;
	std::cout << "AND nodes (MAP):     " << m_spaceBFS->getAndNodes() << std::endl;
	std::cout << "Cache hits (MAP):    " << m_numCacheHits << std::endl;
	std::cout << "SUM evaluations:     " << m_numSumEvals << std::endl;
	std::cout << "OUT expansions:      " << m_numExpansionsOut << std::endl;
	std::cout << "IN expansions:       " << m_numExpansionsIn << std::endl;
	std::cout << "Time elapsed:        " << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
	std::cout << "Preprocessing:       " << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
//	std::cout << "Solution:            " << m_solutionCost << " (" << ELEM_ENCODE(m_solutionCost) << ")" << std::endl;
	std::cout << "Solution:            " << 1/ELEM_DECODE(m_solutionCost) << " (" << -m_solutionCost  << ")" << std::endl;
	std::cout << "-------------------------------" << std::endl;

	// clean up
	if (EmergencyMemory) delete[] EmergencyMemory;

	return res;
}

void AnySAOBF::getPathAssignment(AobfSearchNode* node, std::vector<int>& vars) {

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

// Revise the value of node 'n' based on its children values; log arithmetic
bool AnySAOBF::revise(AobfSearchNode *node) {

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

// revise the f-value of a node (without 'feasible' marking)
bool AnySAOBF::reviseF(AobfSearchNode *node) {

	bool result = true;
	if (node->getType() == NODE_AND) { // AND node

		if (node->isTerminal()) {
			// terminal AND node
			node->setUpperBound(ELEM_ZERO); // feasible value
			node->setFeasible(true);
			node->setSolved(true);
			node->setFringe(false);

			result = true;

		} else {

			double oldFvalue = node->getUpperBound(), fval = 0;
			bool feasible = true;

			std::list<AobfSearchNode*>& succ = node->getChildren();
			std::list<AobfSearchNode*>::iterator li = succ.begin();
			for (; li != succ.end(); ++li) {
				AobfSearchNode *ch = (*li);
				bool sol = ch->isSolved();
				bool fea = ch->isFeasible();
				assert( ch->getType() == NODE_OR );
				feasible &= (fea | sol);

				if (fea || sol) {
					fval += ch->getUpperBound();
				}
			}

			// the node is feasible, if all its OR children
			// are either 'feasible' or 'solved'
			if (feasible) {
				node->setUpperBound(fval);
				node->setFeasible(true);
				node->setFringe(false);
			}

			result = ( feasible || (fval != oldFvalue) );
			//result = true;
		}
	} else { // OR node

		double oldFvalue = node->getUpperBound(), fval = ELEM_INF;
		bool feasible = false;
		std::list<AobfSearchNode*>& succ = node->getChildren();
		std::list<AobfSearchNode*>::iterator li = succ.begin();
		for (; li != succ.end(); ++li) {
			AobfSearchNode *ch = (*li);
			assert( ch->getType() == NODE_AND );
			bool fea = ch->isFeasible();
			feasible |= fea;
			int val = ch->getVal();
			double w = node->getWeight(val);
			double f = ( fea ? (w + ch->getUpperBound()) : ELEM_INF );

			fval = std::min(fval, f);
		}

		// the node is feasible if at least one AND child is feasible or solved
		node->setUpperBound(fval);
		if (feasible) {
			node->setFeasible(true);
			node->setFringe(false);
		}

		result = ( feasible || (fval != oldFvalue) );
	}

	return result;
}

// returns a string of the context instantiation
std::string AnySAOBF::context(int node_type, const std::set<int>& C) {

	std::ostringstream oss;
	if ( C.empty() || node_type == NODE_AND) { // root and dummy node
		oss << "s" << m_globalSearchIndex++;
	} else {
		for (std::set<int>::const_iterator si = C.begin(); si != C.end(); ++si) {
			int cvar = (*si);
			assert(m_assignment[cvar] != UNKNOWN);
			oss << "x" << cvar << "=" << m_assignment[cvar] << ";";
		}
	}

	return oss.str();
}

// expand a search node by generating its successors
bool AnySAOBF::expand(AobfSearchNode *node) {

	// safety checks
	assert(node);
	assert(node->getType() == NODE_AND || node->getType() == NODE_OR);

	// set the path assignment
	std::vector<std::pair<int, int> >& path = node->getPath();
	for (size_t i = 0; i < path.size(); ++i) {
		int var = path[i].first;
		int val = path[i].second;
		m_assignment[var] = val;
	}

	bool noChildren = true;
	if (node->getType() == NODE_AND) { // expand an AND node

		int var = node->getVar();
		int depth = node->getDepth();
		PseudotreeNode* ptNode = m_pseudotree->getNode(var);
		const std::vector<PseudotreeNode*>& children = ptNode->getChildren();
		std::vector<PseudotreeNode*>::const_iterator it = children.begin();
		for (; it != children.end(); ++it) {
			int vChild = (*it)->getVar();
			PseudotreeNode* ptChild = m_pseudotree->getNode(vChild);
			std::string str = context(NODE_OR, ptChild->getFullContext());
			AobfSearchState state(NODE_OR, str);

			if (m_spaceBFS->find(vChild, state) == true) {
				AobfSearchNodeOR* c = (AobfSearchNodeOR*) m_spaceBFS->get(vChild, state);
				node->addChild(c);
				c->addParent(node);
				++m_numCacheHits;
			} else {
				AobfSearchNodeOR* c = new AobfSearchNodeOR(vChild, depth + 1);
				c->addParent(node);
				node->addChild(c);
				c->getPath() = node->getPath(); // inherit the path from parent

				if (m_problem->isMap(vChild)) {
					double h = heuristicOR(c); // log space
					c->setHeur( h );
					c->setValue( m_epsilon*h );
					c->setUpperBound(ELEM_INF);
					c->setFringe(true);
				} else { // conditioned SUM problem rooted at an OR node
					assert(m_problem->isSum(vChild));
					std::vector<int> backup = m_assignment;
					double qval = _sum(vChild, m_assignment);
					c->setHeur(ELEM_NAN);
					c->setValue( -ELEM_ENCODE(qval) );
					c->setUpperBound( -ELEM_ENCODE(qval) );
					c->setTerminal(true);
					c->setSolved(true);
					c->setFringe(false);
					c->setExpanded(true);
					c->setFeasible(true);
					m_assignment = backup;
					++m_numSumEvals;
				}

				// add the OR child to the search space
				m_spaceBFS->add(state, c);
				if (m_problem->isMap(c->getVar())) {
					m_tipsOut.push(c);  // by default add it to the outside tips
										// ignore SUM variables
				}
			}

			noChildren = false;
		}

		// set the flags: expanded, fringe, terminal
		node->setExpanded(true);
		node->setFringe(false);
		node->setTerminal( children.empty() );
		m_spaceBFS->incNodesExpanded(NODE_AND);

	} else { // expand an OR node
		int var = node->getVar();
		int depth = node->getDepth();
		PseudotreeNode* ptNode = m_pseudotree->getNode(var);
		assert( m_problem->isMap(var) == true );
		assert( node->getHeurCache() != NULL );
		for (int val = 0; val < m_problem->getDomainSize(var); ++val) {

			// create the state corresponding to the AND child
			std::string str = context(NODE_AND, ptNode->getFullContext());
			AobfSearchState state(NODE_AND, str);
			AobfSearchNodeAND* c = new AobfSearchNodeAND(var, val, depth + 1);
			node->addChild(c);
			c->addParent(node);
			double *heurCache = node->getHeurCache();
			double h = heurCache[2*val]; // includes the weight
			double w = heurCache[2*val+1];
			c->setHeur(h - w);
			c->setValue(m_epsilon*(h - w));
			c->getPath() = node->getPath(); // inherit the path from parent
			c->getPath().push_back(std::make_pair(var, val)); // extend the path
			c->setFringe(true);
			m_spaceBFS->add(state, c);
			m_tipsOut.push(c); // by default add it to the outside tips

			noChildren = false;
		}

		// set the flags: expanded, fringe, terminal
		node->setExpanded(true);
		node->setFringe(false);
		node->setTerminal(false);
		m_spaceBFS->incNodesExpanded(NODE_OR);
	}

	return noChildren; // "true" if no children, "false" otherwise
}

// DONE: radu
double AnySAOBF::heuristicOR(AobfSearchNode* n) {

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
void AnySAOBF::expandAndRevise(AobfSearchNode *node) {

	// expand node n: Steps 6 and 7 from Nilsson's
	assert( node->isFringe() && !node->isSolved() );
	expand(node);

#ifdef DEBUG
	std::cout << "Expanded "<< node->toString() << std::endl;
#endif

	update(node); // q-values
	updateF(node); // f-values

	// Find next best partial solution tree: Step 4 from Nilsson's
	m_tips.clear();
	findBestPartialTree();

	assert( m_spaceBFS->getRoot()->isSolved() || (m_tips.size() > 0) );
}

// propagate the q-values to ancestors along marked paths
void AnySAOBF::update(AobfSearchNode* node) {

	// keep track of visited ancestors
	std::multiset<AobfSearchNode*, CompNodeIndexAsc> S;

	// start from the current node (ie, just expanded)
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
				// AND nodes have a single parent in the context-minimal AND/OR graph.

				assert(e->getParents().size() == 1);
				AobfSearchNode* p = e->getParents().front();
				assert(m_problem->isMap(p->getVar()));

				// count how many of e's children are still in S (marked visited)
				size_t index = 0;
				std::list<AobfSearchNode*>::iterator ci = p->getChildren().begin();
				for (; ci != p->getChildren().end(); ++ci) {
					AobfSearchNode *c = (*ci);
					if (c->isVisited()) {
						++index;
					}
				}

				AobfSearchNode *best = p->getBestChild();
				bool found = (best == e);	// check if 'e' is the marked AND child of 'p'
				if (p->isVisited()) {
					// if the parent already in S, decrease the index.
					std::multiset<AobfSearchNode*, CompNodeIndexAsc>::iterator si = S.begin();
					for (; si != S.end(); ++si) {
						if (p == *si) {
							if (p->getIndex() > 0) {
								/*si = */S.erase(si);
								p->decIndex();
								S.insert(p);
							}
							break;
						}
					}
				} else if (found) {
					p->setIndex(index);
					p->setVisited(true);
					S.insert(p);
				}
			} else if (e->getType() == NODE_OR) {

				// OR nodes may have multiple parents in the context-minimal AND/OR graph.
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

					if (p->isVisited()) {
						// if the parent already in S, decrease the index.
						std::multiset<AobfSearchNode*, CompNodeIndexAsc>::iterator si = S.begin();
						for (; si != S.end(); ++si) {
							if (p == *si) {
								if (p->getIndex() > 0) {
									/*si = */S.erase(si);
									p->decIndex();
									S.insert(p);
								}
								break;
							}
						}
					} else {
						// new parent through marked connector
						p->setIndex(index);
						p->setVisited(true);
						S.insert(p);
					}
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
							if (p->getIndex() > 0) {
								/*si = */S.erase(si);
								p->decIndex();
								S.insert(p);
							}
							break;
						}
					}
				}
			}
		} // end if
	}

	assert( S.empty() );
}

// propagate the q-values to all ancestors of the current node
void AnySAOBF::updateAll(AobfSearchNode* node) {

	std::multiset<AobfSearchNode*, CompNodeIndexAsc> S;
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
				// AND nodes have a single parent in the context-minimal AND/OR graph.

				assert(e->getParents().size() == 1);
				AobfSearchNode* p = e->getParents().front();
				assert(m_problem->isMap(p->getVar()));

				// count how many of e's children are still in S (marked visited)
				size_t index = 0;
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
							if (p->getIndex() > 0) {
								/*si = */S.erase(si);
								p->decIndex();
								S.insert(p);
							}
							break;
						}
					}
				} else {
					p->setIndex(index);
					p->setVisited(true);
					S.insert(p);
				}
			} else if (e->getType() == NODE_OR) {

				// OR nodes may have multiple parents in the context-minimal AND/OR graph.
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

					if (p->isVisited()) {
						// if the parent already in S, decrease the index.
						std::multiset<AobfSearchNode*, CompNodeIndexAsc>::iterator si = S.begin();
						for (; si != S.end(); ++si) {
							if (p == *si) {
								if (p->getIndex() > 0) {
									/*si = */S.erase(si);
									p->decIndex();
									S.insert(p);
								}
								break;
							}
						}
					} else {
						// new parent through marked connector
						p->setIndex(index);
						p->setVisited(true);
						S.insert(p);
					}
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
							if (p->getIndex() > 0) {
								/*si = */S.erase(si);
								p->decIndex();
								S.insert(p);
							}
							break;
						}
					}
				}
			}
		} // end if
	}

	assert( S.empty() );

	// compute best partial solution tree: Step 4 from Nilsson's
//	m_tips.clear();
//	findBestPartialTree();
//
//	assert( m_spaceBFS->getRoot()->isSolved() || (m_tips.size() > 0) );
}

// propagate the f-values to all ancestors of the current node
void AnySAOBF::updateF(AobfSearchNode* node) {

	std::multiset<AobfSearchNode*, CompNodeIndexAsc> S;
	S.insert(node);
	node->setVisited(true);

	// iterate over S: Step 10 from Nilsson's
	while ( !S.empty() ) {
		AobfSearchNode *e = *S.begin();
		S.erase(S.begin());
		e->setVisited(false);

		assert( e->getIndex() == 0 );
		bool change = reviseF(e);
#ifdef DEBUG
	std::cout << "Revised "<< e->toString() << std::endl;
#endif

		if (change) {

			// update parents: Step 12 from Nilsson's
			if (e->getType() == NODE_AND) {
				// AND nodes have a single parent in the context-minimal AND/OR graph.

				assert(e->getParents().size() == 1);
				AobfSearchNode* p = e->getParents().front();
				assert(m_problem->isMap(p->getVar()));

				// count how many of e's children are still in S (marked visited)
				size_t index = 0;
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
							if (p->getIndex() > 0) {
								/*si = */S.erase(si);
								p->decIndex();
								S.insert(p);
							}
							break;
						}
					}
				} else {
					p->setIndex(index);
					p->setVisited(true);
					S.insert(p);
				}
			} else if (e->getType() == NODE_OR) {

				// OR nodes may have multiple parents in the context-minimal AND/OR graph.
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

					if (p->isVisited()) {
						// if the parent already in S, decrease the index.
						std::multiset<AobfSearchNode*, CompNodeIndexAsc>::iterator si = S.begin();
						for (; si != S.end(); ++si) {
							if (p == *si) {
								if (p->getIndex() > 0) {
									/*si = */S.erase(si);
									p->decIndex();
									S.insert(p);
								}
								break;
							}
						}
					} else {
						// new parent through marked connector
						p->setIndex(index);
						p->setVisited(true);
						S.insert(p);
					}
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
							if (p->getIndex() > 0) {
								/*si = */S.erase(si);
								p->decIndex();
								S.insert(p);
							}
							break;
						}
					}
				}
			}
		} // end if
	}

	assert( S.empty() );
}

// find the new best partial solution tree
bool AnySAOBF::findBestPartialTree() {

	// get the root of the search space
	AobfSearchNode *root = m_spaceBFS->getRoot();

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
void AnySAOBF::arrangeTipNodes() {

	std::sort(m_tips.begin(), m_tips.end(), CompNodeHeurAsc());
}

// select a tip node to expand: with probability p select from the current PST
// and with probability (1-p) select a tip node outside the current PST
AobfSearchNode* AnySAOBF::chooseTipNode(double threshold) {

	// generate the selection probability p
	double p = rand::probability();
	if (p < threshold) { // tip node from from PST
		arrangeTipNodes();
		if (m_tips.empty() == false) {
			++m_numExpansionsIn;
			return *m_tips.begin();
		} else {
			return NULL;
		}
	} else { // tip node outside the PST

		std::vector<AobfSearchNode*> temp;
		bool found = false;
		AobfSearchNode *tip = NULL;
		while (!found) {
			tip = m_tipsOut.pop();
			if (tip == NULL) {
				break;
			} else {
				if (std::find(m_tips.begin(), m_tips.end(), tip) != m_tips.end()) {
					temp.push_back(tip); // tip node from PST so skip it
					continue;
				} else {
					found = true;
				}
			}
		}

		// push back tip nodes that were ignored
		for (size_t i = 0; i < temp.size(); ++i) {
			m_tipsOut.push(temp[i]);
		}

		if (tip != NULL) {
			++m_numExpansionsOut;
		} else { // couldn't find a tip node outside PST so select from PST
			arrangeTipNodes();
			if (m_tips.empty() == false) {
				++m_numExpansionsIn;
				tip = *m_tips.begin();
			} else {
				tip = NULL;
			}
		}

		return tip;
	}
}

// AO* search over the context minimal AND/OR graph restricted to MAP vars
int AnySAOBF::aostar(bool verbose) {

	// create the root node of the search space
	//int numLayers = 2*(m_pseudotree->getHeight() + 2);

	// init the search space
	int varRoot = m_pseudotree->getRoot()->getVar();
	double h = -ELEM_ENCODE( m_heuristic->getGlobalUB() );
	double w = ELEM_ZERO;
	double* dv = new double[2];
	dv[0] = h;
	dv[1] = w;

	m_upperBound = ELEM_INF; // best solution found so far
	m_lowerBound = h; // initial heuristic value

	// OR root node (corresponding to dummy variable)
	AobfSearchNode *root = new AobfSearchNodeOR(varRoot, 0);
	root->setHeur(h);
	root->setValue(m_epsilon*h);
	root->setTerminal(false);
	root->setFringe(true);
	root->setHeurCache(dv);

	AobfSearchState state(NODE_OR, "s-2");
	//m_spaceBFS->initLayers(numLayers);
	m_spaceBFS->setRoot(root);
	m_spaceBFS->add(state, root);
	m_tips.push_back(root);

	// initialize the tips index
	int sz = 2*m_problem->getN() + 1;
	m_tipsOut.init(sz);

	size_t iterations = 0;
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
				double upbo = root->getUpperBound();
				double lowbo = root->getValue();
				m_prevNodeCheckpoint = expansions;
				std::cout << "[*" << std::setw(8) << timestamp << "] u "
					<< std::setw(8)  << -1 << " "
					<< std::setw(12) << m_spaceDFS->nodesOR + m_spaceBFS->getOrNodes() << " "
					<< std::setw(12) << m_spaceDFS->nodesAND + m_spaceBFS->getAndNodes() << " "
					<< std::setw(12) << m_spaceBFS->getOrNodes() << " "
					<< std::setw(12) << m_spaceBFS->getAndNodes() << " "
					<< std::setw(12) << upbo << " (" << 1/ELEM_DECODE(upbo) << ") "
					<< std::setw(12) << lowbo << " (" << 1/ELEM_DECODE(lowbo) << ")"
					<< std::endl;
			}
#else
			double elapsed = m_timer.elapsed();
			if (elapsed - m_prevTimeCheckpoint >= SEARCH_CHECKPOINT_TIME) {

				double timestamp = elapsed;
				double upbo = root->getUpperBound();
				double lowbo = root->getValue();
				m_prevTimeCheckpoint = elapsed;
				std::cout << "[*" << std::setw(8) << timestamp << "] u "
					<< std::setw(8)  << -1 << " "
					<< std::setw(12) << m_spaceDFS->nodesOR + m_spaceBFS->getOrNodes() << " "
					<< std::setw(12) << m_spaceDFS->nodesAND + m_spaceBFS->getAndNodes() << " "
					<< std::setw(12) << m_spaceBFS->getOrNodes() << " "
					<< std::setw(12) << m_spaceBFS->getAndNodes() << " "
					<< std::setw(12) << upbo << " (" << 1/ELEM_DECODE(upbo) << ") "
					<< std::setw(12) << lowbo << " (" << 1/ELEM_DECODE(lowbo) << ")"
					<< std::endl;
			}
#endif
		}

		// AOstar graph search
		double threshold = (iterations == 0) ? 1.0 : m_options->probability;
		AobfSearchNode *n = chooseTipNode(threshold); // first iteration: select from PST
		if (n != NULL) {
			expandAndRevise(n);
			m_tipsOut.erase(n); // remove it from the index

			// Check if we found a new feasible solution
			if (root->getUpperBound() < m_upperBound) {
				m_upperBound = root->getUpperBound();
				double lowbo = root->getValue();
				std::cout << "["
					<< std::setw(9) << m_timer.elapsed() << "] w "
					<< std::setw(8) << m_epsilon << " "
					<< std::setw(12) << m_spaceDFS->nodesOR + m_spaceBFS->getOrNodes() << " "
					<< std::setw(12) << m_spaceDFS->nodesAND + m_spaceBFS->getAndNodes() << " "
					<< std::setw(12) << m_spaceBFS->getOrNodes() << " "
					<< std::setw(12) << m_spaceBFS->getAndNodes() << " "
					<< std::setw(12) << m_upperBound << " (" << 1/ELEM_DECODE(m_upperBound) << ") "
					<< std::setw(12) << lowbo << " (" << 1/ELEM_DECODE(lowbo) << ")"
					<< std::endl;
			}

			++iterations;
		} else {
			break; // search done
		}
	}

	// get the optimal solution (if solved)
	if (root->isSolved()) {
		m_solutionCost = root->getValue();
	} else { // or the best solution found so far if timeout
		m_solutionCost = root->getUpperBound();
	}

	// final output
	std::cout << "[" << std::setw(9) << m_timer.elapsed() << "] u "
		<< std::setw(8)  << -1 << " "
		<< std::setw(12) << m_spaceDFS->nodesOR + m_spaceBFS->getOrNodes() << " "
		<< std::setw(12) << m_spaceDFS->nodesAND + m_spaceBFS->getAndNodes() << " "
		<< std::setw(12) << m_spaceBFS->getOrNodes() << " "
		<< std::setw(12) << m_spaceBFS->getAndNodes() << " "
		<< std::setw(12) << root->getUpperBound() << " (" << 1/ELEM_DECODE(root->getUpperBound()) << ") "
		<< std::setw(12) << root->getValue() << " (" << 1/ELEM_DECODE(root->getValue()) << ")"
		<< std::endl;

	return result;
}

void AnySAOBF::check() {
	// check the tip nodes outside

}

AobbSearchNode* AnySAOBF::_nextNode() {
	AobbSearchNode* n = NULL;

	// check first if we're in a SUM problem
	if (!n && m_stack.size()) {
		n = m_stack.top();
		m_stack.pop();
		return n;
	}

	while (m_stacks.size()) {
		MyStack* st = m_stacks.front();
		if (st->getChildren()) {
			m_stacks.push(st);
			m_stackCount = 0;
		} else if (st->empty()) {
			if (st->getParent())
				st->getParent()->delChild();
			delete st;
			m_stackCount = 0;
		} else if (m_stackLimit && m_stackCount++ == m_stackLimit) {
			m_stacks.push(st);
			m_stackCount = 0;
		} else {
			n = st->top();
			st->pop();
			break; // while loop
		}
		m_stacks.pop();
	}
	return n;
}

bool AnySAOBF::_doExpand(AobbSearchNode* n) {
	assert(n);
	m_expand.clear();

	MyStack* stack = m_stacks.front();
	if (n->getType() == NODE_AND) {  // AND node
		if (_generateChildrenAND(n, m_expand))
			return true; // no children

		if (m_problem->isSum(m_expand.at(0)->getVar())) {
			for (vector<AobbSearchNode*>::reverse_iterator it = m_expand.rbegin();
					it != m_expand.rend(); ++it) {
				m_stack.push(*it);
			}
		} else { // MAP variables
//			std::cout << "[DEBUG] expanded " << n->toString() << endl;
			if (m_expand.size() == 1) {  // no decomposition
				stack->push(m_expand.at(0));
			} else {  // decomposition, split stacks
				// reverse iterator needed since new stacks are put in queue (and not depth-first stack)
				for (vector<AobbSearchNode*>::reverse_iterator it = m_expand.rbegin();
						it != m_expand.rend(); ++it) {
					MyStack* s = new MyStack(stack);
					stack->addChild();
					s->push(*it);
					m_stacks.push(s);
				}
			}
		}
	} else {  // OR node
		if (_generateChildrenOR(n, m_expand))
			return true; // no children

		if (m_problem->isSum(m_expand.at(0)->getVar())) {
			for (vector<AobbSearchNode*>::reverse_iterator it = m_expand.rbegin();
					it != m_expand.rend(); ++it) {
				m_stack.push(*it);
			}  // for loop
		} else { // MAP variables
//			std::cout << "[DEBUG] expanded " << n->toString() << endl;
			for (vector<AobbSearchNode*>::reverse_iterator it = m_expand.rbegin();
					it != m_expand.rend(); ++it) {
				stack->push(*it);
			}  // for loop
		}
	}
	// if over node type

	return false; // default false
}

// Solve the conditioned summation problem (SUM) via DFS search
double AnySAOBF::_sum(int var, const std::vector<int>& assignment) {

	// safety checks
	assert(m_problem->isSum(var));

	// conditional likelihood value
	double result = ELEM_NAN;

	// reinit the search space
	m_spaceDFS->root = NULL;
	m_assignment = assignment;

	// init root of the search space
	AobbSearchNode* r = _initSearchSpace(var);
	assert(m_stack.empty());
	m_stack.push(r);

	AobbSearchNode* n = _nextLeaf();
	while (n) { // throws timeout
		m_propagator->propagate(n, m_timer, m_currentDomains, false); // true = report solutions
		n = _nextLeaf();
	}

	// Conditioned SUM problem solved
	result = m_spaceDFS->root->getValue();
	delete m_spaceDFS->root;
	m_spaceDFS->root = NULL;

	return result;
}

// DFS probe below the node (returning an upper bound). Stop at the first bound.
double AnySAOBF::_probe(int var, const std::vector<int>& assignment) {

	// Upper bound
	double upperBound = ELEM_INF;

	// Re-init the search space
	m_spaceDFS->root = NULL;
	m_assignment = assignment;
	AobbSearchNode* first = _initSearchSpace(var);
	if (first) {
		m_rootStack = new MyStack(NULL);
		m_rootStack->push(first);
		m_stacks.push(m_rootStack);
		m_stackCount = 0;
	}

	//std::cout << "[DEBUG] first node is " << first->toString() << std::endl;

	AobbSearchNode* n = _nextLeaf();
	while (n) { // throws timeout
		m_propagator->propagate(n, m_timer, m_currentDomains, false); // true = report solutions
		double lb = m_propagator->getLowerBound(); // best solution found so far
		if (ISNAN(lb) == false) {
			upperBound = -ELEM_ENCODE(lb);
			break; // do something more sophisticated here
		}

		n = _nextLeaf();
	}

	// clean-up the stacks
	delete m_spaceDFS->root;
	m_spaceDFS->root = NULL;

	while (m_stacks.size()) {
		delete m_stacks.front();
		m_stacks.pop();
	}

	return upperBound;
}

// DONE: radu
AobbSearchNode* AnySAOBF::_initSearchSpace(int var) {
	assert(m_spaceDFS->root == NULL);

	// create root OR node
	AobbSearchNode* node = new AobbSearchNodeOR(NULL, var, 0);
	m_spaceDFS->root = node;
	_heuristicOR(node);

	return node;
}

// DONE: radu
double AnySAOBF::_heuristicOR(AobbSearchNode* n) {

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
bool AnySAOBF::_doProcess(AobbSearchNode* node) {
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
void AnySAOBF::_addCacheContext(AobbSearchNode* node, const set<int>& ctxt) const {

	context_t sig;
	sig.reserve(ctxt.size());
	for (set<int>::const_iterator itC = ctxt.begin(); itC != ctxt.end(); ++itC) {
		sig.push_back(m_assignment[*itC]);
	}

	node->setCacheContext(sig);
}

// DONE: radu
bool AnySAOBF::_doCaching(AobbSearchNode* node) {

	assert(node);
	int var = node->getVar();
	PseudotreeNode* ptnode = m_pseudotree->getNode(var);

	if (node->getType() == NODE_AND) { // AND node -> reset associated adaptive cache tables

		// if adaptive caching then purge corresponding caches
	    const list<int>& resetList = ptnode->getCacheReset();
	    for (list<int>::const_iterator it = resetList.begin();
	    		it != resetList.end(); ++it) {
	    	m_spaceDFS->cache->reset(*it);
	    }

	} else { // OR node, try actual caching

		if (!ptnode->getParent())
			return false;

		if (ptnode->getFullContext().size() <=
			ptnode->getParent()->getFullContext().size()) {

			// add cache context information
			//addCacheContext(node, ptnode->getFullContext());
			_addCacheContext(node, ptnode->getCacheContext());

			// try to get value from cache
			try {
				// will throw int(UNKNOWN) if not found
#ifndef NO_ASSIGNMENT
				pair<double,vector<val_t> > entry = m_spaceDFS->cache->read(var, node->getCacheContext());
				node->setValue( entry.first ); // set value
				node->setOptAssig( entry.second ); // set assignment
#else
				double entry = m_spaceDFS->cache->read(var, node->getCacheContext());
				node->setValue(entry); // set value
#endif
				node->setLeaf(); // mark as leaf
				//++m_cacheHits;

				return true;
			} catch (...) { // cache lookup failed
				node->setCachable(); // mark for caching later
			}
		}
	} // if on node type

	return false; // default, no caching applied

} // AOBB::doCaching

// DONE: radu
bool AnySAOBF::_doPruning(AobbSearchNode* node) {

	assert(node);

	return false; // default false
} // AOBB::doPruning

// DONE: radu
AobbSearchNode* AnySAOBF::_nextLeaf() {

	AobbSearchNode* node = this->_nextNode();
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
				double upbo = m_spaceBFS->getRoot()->getUpperBound();
				double lowbo = m_spaceBFS->getRoot()->getValue();
				m_prevNodeCheckpoint = expansions;
				std::cout << "[*"
					<< std::setw(8) << timestamp << "] u "
					<< std::setw(8) << -1 << " "
					<< std::setw(12) << m_spaceDFS->nodesOR + m_spaceBFS->getOrNodes() << " "
					<< std::setw(12) << m_spaceDFS->nodesAND + m_spaceBFS->getAndNodes() << " "
					<< std::setw(12) << m_spaceBFS->getOrNodes() << " "
					<< std::setw(12) << m_spaceBFS->getAndNodes() << " "
					<< std::setw(12) << upbo << " (" << 1/ELEM_DECODE(upbo) << ") "
					<< std::setw(12) << lowbo << " (" << 1/ELEM_DECODE(lowbo) << ")"
					<< std::endl;
			}
#else
			double elapsed = m_timer.elapsed();
			if (elapsed - m_prevTimeCheckpoint >= SEARCH_CHECKPOINT_TIME) {

				double timestamp = elapsed;
				double upbo = m_spaceBFS->getRoot()->getUpperBound();
				double lowbo = m_spaceBFS->getRoot()->getValue();
				m_prevTimeCheckpoint = elapsed;
				std::cout << "[*"
					<< std::setw(8) << timestamp << "] u "
					<< std::setw(8) << -1 << " "
					<< std::setw(12) << m_spaceDFS->nodesOR + m_spaceBFS->getOrNodes() << " "
					<< std::setw(12) << m_spaceDFS->nodesAND + m_spaceBFS->getAndNodes() << " "
					<< std::setw(12) << m_spaceBFS->getOrNodes() << " "
					<< std::setw(12) << m_spaceBFS->getAndNodes() << " "
					<< std::setw(12) << upbo << " (" << 1/ELEM_DECODE(upbo) << ") "
					<< std::setw(12) << lowbo << " (" << 1/ELEM_DECODE(lowbo) << ")"
					<< std::endl;
			}
#endif
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

void AnySAOBF::_getPathAssignment(AobbSearchNode* node, std::vector<int>& vars) {

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
bool AnySAOBF::_generateChildrenAND(AobbSearchNode* n, vector<AobbSearchNode*>& chi) {

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
	m_spaceDFS->nodesAND += 1;
	if (m_problem->isMap(var))
		m_spaceDFS->nodesANDmap += 1;

	// create new OR children (going in reverse due to reversal on stack)
	for (vector<PseudotreeNode*>::const_reverse_iterator it =
			ptnode->getChildren().rbegin(); it != ptnode->getChildren().rend();
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

	// order subproblems in decreasing order of their heuristic - largest UB first
	// (use reverse iterator due to stack reversal)
	if (m_problem->isMap(n->getVar())) {
		std::sort(chi.begin(), chi.end(), AobbSearchNode::heurGreater);
	}

	n->addChildren(chi);

	return false; // default

} // AOBB::generateChildrenAND

// DONE: radu
bool AnySAOBF::_generateChildrenOR(AobbSearchNode* n, vector<AobbSearchNode*>& chi) {

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
	m_spaceDFS->nodesOR += 1;
	if (m_problem->isMap(var))
		m_spaceDFS->nodesORmap += 1;

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

