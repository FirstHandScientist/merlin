/*
 * any_anaobf.cpp
 *
 *  Created on: 9 Oct 2015
 *      Author: radu
 */

#include "any_aaobf.h"

// AO* search over the context minimal AND/OR graph restricted to MAP vars
int AnyAAOBF::aostar(bool verbose) {

	// create the root node of the search space
	int numLayers = 2*(m_pseudotree->getHeight() + 2);

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
	m_space->initLayers(numLayers);
	m_space->setRoot(root);
	m_space->add(state, root);
	m_tipsF.push_back(root); // start from the FST first
	m_assignmentF = m_assignment;

	bool doRepair = false;
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
				std::cout << "[*"
					<< std::setw(8) << timestamp << "] w "
					<< std::setw(8) << m_epsilon << " "
					<< std::setw(12) << m_space->getOrNodes() << " "
					<< std::setw(12) << m_space->getAndNodes() << " "
					<< std::setw(12) << m_upperBound << " (" << 1/ELEM_DECODE(m_upperBound) << ") "
					<< std::setw(12) << m_lowerBound << " (" << 1/ELEM_DECODE(m_lowerBound) << ")"
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
					<< std::setw(12) << m_space->getAndNodes() << " "
					<< std::setw(12) << m_numExpansionsFST << " "
					<< std::setw(12) << m_numExpansionsBST << " "
					<< std::setw(12) << m_upperBound << " (" << 1/ELEM_DECODE(m_upperBound) << ") "
					<< std::setw(12) << m_lowerBound << " (" << 1/ELEM_DECODE(m_lowerBound) << ")"
					<< std::endl;
			}
#endif
		}

		// AOstar graph search
		AobfSearchNode* n = NULL;
		if ( !m_tipsF.empty() ) { // current FST is a partial tree
			m_assignment = m_assignmentF;	// set the current assignment
			arrangeTipNodesF(); 			// select a tip node from FST
			n = chooseTipNodeF();
			expandAndReviseF(n);			// expand and propagate f-values (upper bounds)
			findFeasiblePartialTree(); 		// find the next FST (updates assignmnent)
			doRepair = true;				// set the flag for full repair
			m_numExpansionsFST++;
		} else { // current FST is a solution tree or the current upper bound hasn't improved
			// check if we have a new solution
			if (root->getUpperBound() < m_upperBound) {
				m_upperBound = root->getUpperBound();
				std::cout << "["
					<< std::setw(9) << m_timer.elapsed() << "] w "
					<< std::setw(8) << m_epsilon << " "
					<< std::setw(12) << m_space->getOrNodes() << " "
					<< std::setw(12) << m_space->getAndNodes() << " "
					<< std::setw(12) << m_numExpansionsFST << " "
					<< std::setw(12) << m_numExpansionsBST << " "
					<< std::setw(12) << m_upperBound << " (" << 1/ELEM_DECODE(m_upperBound) << ") "
					<< std::setw(12) << m_lowerBound << " (" << 1/ELEM_DECODE(m_lowerBound) << ")"
					<< std::endl;
				// write solution to output file
#ifdef UAI_COMPETITION
				writeSolution();
#endif
			}

			// repair the q-values, f-values, and markings (upon finding a solution)
			if (doRepair) {
				repair();				// propagate f-values, q-values
				doRepair = false;		//
				findBestPartialTree(); 	// find the BST (and set m_assignment)
				m_numRepairs++;
			}

			if (m_tips.empty()) {
				break; // solved, so we're done!
			}

			// expand from BST, hoping that we move to a better FST
			// by updating the f-values, q-values and feasible marks
			arrangeTipNodes();			// select a tip node from BST
			n = chooseTipNode();
			expandAndRevise(n);	 		// expand and update q-values, and finds new BST
			findFeasiblePartialTree();  // if still empty, expand from BST,
										// else expand from FST
			m_numExpansionsBST++;
		}
	} // end while (search)

	// Collect final solution (upper and lower bound)
	m_solutionCost = m_upperBound = root->getUpperBound();
	m_lowerBound = root->getValue();
	double timestamp = m_timer.elapsed();
	std::cout << "["
		<< std::setw(9) << timestamp << "] w "
		<< std::setw(8) << m_epsilon << " "
		<< std::setw(12) << m_space->getOrNodes() << " "
		<< std::setw(12) << m_space->getAndNodes() << " "
		<< std::setw(12) << m_numExpansionsFST << " "
		<< std::setw(12) << m_numExpansionsBST << " "
		<< std::setw(12) << m_upperBound << " (" << 1/ELEM_DECODE(m_upperBound) << ") "
		<< std::setw(12) << m_lowerBound << " (" << 1/ELEM_DECODE(m_lowerBound) << ")"
		<< std::endl;

	return result;
}

// solve the problem
int AnyAAOBF::solve() {

#ifndef UAI_COMPETITION
	// check if solved during initialization
	if (m_solved) {
		std::cout << "--------- Solved during initialization ---------" << std::endl;
		std::cout << "Problem name:        " << m_problem->getName() << std::endl;
		std::cout << "Status:              " << solver_status[0] << std::endl;
		std::cout << "OR nodes:            " << 0 << std::endl;
		std::cout << "AND nodes:           " << 0 << std::endl;
		std::cout << "OR nodes (MAP):      " << 0 << std::endl;
		std::cout << "AND nodes (MAP):     " << 0 << std::endl;
		std::cout << "FST expansions:      " << 0 << std::endl;
		std::cout << "BST expansions:      " << 0 << std::endl;
		std::cout << "Repairs:             " << 0 << std::endl;
		std::cout << "Deadends (CP):       " << 0 << std::endl;
		std::cout << "SUM evaluations:     " << 0 << std::endl;
		std::cout << "Time elapsed:        " << (m_tmLoad + m_tmHeuristic + m_tmSolve) << " seconds" << std::endl;
		std::cout << "Preprocessing:       " << (m_tmLoad + m_tmHeuristic) << " seconds" << std::endl;
		std::cout << "Solution:            " << m_solutionCost << " (" << ELEM_ENCODE(m_solutionCost) << ")" << std::endl;
		std::cout << "-------------------------------" << std::endl;

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
	m_space.reset(new AobfSearchSpace2(m_problem->getN()));
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
	std::cout << "OR nodes:            " << ao_space->nodesOR + m_space->getOrNodes() << std::endl;
	std::cout << "AND nodes:           " << ao_space->nodesAND + m_space->getAndNodes() << std::endl;
	std::cout << "OR nodes (MAP):      " << m_space->getOrNodes() << std::endl;
	std::cout << "AND nodes (MAP):     " << m_space->getAndNodes() << std::endl;
	std::cout << "Cache hits (MAP):    " << m_numCacheHits << std::endl;
	std::cout << "FST expansions:      " << m_numExpansionsFST << std::endl;
	std::cout << "BST expansions:      " << m_numExpansionsBST << std::endl;
	std::cout << "Repairs:             " << m_numRepairs << std::endl;
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


// expand a search node by generating its successors
bool AnyAAOBF::expand(AobfSearchNode *node) {

	// safety checks
	assert(node);
	assert(node->getType() == NODE_AND || node->getType() == NODE_OR);

	//bool noChildren = true;
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

			//std::cout << vChild << " OR context: " << str.c_str() << std::endl;
			//std::cout.flush();

			if (m_space->find(vChild, state) == true) {
				AobfSearchNodeOR* c = (AobfSearchNodeOR*) m_space->get(vChild, state);
				node->addChild(c);
				c->addParent(node);
				++m_numCacheHits;
			} else {
				AobfSearchNodeOR* c = new AobfSearchNodeOR(vChild, depth + 1);
				c->addParent(node);
				node->addChild(c);

				if (m_problem->isMap(vChild)) {
					double h = heuristicOR(c); // log space
					c->setHeur( h );
					c->setValue( m_epsilon*h );
					c->setUpperBound(ELEM_INF);
				} else { // conditioned SUM problem rooted at an OR node
					assert(m_problem->isSum(vChild));
					std::vector<int> backup = m_assignment;
					double qval = _ao(vChild, m_assignment);
					c->setHeur(ELEM_NAN);
					c->setValue( -ELEM_ENCODE(qval) );
					c->setUpperBound( -ELEM_ENCODE(qval) );
					c->setTerminal(true);
					c->setSolved(true);
					c->setFeasible(true);
					c->setFringe(false);
					c->setExpanded(true);
					m_assignment = backup;
					++m_numSumEvals;
				}

				// add the OR child to the search space
				m_space->add(state, c);
			}

			//noChildren = false;
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
//			if (isDeadend(node, var, val)) {
//				++m_deadendsCP;
//				continue;
//			}

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
			m_space->add(state, c);

			//noChildren = false;
		}

		// set the flags: expanded, fringe, terminal
		node->setExpanded(true);
		node->setFringe(false);
		node->setTerminal(false);
		m_space->incNodesExpanded(NODE_OR);
	}

	//return noChildren; // "true" if no children, "false" otherwise
	return true;
}

// Expand and revise (from the FST)
void AnyAAOBF::expandAndReviseF(AobfSearchNode* node) {
	// expand node n
	assert( node->isFringe() );
	assert( node->getChildren().empty() );

	node->setFringe(false);
	expand(node);
	node->setExpanded(true);
	m_space->incNodesExpanded(node->getType());
	if (node->getType() == NODE_OR) {
		markFeasibleChild(node);
	}

#ifdef DEBUG
//	if (node->getType() == NODE_AND) {
		std::cout << "Expanded (FST) "<< node->toString() << std::endl;
		std::cout.flush();
//	}
#endif

	// update the q-values and corresponding q-markings
//	update(node);

	// update the f-values (ie upper-bounds)
	updateF(node);

	// dump the node
	//std::cout << "Expanded " << SearchNode::toString(node) << std::endl;
}

// Expand and revise (from the BST)
void AnyAAOBF::expandAndRevise(AobfSearchNode* node) {
	// expand node n
	assert( node->isFringe() );
	assert( node->getChildren().empty() );

	node->setFringe(false);
	expand(node);
	node->setExpanded(true);
	m_space->incNodesExpanded(node->getType());
	if (node->getType() == NODE_OR) {
		markFeasibleChild(node);
	}

#ifdef DEBUG
	if (node->getType() == NODE_AND) {
		std::cout << "Expanded (BST) "<< node->toString() << std::endl;
	}
#endif

	// update the q-values and corresponding q-markings
	update(node);

	// update the f-values (ie upper-bounds)
	updateF(node);

	// dump the node
	//std::cout << "Expanded " << SearchNode::toString(node) << std::endl;
	m_lowerBound = m_space->getRoot()->getValue();
}

// propagate the q-values
void AnyAAOBF::update(AobfSearchNode* node) {

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

	// compute best partial solution tree: Step 4 from Nilsson's
	m_tips.clear();
	findBestPartialTree();

	assert( m_space->getRoot()->isSolved() || (m_tips.size() > 0) );
}

// propagate the f-values
void AnyAAOBF::updateF(AobfSearchNode* node) {

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

	// compute best partial solution tree: Step 4 from Nilsson's
	//m_tipsF.clear();
	//findFeasiblePartialTree();

	// dump tips
//	std::cout << "tipsF: (" << m_tipsF.size() << "):" << std::endl;
//	for (size_t ii = 0; ii < m_tipsF.size(); ++ii) {
//		std::cout << " " << m_tipsF[ii]->toString() << std::endl;
//	}
//	std::cout.flush();
}

// Find current feasible partial solution tree
// compute the best partial solution tree by tracing down the max e-values at
// each of the OR nodes (we should also ignore already solved nodes)
bool AnyAAOBF::findFeasiblePartialTree() {

	// clear existing tip nodes
	m_tipsF.clear();

	// get the root of the search space
	AobfSearchNode *root = m_space->getRoot();
	if ( !root->isSolved() ) {

		double g = 0.0, h = 0.0;

		// clear the previous value assignment
		std::fill(m_assignmentF.begin(), m_assignmentF.end(), UNKNOWN);

		std::stack<AobfSearchNode*> s;
		s.push(root);

		while (!s.empty()) {

			AobfSearchNode *e = s.top();
			s.pop();

			switch (e->getType()) {
			case NODE_AND:
				{
					// add all its pseudo-tree children (if any)
					int var = e->getVar();
					int val = e->getVal();
					m_assignmentF[var] = val;

					std::list<AobfSearchNode*>& succ = e->getChildren();
					if ( succ.empty() == false ) { // expanded AND node
						std::list<AobfSearchNode*>::iterator ci = succ.begin();
						for (; ci != succ.end(); ++ci) {
							AobfSearchNode *c = (*ci);
							s.push(c);
						}
					} else { // unexpanded AND node
						if (!e->isSolved()) {
							e->setFringe(true);
							m_tipsF.push_back(e); // sort tips by depth
							h += e->getHeur();
						} else {
							g += e->getUpperBound();
						}
					}

					break;
				}
			case NODE_OR:
				{
					std::list<AobfSearchNode*>& succ = e->getChildren();
					if ( succ.empty() == false ) { // expanded OR node
						AobfSearchNode* m = e->getFeasibleChild();
						//assert(m != NULL);
						//s.push(m);
						if (m == NULL) {
							m_tipsF.clear();
							return false;
						} else {
							s.push(m);
							g += e->getWeight(m->getVal());
						}
					} else {
						if (!e->isSolved()) { // unexpanded OR node
							e->setFringe(true);
							m_tipsF.push_back(e); // sort tips by depth
							h += e->getHeur();
						} else {
							g += e->getUpperBound();
						}
					}

					break;
				}
			} // end switch.
		} // end while.

		if (g+h >= m_upperBound) {
			m_tipsF.clear();
			return false;
		}

		return true;
	} else {
		return false;
	}
}

// sort the tip nodes (ascending, descending, random shuffle)
void AnyAAOBF::arrangeTipNodesF() {

	std::sort(m_tipsF.begin(), m_tipsF.end(), CompNodeHeurAsc());
}

// select a tip node to expand
AobfSearchNode* AnyAAOBF::chooseTipNodeF() {

	if (m_tipsF.empty() == false) {
		return *m_tipsF.begin();
	} else {
		return NULL;
	}
}

// Revise the value of node 'n' based on its children values; log arithmetic
bool AnyAAOBF::revise(AobfSearchNode *node) {

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
bool AnyAAOBF::reviseF(AobfSearchNode *node) {
	bool result = true;
	switch (node->getType()) {
	case NODE_AND: // AND node
		{

			if (node->isTerminal()) {
				// terminal AND node
				node->setUpperBound(ELEM_ZERO);
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

			break;
		}
	case NODE_OR: // OR node
		{
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
			//result = true;

			break;
		}
	}

	return result;

}

// Recompute q-values of the nodes (lower bounds) starting from the leaf
// nodes in the explicated search graph. It also corrects the best-first
// markings on the arcs.
bool AnyAAOBF::repair() {

//	std::cout << "BEGIN REPAIR:\n";
	bool ok = true;
	LAYERS& S = m_space->getLayers();
	LAYERS::reverse_iterator ri = S.rbegin();
	for (; ri != S.rend(); ++ri) {
		std::list<AobfSearchNode*>& L = (*ri);
		std::list<AobfSearchNode*>::iterator li = L.begin();
		for (; li != L.end(); ++li) {
			AobfSearchNode* node = (*li);
			int var = node->getVar();
			if (m_problem->isMap(var)) {
				node->setSolved(false);
				if (node->isExpanded() == false) { // unexpanded (tip) MAP node
					double heur = node->getHeur();
					node->setValue(m_epsilon * heur);
					if (node->isTerminal() == false) node->setSolved(false);
					node->setFringe(true);
				} else { // expanded node
					revise(node);
					reviseF(node);
				}

				// correct the feasible marking
				if (node->getType() == NODE_OR && !node->isSolved()) {
					markFeasibleChild(node);
				}

				node->setRepaired(true);
//				std::cout << "  " << node->toString() << "\n";
			}
		}
	}
//	std::cout << " - lower bound: " << m_space->getRoot()->getValue() << "\n";
//	std::cout << "END REPAIR.\n";

	return ok;
}


// redo feasible markings to get the most promissing partial solution tree
bool AnyAAOBF::repairF() {
	// all nodes in the search graph
	bool ok = true;
	LAYERS& S = m_space->getLayers();
	LAYERS::reverse_iterator ri = S.rbegin();
	for (; ri != S.rend(); ++ri) {
		std::list<AobfSearchNode*>& L = (*ri);
		std::list<AobfSearchNode*>::iterator li = L.begin();
		for (; li != L.end(); ++li) {
			AobfSearchNode* node = (*li);
			if (node->getType() == NODE_OR && !node->isSolved()) {
				ok &= markFeasibleChild(node);
			}
		}
	}

	return ok;
}


// select the next best feasible AND child of an OR node.
// return 'false' if cannot find a new feasible child, otherwise return true.
bool AnyAAOBF::markFeasibleChild(AobfSearchNode* node) {

	assert(node->getType() == NODE_OR);
	assert(m_problem->isMap(node->getVar()));

	// check if OR node has an upper bound already
	double upbo = node->getUpperBound();
	if (upbo != ELEM_INF) {
		AobfSearchNode* m = NULL;
		double max_eval = -ELEM_INF;
		std::list<AobfSearchNode*>& succ = node->getChildren();
		std::list<AobfSearchNode*>::iterator ci = succ.begin();
		for (; ci != succ.end(); ++ci) {

			AobfSearchNode* ch = (*ci);
			int val = ch->getVal();
			double w = node->getWeight(val);
			double qval = ch->getValue();

			// check if we can prune this node
			double lowbo = qval + w;
			if ( lowbo >= upbo ) {
				continue; // skip this child (won't guarantee a better solution)
			}

			// check if non-terminal AND child
			double eval = ELEM_INF;
			if (ch->isTerminal() == false) {
				eval = (qval != 0 ? (upbo - w)/qval : ELEM_INF);
			} else {
				eval = (upbo - w);
			}

			// keep track of the max e-value (ties broken lexicographically)
			if (eval + DOUBLE_PRECISION > max_eval) {
				max_eval = eval;
				m = ch;
			}
		}

		//assert(m != NULL);
		node->setFeasibleChild(m); // m can be NULL in which case the parent should be considered deadend
		return true;

	} else { // no upper bound

		AobfSearchNode* m = NULL;
		std::list<AobfSearchNode*>& succ = node->getChildren();
		if (succ.empty() == false) {
			double min_qval = ELEM_INF;
			std::list<AobfSearchNode*>::iterator ci = succ.begin();
			for (; ci != succ.end(); ++ci) {

				AobfSearchNode* ch = (*ci);
				int val = ch->getVal();
				double w = node->getWeight( val );
				double qval = ch->getValue();

				if ( (qval + w + DOUBLE_PRECISION) < min_qval ) {
					min_qval = (w + qval);
					m = ch;
				}
			}

			//assert(m != NULL);
			node->setFeasibleChild(m);
		}

		return true;
	}
}

/*
// select the next best feasible AND child of an OR node
bool AnyAOBF::markFeasibleChild(AobfSearchNode* node) {

	assert(node->getType() == NODE_OR);
	assert(m_problem->isMap(node->getVar()));

	if (m_selection == MARK_FEASIBLE_TIGHTEST) {
		if (node->getUpperBound() != INFINITY) {
			double U = node->getUpperBound();
			AobfSearchNode* m = NULL;
			double max_e_val = -ELEM_INF;
			std::list<AobfSearchNode*>& succ = node->getChildren();
			std::list<AobfSearchNode*>::iterator ci = succ.begin();
			for (; ci != succ.end(); ++ci) {

				AobfSearchNode* ch = (*ci);
				int val = ch->getVal();
				double w = node->getWeight(val);

				// check if solved AND child
				if (ch->isSolved()) {
					m = ch;
					break;
				}

				// check if non-terminal AND child
				double q = ch->getValue();
				double e_val = w + q;
				if ( e_val < U ) {
					if ( e_val > max_e_val ) {
						max_e_val = e_val;
						m = ch;
					}
				}
			}

			if (m != NULL) {
				node->setFeasibleChild(m);
			}

			assert( node->getFeasibleChild() != NULL );

			return true;
		} else { // no upper bound

			AobfSearchNode* m = NULL;
			std::list<AobfSearchNode*>& succ = node->getChildren();
			if (succ.empty() == false) {
				double min_h_val = INFINITY;
				std::list<AobfSearchNode*>::iterator ci = succ.begin();
				for (; ci != succ.end(); ++ci) {

					AobfSearchNode* ch = (*ci);
					int val = ch->getVal();
					double w = node->getWeight( val );
					double q = ch->getValue();

					if ( (w + q + DOUBLE_PRECISION) < min_h_val ) {
						min_h_val = (w + q);
						m = ch;
					}
				}

				assert(m != NULL);
				node->setFeasibleChild(m);
			}

			return true;
		}
	} else if (m_selection == MARK_FEASIBLE_EPSILON) { // default

		if (node->getUpperBound() != ELEM_INF) {
			AobfSearchNode* m = NULL, *argmin_f_val = NULL;
			double max_e_val = -ELEM_INF;
			double min_f_val = ELEM_INF;
			std::list<AobfSearchNode*>& succ = node->getChildren();
			std::list<AobfSearchNode*>::iterator ci = succ.begin();
			for (; ci != succ.end(); ++ci) {

				AobfSearchNode* ch = (*ci);
				double e_val = ELEM_INF;
				int val = ch->getVal();
				double w = node->getWeight(val);

				// check if solved AND child
				if (ch->isSolved()) {
					double fval = ch->getUpperBound() + w;
					if (fval < min_f_val) {
						min_f_val = fval;
						argmin_f_val = ch;
					}
					continue;
				}

				// check if non-terminal AND child
				if (ch->isTerminal() == false) {

					double q = ch->getValue();
					double G = node->getUpperBound();
					e_val = (q != 0 ? (G - w)/q : ELEM_INF);
	//				e_val = (ch->q() != 0 ? (n->f() - w) / ch->q() : PLUS_INFINITY);
	//				if (ch->feasible() == true) {
	//					e_val = (n->f() - w) / ch->q();
	//				} else {
	//					e_val = (n->f() - w) / ch->q();
	//				}
				} else {
					e_val = (node->getUpperBound() - w);
				}

				// keep track of the max e-value (ties broken lexicographically)
				if (e_val + DOUBLE_PRECISION > max_e_val) {
					max_e_val = e_val;
					m = ch;
				}
			}

			//assert(m != NULL);
			if (argmin_f_val != NULL) node->setFeasibleChild(argmin_f_val);
			else node->setFeasibleChild(m);

			// check if lower-upper bound violation
//			double LB = node->getFeasibleChild()->getValue();
//			int ii = node->getFeasibleChild()->getVal();
//			double ww = node->getWeight(ii);
//			LB += ww;
//			double UB = node->getUpperBound();
//			if (LB >= UB) {
//				std::cout << "*** LB >= UB " << LB << " " << UB << "\n";
//				std::list<AobfSearchNode*>::const_iterator it = node->getChildren().begin();
//				for (; it != node->getChildren().end(); ++it) {
//					AobfSearchNode* ch = *it;
//					int jj = ch->getVal();
//					std::cout << " weight=" << node->getWeight(jj) << " ch: " << ch->toString() << "\n";
//				}
//			}

			return true;
		} else { // no upper bound

			AobfSearchNode* m = NULL;
			std::list<AobfSearchNode*>& succ = node->getChildren();
			if (succ.empty() == false) {
				double min_h_val = ELEM_INF;
				std::list<AobfSearchNode*>::iterator ci = succ.begin();
				for (; ci != succ.end(); ++ci) {

					AobfSearchNode* ch = (*ci);
					int val = ch->getVal();
					double w = node->getWeight( val );
					double q = ch->getValue();

					if ( (w + q + DOUBLE_PRECISION) < min_h_val ) {
						min_h_val = (w + q);
						m = ch;
					}
				}

				assert(m != NULL);
				node->setFeasibleChild(m);
			}

			return true;
		}
	} else {
		throw "Selection strategy not implemented.";
	}
}
*/
