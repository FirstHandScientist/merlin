/*
 * any_lrtao.cpp
 *
 *  Created on: 4 Mar 2016
 *      Author: radu
 */


#include "any_ldfs_bounded.h"


// AO* search over the context minimal AND/OR graph restricted to MAP vars
int AnyLDFS_Bounded::aostar(bool verbose) {

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
	m_space->setRoot(root);
	m_space->add(state, root);

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
					<< std::setw(12) << m_upperBound << " (" << 1/ELEM_DECODE(m_upperBound) << ") "
					<< std::setw(12) << m_lowerBound << " (" << 1/ELEM_DECODE(m_lowerBound) << ")"
					<< std::endl;
			}
#endif
		}

		// AOstar graph search
		AobfSearchNode* n = NULL;
		bool ok = findBestPartialTree();
		if (ok) {

			// check if solution tree;
			if (m_tips.empty()) {
				// propagate the values of the leaf nodes to get the solution
				for (std::vector<AobfSearchNode*>::iterator li = m_leaves.begin();
						li != m_leaves.end(); ++li) {
					AobfSearchNode* leaf = *li; // could also be OR sum nodes
					update(leaf);
				}
			} else {
				arrangeTipNodes(); 			// order the tip nodes by q-values
				n = chooseTipNode();		// select a tip node
				expandAndRevise(n);			// expand and revise the q-value of the node
											// also mark the best child node (if OR)
			}
		} else {
			// the partial solution tree cannot be extended
			// propagate the values of the leaf nodes to get the solution
			for (std::vector<AobfSearchNode*>::iterator li = m_leaves.begin();
					li != m_leaves.end(); ++li) {
				AobfSearchNode* leaf = *li; // could also be OR sum nodes
				update(leaf);
			}

			// pick a tip node, expand it and update the ancestors values
			if (m_tips.empty() == false) {
				arrangeTipNodes();
				n = chooseTipNode();
				expandAndUpdate(n);
			}
		}

		// check the new upper and lower bounds
		m_lowerBound = root->getValue();
		if (root->getUpperBound() < m_upperBound) {
			m_upperBound = root->getUpperBound();
			std::cout << "["
				<< std::setw(9) << m_timer.elapsed() << "] w "
				<< std::setw(8) << m_epsilon << " "
				<< std::setw(12) << m_space->getOrNodes() << " "
				<< std::setw(12) << m_space->getAndNodes() << " "
				<< std::setw(12) << m_upperBound << " (" << 1/ELEM_DECODE(m_upperBound) << ") "
				<< std::setw(12) << m_lowerBound << " (" << 1/ELEM_DECODE(m_lowerBound) << ")"
				<< std::endl;
		}

		// check the termination condition
		if (m_upperBound == m_lowerBound) {
			root->setSolved(true);
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
		<< std::setw(12) << m_upperBound << " (" << 1/ELEM_DECODE(m_upperBound) << ") "
		<< std::setw(12) << m_lowerBound << " (" << 1/ELEM_DECODE(m_lowerBound) << ")"
		<< std::endl;

	return result;
}



// solve the problem
int AnyLDFS_Bounded::solve() {

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
bool AnyLDFS_Bounded::expand(AobfSearchNode *node) {

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

// Expand and revise the node value
void AnyLDFS_Bounded::expandAndRevise(AobfSearchNode* node) {
	// expand node n
	assert( node->isFringe() );
	assert( node->getChildren().empty() );

	node->setFringe(false);
	expand(node);
	node->setExpanded(true);
	m_space->incNodesExpanded(node->getType());

	// revise the q-value of this node (don't propagate yet)
	revise(node);

#ifdef DEBUG
	std::cout << "Expanded "<< node->toString() << std::endl;
#endif

}

// Expand and update the node values of the ancestors
void AnyLDFS_Bounded::expandAndUpdate(AobfSearchNode* node) {
	// expand node n
	assert( node->isFringe() );
	assert( node->getChildren().empty() );

	node->setFringe(false);
	expand(node);
	node->setExpanded(true);
	m_space->incNodesExpanded(node->getType());

	// revise the q-value of this node and propagate it to the ancestors
	update(node);

#ifdef DEBUG
	std::cout << "Expanded "<< node->toString() << std::endl;
#endif

}

// propagate the q-values
void AnyLDFS_Bounded::update(AobfSearchNode* node) {

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
					// new parent through marked connector
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

// Find current best partial solution tree by following down marked arcs.
// Returns true if the f-value of the best tree is smaller than the current U,
// otherwise it returns false, namely the current partial tree is a deadend.
bool AnyLDFS_Bounded::findBestPartialTree() {

	// clear existing tip nodes
	m_tips.clear();
	m_leaves.clear();

	// get the root of the search space
	AobfSearchNode *root = m_space->getRoot();
	double g = ELEM_ZERO, h = ELEM_ZERO;

	// clear the previous value assignment
	std::fill(m_assignment.begin(), m_assignment.end(), UNKNOWN);

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
				m_assignment[var] = val;

				std::list<AobfSearchNode*>& succ = e->getChildren();
				if ( succ.empty() == false ) { // expanded AND node
					std::list<AobfSearchNode*>::iterator ci = succ.begin();
					for (; ci != succ.end(); ++ci) {
						AobfSearchNode *c = (*ci);
						s.push(c);
						c->setCurrentParent(e);
					}
				} else { // unexpanded AND node
					if (!e->isSolved()) {
						e->setFringe(true);
						m_tips.push_back(e); // sort tips by depth
						h += e->getHeur();
					} else { // terminal AND node
						m_leaves.push_back(e);
						g += e->getValue();
					}
				}

				break;
			}
		case NODE_OR:
			{
				std::list<AobfSearchNode*>& succ = e->getChildren();
				if ( succ.empty() == false ) { // expanded OR node
					AobfSearchNode* m = e->getBestChild();
					assert(m != NULL);
					s.push(m);
					g += e->getWeight(m->getVal());
					m->setCurrentParent(e);
				} else { // unexpanded OR node
					if (!e->isSolved()) { // unexpanded OR node
						e->setFringe(true);
						m_tips.push_back(e); // sort tips by depth
						h += e->getValue();
					} else { // terminal OR node labeled by SUM variable
						m_leaves.push_back(e);
						g += e->getUpperBound();
					}
				}

				break;
			}
		} // end switch.
	} // end while.

	// Check if current partial tree is can be extended
	bool ok = true;
	if (g + h >= m_upperBound) {
		ok = false;
	}

	return ok;
}

// Check if a tip node is deadend (node must not be expanded yet); it is supposed
// to be called after a best partial solution tree has been selected/found.
bool AnyLDFS_Bounded::isDeadend(AobfSearchNode* tip) {

	assert(tip->isExpanded() == false);
	assert(m_problem->isMap(tip->getVar()));

	bool deadend = false;
	if (tip->getType() == NODE_OR) { // tip node is an OR node

		// get the current lower bound (q-value)
		double lb = tip->getValue();
		AobfSearchNode *prevOR = tip; // OR level (current)
		AobfSearchNode *prevAND = tip->getCurrentParent(); // AND level
		while (prevAND != NULL) { // loop until we reach the root level
			std::list<AobfSearchNode*>& children = prevAND->getChildren();
			std::list<AobfSearchNode*>::const_iterator ci = children.begin();
			for (; ci != children.end(); ++ci) {
				if (*ci != prevOR) {
					lb += (*ci)->getValue(); // the q-value
				}
			}

			prevOR = prevAND->getCurrentParent(); // OR level
			double w = prevOR->getWeight(prevAND->getVal());
			lb += w; // add the current arc weight
			double ub = prevOR->getUpperBound(); // upper bound at the current OR level
			if (lb > ub + DOUBLE_PRECISION) {
				deadend = true;
				break;
			}

			// move up another level
			prevAND = prevOR->getCurrentParent();
		}

	} else { // tip node is an AND node

		// get current lower bound (q-value)
		double lb = tip->getValue();

		// otherwise, go upward to the root and check the upper bound at each
		// OR ancestor along the current best partial solution tree
		AobfSearchNode *prevOR = tip->getCurrentParent(); // OR level (current)
		AobfSearchNode *prevAND = tip; // AND level
		while (prevAND != NULL) { // loop until we reach the root level
			std::list<AobfSearchNode*>& children = prevAND->getChildren();
			std::list<AobfSearchNode*>::const_iterator ci = children.begin();
			for (; ci != children.end(); ++ci) {
				if (*ci != prevOR) {
					lb += (*ci)->getValue(); // the q-value
				}
			}

			double w = prevOR->getWeight(prevAND->getVal());
			lb += w; // add the current arc weight
			double ub = prevOR->getUpperBound(); // upper bound at the current OR level
			if (lb > ub + DOUBLE_PRECISION) {
				deadend = true;
				break;
			}

			// move up another level
			prevAND = prevOR->getCurrentParent(); // AND level
			if (prevAND != NULL)
				prevOR = prevAND->getCurrentParent(); // OR level
									// otherwise AND level is NULL (root)
		}

	}

	return deadend; // true if the tip node is a deadend
}

// Revise the q-value and the u-value of a node based on its children values; log arithmetic
bool AnyLDFS_Bounded::revise(AobfSearchNode *node) {

	// safety checks
	assert(node);
	assert(node->getType() == NODE_AND || node->getType() == NODE_OR);
	//assert(m_problem->isMap(node->getVar()));

	bool change = true;
	if (node->getType() == NODE_AND) { // revise an AND node
		if (node->isTerminal()) { // parent OR is a MAP variable that is leaf in pseudo tree

			// terminal AND node
			node->setValue(ELEM_ZERO);
			node->setUpperBound(ELEM_ZERO);
			node->setSolved(true);
			node->setFringe(false);

		} else {

			// non-terminal AND node
			bool solved = true;
			double qval = ELEM_ZERO, uval = ELEM_ZERO;
			std::list<AobfSearchNode*>::const_iterator li =
					node->getChildren().begin();
			for (; li != node->getChildren().end(); ++li) {
				AobfSearchNode *ch = (*li);
				bool sol = ch->isSolved();
				solved &= sol;
				qval += ch->getValue();
				uval += ch->getUpperBound();
			}

			// label the node 'solved' if all its OR children are 'solved'
			node->setValue(qval);
			node->setUpperBound(uval);
			if (solved) {
				node->setSolved(true);
				node->setFringe(false);
			}
		}
	} else { // revise the cost of an OR node

		if (node->isTerminal()) {

			// terminal OR node ie. SUM variable (q-value and u-value set during its parent expansion)
			node->setSolved(true);
			node->setFringe(false);

		} else {

			double qval = INFINITY, uval = INFINITY;
			bool solved = true;
			std::list<AobfSearchNode*>::const_iterator li = node->getChildren().begin();
			for (; li != node->getChildren().end(); ++li) {
				AobfSearchNode *ch = (*li);
				int val = ch->getVal();
				double w = node->getWeight(val);
				double q = w + ch->getValue(); // could be infinity
				double u = w + ch->getUpperBound();
				solved &= ch->isSolved();
				qval = std::min(qval, q);
				uval = std::min(uval, u);
			}

			// mark the best child
			setBestChild(node);

			// label the node 'solved' if all its children are 'solved'
			node->setValue(qval);
			node->setUpperBound(uval);
			if (solved) {
				node->setSolved(true);
				node->setFringe(false);
			}
		}
	}

	return change; // should be true to propagate all the way up to the root
}

// Find the best child of an OR node (ie, minimizing w+q). Use --positive to
// avoid numerical issues when dealing with zeros (log 0 = inf)
bool AnyLDFS_Bounded::setBestChild(AobfSearchNode* node) {

	// safety checks: OR node labeled by MAP variable
	assert(node->getType() == NODE_OR);
	assert(m_problem->isMap(node->getVar()));

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

		assert(m != NULL);
		node->setBestChild(m); // keep track of the best child (ie, min q-value)
	}

	return true;
}
