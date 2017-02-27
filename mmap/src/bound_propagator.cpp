/*
 * bound_propagator.cpp
 *
 *  Created on: Jun 7, 2013
 *      Author: radu
 */

#include "bound_propagator.h"
#include "timer.h"

//#undef DEBUG

/**
 * For MAP computation, AND nodes do multiplication, while OR nodes do either
 * summation (if VAR_SUM) or maximization (if VAR_MAX -- MAP variable)
 */

AobbSearchNode* BoundPropagator::propagate(AobbSearchNode* n, Timer& tm,
		std::vector<std::vector<bool> >& currentDomains,
		bool reportSolution, AobbSearchNode* upperLimit) {

	// 'n' is a leaf node (which triggers the propagation)
	assert(n->getChildCountAct() == 0);

	// these two pointers move upward in the search space, always one level
	// apart s.t. cur is the parent node of prev
	AobbSearchNode* cur = n->getParent(), *prev = n;
	if (cur == NULL) { // n is actually the root
		m_lowerBound = n->getValue();
		return NULL;
	}

	// Keeps track of the highest node to be deleted during cleanup,
	// where .second will be deleted as a child of .first
	pair<AobbSearchNode*, AobbSearchNode*> highestDelete(NULL, NULL);
	int numRollbacks = 0, numBacktracks = 0;
	list<pair<AobbSearchNode*, AobbSearchNode*> > rollbacks, backtracks;

	// 'prop' signals whether we are still propagating values in this call
	bool prop = true;
	// 'del' signals whether we are still deleting nodes in this call
	bool del = (upperLimit != n) ? true : false;

//#ifdef DEBUG
//	cout << " --- begin propagation:" << endl;
//#endif

	// going all the way to the root, if we have to
	do {

		if (cur->getType() == NODE_AND) {
			// =================================================================

			// propagate and update node values
			if (prop) {
				// optimal solution to previously solved and deleted child OR nodes
				double d = cur->getSubSolved();
				// current best solution to yet-unsolved OR child nodes
				NodeP* children = cur->getChildren();
				for (size_t i = 0; i < cur->getChildCountFull(); ++i) {
					if (children[i]) {
						d *= children[i]->getValue();
					}
				}

				// store into value (thus includes cost of subSolved)
				cur->setValue(d);

//#ifdef DEBUG
//			cout << "   current AND node updated: " << cur->toString() << endl;
//#endif

				if ( ISNAN(d) ) { // not all OR children solved yet, propagation stops here
					prop = false;
#ifndef NO_ASSIGNMENT
					propagateTuple(n, cur); // save (partial) opt. subproblem solution at current AND node
#endif
				}

			}

			// clean up fully propagated nodes (i.e. no children or only one (=prev))
			if (del) {
				if (prev->getChildCountAct() <= 1) {

					// prev is OR node, try to cache
					if (m_doCaching && prev->isCachable() /*&& !m_problem->isMap(prev->getVar())*/
							&& !prev->isNotOpt()) {
						try {
							//assert( ISNAN(prev->getValue()) == false );
//							m_space->cache->write(prev->getVar(),
//								prev->getCacheContext(), prev->getValue());
#ifndef NO_ASSIGNMENT
							m_space->cache->write(prev->getVar(),
								prev->getCacheContext(), prev->getValue(),
								prev->getOptAssig() );
#else
							m_space->cache->write(prev->getVar(),
								prev->getCacheContext(), prev->getValue() );
#endif
						} catch (...) { /* tried to cache NaN value */
						}
					}

					highestDelete = make_pair(cur, prev);
				} else {
					del = false; // there are unsolved children
				}

			}

			// ===========================================================================
		} else { // cur is OR node
			// ===========================================================================

			if (prop) {
				double d = prev->getValue() * prev->getLabel(); // getValue includes subSolved

				if (m_problem->isMap(cur->getVar())) { // MAP variable
					if (ISNAN(cur->getValue()) || d > cur->getValue()) {
						cur->setValue(d); // update max value
					} else {
						prop = false; // no more value propagation upwards in this call
					}
				} else { // SUM var
					if (m_pseudotree->isSumRoot(cur->getVar())) {
						m_space->numSumEvals++;
					}
					if (ISNAN(cur->getValue())) {
						cur->setValue(d); // update sum value
					} else {
						cur->setValue( d + cur->getValue() );
					}

					if (cur->getChildCountAct() > 1) {
						prop = false; // SUM nodes must wait for all children to
									  // be evaluated before propagating their
									  // values upwards in the search graph.
					}
				}
			}

//#ifdef DEBUG
//			cout << "   current OR node updated: " << cur->toString() << endl;
//#endif

			if (del) {
				if (prev->getChildCountAct() <= 1) { // prev has no or one children?
					highestDelete = make_pair(cur, prev);

					if (cur->getType() == NODE_OR &&
						prev->getType() == NODE_AND) {
						numBacktracks++;
						backtracks.push_back(make_pair(cur, prev));
					}

					if (cur->getType() == NODE_OR &&
						prev->getType() == NODE_AND &&
						m_problem->isMap(cur->getVar()) &&
						m_heuristic != NULL) {

						if (prev->isExpanded()) {
							numRollbacks++;
							rollbacks.push_back(make_pair(cur, prev));
						}
					}
				} else {
					del = false;
				}
			}

#ifndef NO_ASSIGNMENT
		  // save opt. tuple, will be needed for caching later
		  if ( prop && cur->isCachable() && !cur->isNotOpt() ) {
			  //DIAG(myprint("< Cachable OR node found\n"));
			  propagateTuple(n, cur);
		  }
#endif

		  // ===========================================================================
		}

		// don't delete anything higher than upperLimit
		if (upperLimit == cur)
			del = false;

		// move pointers up in search space
		if (prop || del) {
			prev = cur;
			cur = cur->getParent();
		} else {
			break;
		}

	} while (cur); // until cur==NULL, i.e. 'parent' of root

//#ifdef DEBUG
//	cout << " --- end propagation." << endl;
//#endif

	// propagated up to root node, update tuple as well
	if (prop && !cur) {
#ifndef NO_ASSIGNMENT
		propagateTuple(n, prev);
#endif
		if (reportSolution) {
			double timestamp = tm.elapsed();
#ifndef NO_ASSIGNMENT
			m_solutions.push_back(prev->getOptAssig());
			updateSolution(timestamp,
					prev->getValue(),
					prev->getOptAssig(),
					std::make_pair(m_space->nodesORmap, m_space->nodesANDmap),
					std::make_pair(m_space->nodesOR, m_space->nodesAND));
#else
			updateSolution(timestamp,
					prev->getValue(),
					std::make_pair(m_space->nodesORmap, m_space->nodesANDmap),
					std::make_pair(m_space->nodesOR, m_space->nodesAND));
#endif
		}
	}

	if (highestDelete.first) {
		AobbSearchNode* parent = highestDelete.first;
		AobbSearchNode* child = highestDelete.second;
		if (parent->getType() == NODE_AND) {
			// Store value of OR node to be deleted into AND parent
			parent->addSubSolved(child->getValue());
		}

		// rollback the heuristic
		if (numRollbacks > 0) {
			assert(numRollbacks == (int) rollbacks.size());
			list<pair<AobbSearchNode*, AobbSearchNode*> >::iterator li = rollbacks.begin();
			for (; li != rollbacks.end(); ++li) {
				AobbSearchNode* c = (*li).second;
//				std::cout << "Rollback " << c->toString() << std::endl;
				m_heuristic->check(c->getDepth());
				m_heuristic->rollback();
			}
		}

		// backtrack the SAT solver (manually)
		if (numBacktracks > 0) {
			assert(numBacktracks == (int) backtracks.size());
			list<pair<AobbSearchNode*, AobbSearchNode*> >::iterator li = backtracks.begin();
			for (; li != backtracks.end(); ++li) {
//				std::cout << "Backtrack " << (*li).second->toString() << std::endl;
				if (m_zchaff) {
					if (!currentDomains.empty()) {
						std::list<std::pair<int,int> >::iterator ci = li->second->changes().begin();
						for (; ci != li->second->changes().end(); ++ci) {
							std::pair<int,int> &p = (*ci);
							currentDomains[p.first][p.second] = true;
						}
					}

					m_zchaff->release_assignments();
					m_zchaff->dlevel()--;

				}
			}

		}

		// finally clean up, delete subproblem with unnecessary nodes from memory
		parent->eraseChild(child);
	}

	// keep the most recent best solution found so far
	m_lowerBound = m_space->root->getValue();

	return highestDelete.first;
}

/**
 * For MAP computation, AND nodes do multiplication, while OR nodes do either
 * summation (if VAR_SUM) or maximization (if VAR_MAX -- MAP variable)
 */

AobbSearchNode* BoundPropagator::propagate(AobbSearchNode* n, Timer& tm,
		bool reportSolution, AobbSearchNode* upperLimit) {

	// 'n' is a leaf node (which triggers the propagation)
	assert(n->getChildCountAct() == 0);

	// these two pointers move upward in the search space, always one level
	// apart s.t. cur is the parent node of prev
	AobbSearchNode* cur = n->getParent(), *prev = n;

	// Keeps track of the highest node to be deleted during cleanup,
	// where .second will be deleted as a child of .first
	pair<AobbSearchNode*, AobbSearchNode*> highestDelete(NULL, NULL);

	// 'prop' signals whether we are still propagating values in this call
	bool prop = true;
	// 'del' signals whether we are still deleting nodes in this call
	bool del = (upperLimit != n) ? true : false;

//#ifdef DEBUG
//	cout << " --- begin propagation:" << endl;
//#endif

	// going all the way to the root, if we have to
	do {

		if (cur->getType() == NODE_AND) {
			// =================================================================

			// propagate and update node values
			if (prop) {
				// optimal solution to previously solved and deleted child OR nodes
				double d = cur->getSubSolved();
				// current best solution to yet-unsolved OR child nodes
				NodeP* children = cur->getChildren();
				for (size_t i = 0; i < cur->getChildCountFull(); ++i) {
					if (children[i]) {
						d *= children[i]->getValue();
					}
				}

				// store into value (thus includes cost of subSolved)
				cur->setValue(d);

//#ifdef DEBUG
//			cout << "   current AND node updated: " << cur->toString() << endl;
//#endif

				if ( ISNAN(d) ) { // not all OR children solved yet, propagation stops here
					prop = false;
				}

			}

			// clean up fully propagated nodes (i.e. no children or only one (=prev))
			if (del) {
				if (prev->getChildCountAct() <= 1) {

					// prev is OR node, try to cache
					if (m_doCaching && prev->isCachable() /*&& !m_problem->isMap(prev->getVar())*/
							&& !prev->isNotOpt()) {
						try {
							//assert( ISNAN(prev->getValue()) == false );
//							m_space->cache->write(prev->getVar(),
//								prev->getCacheContext(), prev->getValue() );

#ifndef NO_ASSIGNMENT
							m_space->cache->write(prev->getVar(),
								prev->getCacheContext(), prev->getValue(),
								prev->getOptAssig() );
#else
							m_space->cache->write(prev->getVar(),
								prev->getCacheContext(), prev->getValue() );
#endif

						} catch (...) { /* tried to cache NaN value */
						}
					}

					highestDelete = make_pair(cur, prev);
				} else {
					del = false; // there are unsolved children
				}

			}

			// ===========================================================================
		} else { // cur is OR node
			// ===========================================================================

			if (prop) {
				double d = prev->getValue() * prev->getLabel(); // getValue includes subSolved

				if (m_problem->isMap(cur->getVar())) { // MAP variable
					if (ISNAN(cur->getValue()) || d > cur->getValue()) {
						cur->setValue(d); // update max value
					} else {
						prop = false; // no more value propagation upwards in this call
					}
				} else { // SUM var
					if (ISNAN(cur->getValue())) {
						cur->setValue(d); // update sum value
					} else {
						cur->setValue( d + cur->getValue() );
					}

					if (cur->getChildCountAct() > 1) {
						prop = false; // SUM nodes must wait for all children to
									  // be evaluated before propagating their
									  // values upwards in the search graph.
					}
				}
			}

//#ifdef DEBUG
//			cout << "   current OR node updated: " << cur->toString() << endl;
//#endif

			if (del) {
				if (prev->getChildCountAct() <= 1) { // prev has no or one children?
					highestDelete = make_pair(cur, prev);
				} else {
					del = false;
				}
			}

		  // ===========================================================================
		}

		// don't delete anything higher than upperLimit
		if (upperLimit == cur)
			del = false;

		// move pointers up in search space
		if (prop || del) {
			prev = cur;
			cur = cur->getParent();
		} else {
			break;
		}

	} while (cur); // until cur==NULL, i.e. 'parent' of root

//#ifdef DEBUG
//	cout << " --- end propagation." << endl;
//#endif

	// propagated up to root node, update tuple as well
//	if (prop && !cur) {
//		if (reportSolution) {
//			double timestamp = tm.elapsed();
//			updateSolution(timestamp,
//					prev->getValue(),
//					std::make_pair(m_space->nodesORmap, m_space->nodesANDmap),
//					std::make_pair(m_space->nodesOR, m_space->nodesAND));
//		}
//	}
	// propagated up to root node, update tuple as well
	if (prop && !cur) {
#ifndef NO_ASSIGNMENT
		propagateTuple(n, prev);
#endif
		if (reportSolution) {
			double timestamp = tm.elapsed();
#ifndef NO_ASSIGNMENT
			m_solutions.push_back(prev->getOptAssig());
			updateSolution(timestamp,
					prev->getValue(),
					prev->getOptAssig(),
					std::make_pair(m_space->nodesORmap, m_space->nodesANDmap),
					std::make_pair(m_space->nodesOR, m_space->nodesAND));
#else
			updateSolution(timestamp,
					prev->getValue(),
					std::make_pair(m_space->nodesORmap, m_space->nodesANDmap),
					std::make_pair(m_space->nodesOR, m_space->nodesAND));
#endif
		}
	}

	if (highestDelete.first) {
		AobbSearchNode* parent = highestDelete.first;
		AobbSearchNode* child = highestDelete.second;
		if (parent->getType() == NODE_AND) {
			// Store value of OR node to be deleted into AND parent
			parent->addSubSolved(child->getValue());
		}

		// finally clean up, delete subproblem with unnecessary nodes from memory
		parent->eraseChild(child);
	}

	return highestDelete.first;
}

//// Node values propagation for the MIN-SUM formulation (ie, log space)
//AobbSearchNode* BoundPropagator::propagateLog(AobbSearchNode* n, Timer& tm,
//		bool reportSolution, AobbSearchNode* upperLimit) {
//
//	// 'n' is a leaf node (which triggers the propagation)
//	assert(n->getChildCountAct() == 0);
//
//	// these two pointers move upward in the search space, always one level
//	// apart s.t. cur is the parent node of prev
//	AobbSearchNode* cur = n->getParent(), *prev = n;
//
//	// Keeps track of the highest node to be deleted during cleanup,
//	// where .second will be deleted as a child of .first
//	pair<AobbSearchNode*, AobbSearchNode*> highestDelete(NULL, NULL);
//
//	// 'prop' signals whether we are still propagating values in this call
//	bool prop = true;
//	// 'del' signals whether we are still deleting nodes in this call
//	bool del = (upperLimit != n) ? true : false;
//
////#ifdef DEBUG
////	cout << " --- begin propagation:" << endl;
////#endif
//
//	// going all the way to the root, if we have to
//	do {
//
//		if (cur->getType() == NODE_AND) {
//			// =================================================================
//
//			// propagate and update node values
//			if (prop) {
//				// optimal solution to previously solved and deleted child OR nodes
//				double d = cur->getSubSolved();
//				// current best solution to yet-unsolved OR child nodes
//				NodeP* children = cur->getChildren();
//				for (size_t i = 0; i < cur->getChildCountFull(); ++i) {
//					if (children[i]) {
//						d += children[i]->getValue();
//					}
//				}
//
//				// store into value (thus includes cost of subSolved)
//				cur->setValue(d);
//
////#ifdef DEBUG
////			cout << "   current AND node updated: " << cur->toString() << endl;
////#endif
//
//				if ( ISNAN(d) ) { // not all OR children solved yet, propagation stops here
//					prop = false;
//				}
//
//			}
//
//			// clean up fully propagated nodes (i.e. no children or only one (=prev))
//			if (del) {
//				if (prev->getChildCountAct() <= 1) {
//
//					// prev is OR node, try to cache
//					if (m_doCaching && prev->isCachable() /*&& !m_problem->isMap(prev->getVar())*/
//							&& !prev->isNotOpt()) {
//						try {
//							//assert( ISNAN(prev->getValue()) == false );
//							m_space->cache->write(prev->getVar(),
//								prev->getCacheContext(), prev->getValue() );
//						} catch (...) { /* tried to cache NaN value */
//						}
//					}
//
//					highestDelete = make_pair(cur, prev);
//				} else {
//					del = false; // there are unsolved children
//				}
//
//			}
//
//			// ===========================================================================
//		} else { // cur is OR node
//			// ===========================================================================
//
//			if (prop) {
//				double d = prev->getValue() + prev->getLabel(); // getValue includes subSolved
//
//				assert(m_problem->isMap(cur->getVar()));
//				if (ISNAN(cur->getValue()) || d < cur->getValue()) {
//					cur->setValue(d); // update max value
//				} else {
//					prop = false; // no more value propagation upwards in this call
//				}
//			}
//
////#ifdef DEBUG
////			cout << "   current OR node updated: " << cur->toString() << endl;
////#endif
//
//			if (del) {
//				if (prev->getChildCountAct() <= 1) { // prev has no or one children?
//					highestDelete = make_pair(cur, prev);
//				} else {
//					del = false;
//				}
//			}
//
//		  // ===========================================================================
//		}
//
//		// don't delete anything higher than upperLimit
//		if (upperLimit == cur)
//			del = false;
//
//		// move pointers up in search space
//		if (prop || del) {
//			prev = cur;
//			cur = cur->getParent();
//		} else {
//			break;
//		}
//
//	} while (cur); // until cur==NULL, i.e. 'parent' of root
//
////#ifdef DEBUG
////	cout << " --- end propagation." << endl;
////#endif
//
//	// propagated up to root node, update tuple as well
//	if (prop && !cur) {
//		if (reportSolution) {
//			double timestamp = tm.elapsed();
//			updateSolutionLog(timestamp,
//					prev->getValue(),
//					prev->getHeurUpdated(),
//					std::make_pair(m_space->nodesORmap, m_space->nodesANDmap));
//		}
//	}
//
//	if (highestDelete.first) {
//		AobbSearchNode* parent = highestDelete.first;
//		AobbSearchNode* child = highestDelete.second;
//		if (parent->getType() == NODE_AND) {
//			// Store value of OR node to be deleted into AND parent
//			parent->setSubSolved( parent->getSubSolved() + child->getValue() );
//			//parent->addSubSolved(child->getValue());
//		}
//
//		// finally clean up, delete subproblem with unnecessary nodes from memory
//		parent->eraseChild(child);
//	}
//
//	return highestDelete.first;
//}

//// Node values propagation for the MIN-SUM formulation - up to the root!
//AobbSearchNode* BoundPropagator::updateLog(AobbSearchNode* n) {
//
//	// revise n first; n was just expanded
//	if (n->getChildCountFull() > 0) {
//		if (n->getType() == NODE_OR) { // OR node
//			double q = ELEM_INF;
//			NodeP* children = n->getChildren();
//			for (size_t i = 0; i < n->getChildCountFull(); ++i) {
//				if (children[i]) {
//					q = std::min(q, children[i]->getHeurUpdated());
//				}
//			}
//
//			n->setHeurUpdated(q);
//		} else { // AND node
//			double q = ELEM_ZERO;
//			NodeP* children = n->getChildren();
//			for (size_t i = 0; i < n->getChildCountFull(); ++i) {
//				if (children[i]) {
//					q += children[i]->getHeurUpdated();
//				}
//			}
//
//			n->setHeurUpdated(q);
//		}
//	}
//
//	// these two pointers move upward in the search space, always one level
//	// apart s.t. cur is the parent node of prev
//	AobbSearchNode* cur = n->getParent(), *prev = n;
//	if (cur == NULL) return NULL;
//
//	// Keeps track of the highest node to be deleted during cleanup,
//	// where .second will be deleted as a child of .first
//	pair<AobbSearchNode*, AobbSearchNode*> highestDelete(NULL, NULL);
//
//	// 'prop' signals whether we are still propagating values in this call
//	bool prop = true;
//
////#ifdef DEBUG
////	cout << " --- begin propagation:" << endl;
////#endif
//
//	// going all the way to the root, if we have to
//	do {
//
//		if (cur->getType() == NODE_AND) {
//			// =================================================================
//
//			// propagate and update node values
//			if (prop) {
//				// optimal solution to previously solved and deleted child OR nodes
//				double d = cur->getSubSolved();
//				// current best solution to yet-unsolved OR child nodes
//				NodeP* children = cur->getChildren();
//				for (size_t i = 0; i < cur->getChildCountFull(); ++i) {
//					if (children[i]) {
//						d += children[i]->getHeurUpdated();
//					}
//				}
//
//				// store into value (thus includes cost of subSolved)
//				cur->setHeurUpdated(d);
//
////#ifdef DEBUG
////			cout << "   current AND node updated: " << cur->toString() << endl;
////#endif
//
//			}
//
//			// clean up fully propagated nodes (i.e. no children or only one (=prev))
//			if (prev->getChildCountAct() <= 1) {
//				// prev is OR and is solved
//				prev->setHeurUpdated(prev->getValue());
//				highestDelete = make_pair(cur, prev);
//			}
//
//			// ===========================================================================
//		} else { // cur is OR node
//			// ===========================================================================
//
//			if (prop) {
//				double d;
////				if (prev->getChildCountAct() <= 1) {
////					d = prev->getValue() + prev->getLabel();
////				} else {
//					d = prev->getHeurUpdated();
////				}
//
//				NodeP* children = cur->getChildren();
//				for (size_t i = 0; i < cur->getChildCountFull(); ++i) {
//					if (children[i] && children[i] != prev) {
//						d = std::min(d, children[i]->getHeurUpdated());
//					}
//				}
//
//				cur->setHeurUpdated(d);
//			}
//
////#ifdef DEBUG
////			cout << "   current OR node updated: " << cur->toString() << endl;
////#endif
//
//			if (prev->getChildCountAct() <= 1) { // prev has no or one children?
//				highestDelete = make_pair(cur, prev);
//			}
//
//		  // ===========================================================================
//		}
//
//		// move pointers up in search space
//		if (prop) {
//			prev = cur;
//			cur = cur->getParent();
//		}
//
//	} while (cur); // until cur==NULL, i.e. 'parent' of root
//
////#ifdef DEBUG
////	cout << " --- end propagation." << endl;
////#endif
//
//	// propagated up to root node, update tuple as well
//	if (prop && !cur) {
//		m_lowerBound = prev->getHeurUpdated();
//	}
//
////	if (highestDelete.first) {
////		AobbSearchNode* parent = highestDelete.first;
////		AobbSearchNode* child = highestDelete.second;
////		if (parent->getType() == NODE_AND) {
////			// Store value of OR node to be deleted into AND parent
////			parent->setSubSolved( parent->getSubSolved() + child->getValue() );
////			//parent->addSubSolved(child->getValue());
////		}
////
////		// finally clean up, delete subproblem with unnecessary nodes from memory
////		parent->eraseChild(child);
////	}
//
//	return highestDelete.first;
//}


#ifndef NO_ASSIGNMENT

/* collects the joint assignment from 'start' upwards until 'end' and
 * records it into 'end' for later use */
void BoundPropagator::propagateTuple(AobbSearchNode* start, AobbSearchNode* end) {
	assert(start && end);

	int endVar = end->getVar();
	const set<int>& endSubprob = m_pseudotree->getNode(endVar)->getSubprobVars();

	// get variable map for end node
	vector<int> endVarMap = m_pseudotree->getNode(endVar)->getSubprobVarMap();
	// allocate assignment in end node
	vector<val_t>& assig = end->getOptAssig();
	assig.resize(endSubprob.size(), UNKNOWN);

	int curVar = UNKNOWN, curVal = UNKNOWN;
	for (AobbSearchNode* cur=start; cur!=end; cur=cur->getParent()) {
		curVar = cur->getVar();

		if (cur->getType() == NODE_AND) {
			curVal = cur->getVal();
			if (curVal!=UNKNOWN)
			assig.at(endVarMap.at(curVar)) = curVal;
		}

		if (cur->getOptAssig().size()) {
			// check previously saved partial assignment
			const set<int>& curSubprob = m_pseudotree->getNode(cur->getVar())->getSubprobVars();
			set<int>::const_iterator itVar = curSubprob.begin();
			vector<val_t>::const_iterator itVal = cur->getOptAssig().begin();

			for(; itVar!= curSubprob.end(); ++itVar, ++itVal ) {
				if (*itVal != UNKNOWN)
				assig[endVarMap[*itVar]] = *itVal;
			}

			// clear optimal assignment of AND node, since now propagated upwards
			if (cur->getType() == NODE_AND)
			cur->clearOptAssig();// TODO correct ?
		}
	} // end for

}

#endif // NO_ASSIGNMENT

#ifndef NO_ASSIGNMENT
void BoundPropagator::updateSolution(double timestamp, double cost,
		std::vector<val_t>& sol,
		std::pair<size_t, size_t> nodesMap, std::pair<size_t, size_t> nodesAll) {

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

	if ( (ISNAN(m_lowerBound) || cost > m_lowerBound) ) {
		m_lowerBound = cost;
		std::cout << "[" << std::setw(9) << timestamp << "] u "
			<< std::setw(8)  << -1 << " "
			<< std::setw(12) << nodesMap.first << " "
			<< std::setw(12) << nodesMap.second << " "
			<< std::setw(12) << nodesAll.first << " "
			<< std::setw(12) << nodesAll.second << " "
			<< std::setw(12) << m_lowerBound << std::endl;
	} else {
		return;
	}
}
#else
void BoundPropagator::updateSolution(double timestamp, double cost,
		std::pair<size_t, size_t> nodesMap, std::pair<size_t, size_t> nodesAll) {

	if ( (ISNAN(m_lowerBound) || cost > m_lowerBound) ) {
		m_lowerBound = cost;
		std::cout << "[" << std::setw(9) << timestamp << "] u "
			<< std::setw(8)  << -1 << " "
			<< std::setw(12) << nodesMap.first << " "
			<< std::setw(12) << nodesMap.second << " "
			<< std::setw(12) << nodesAll.first << " "
			<< std::setw(12) << nodesAll.second << " "
			<< std::setw(12) << m_lowerBound << " (" << -(ELEM_ENCODE(m_lowerBound)) << ")"
			<< std::endl;
	} else {
		return;
	}
}
#endif

//void BoundPropagator::updateSolutionLog(double timestamp, double cost,
//		double lowerBound, std::pair<size_t, size_t> nodesMap) {
//
//	if ( (ISNAN(m_upperBound) || cost < m_upperBound) ) {
//		m_upperBound = cost;
////		std::cout << "[" << std::setw(9) << timestamp << "] u "
////			<< std::setw(8)  << -1 << " "
////			<< std::setw(12) << nodesMap.first << " "
////			<< std::setw(12) << nodesMap.second << " "
////			<< std::setw(12) << m_upperBound << " (" << 1/(ELEM_DECODE(m_upperBound)) << ") "
////			<< std::setw(12) << lowerBound << " (" << 1/(ELEM_DECODE(lowerBound)) << ") "
////			<< std::endl;
//	} else {
//		return;
//	}
//}
