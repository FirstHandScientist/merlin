/*
 * any_wraobf.cpp
 *
 *  Created on: Oct 29, 2014
 *      Author: radu
 */

#include "any_wraobf2.h"

// repair the search graph with the new weight - traverse the graph bottom up,
// starting from the current tip nodes and update node values using the new
// weighted h-value; redo the markings and solved labels as well
// return TRUE if new partial solution tree found; return FALSE if no new
// partial solution tree found
void AnyWRAOBF2::repair() {

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
				}

				node->setRepaired(true);
			}
		}
	}

	// find the best partial solution tree
	m_tips.clear();
	findBestPartialTree();
}

// search routine
int AnyWRAOBF2::aostar(bool verbose) {

	// initialize the weight
	m_epsilon = m_schedule->getInitialWeight();
	assert(m_epsilon >= 1);

	// create the root node of the search space
	int numLayers = 2*(m_pseudotree->getHeight() + 2);

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
	m_space->initLayers(numLayers);
	m_space->setRoot(root);
	m_space->add(state, root);
	m_tips.push_back(root);

	int result = SEARCH_SUCCESS;
	bool timeout = false;
	while (m_epsilon >= 1) {

		// perform a weighted AO* search
		while ( !root->isSolved() ) {

			// check for time limit violation
			if (m_timer.timeout()) {
				timeout = true;
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
						<< std::setw(12) << m_upperBound << " (" << 1/ELEM_DECODE(m_upperBound) << ")"
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
						<< std::setw(12) << m_upperBound << " (" << 1/ELEM_DECODE(m_upperBound) << ")"
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

		// check for timeout
		if (timeout) {
			result = SEARCH_TIMEOUT;
			break;
		}

		// get the current (sub)optimal solution
		if (root->isSolved()) {
			m_solutionCost = root->getValue();
			m_solutionCost += -ELEM_ENCODE(m_problem->getGlobalConstant());
			double timestamp = m_timer.elapsed();
			if ( ISNAN(m_upperBound) ) {
				m_upperBound = m_solutionCost;
			} else {
				m_upperBound = std::min(m_upperBound, m_solutionCost);
			}

			// output solution
			std::cout << "["
				<< std::setw(9) << timestamp << "] w "
				<< std::setw(8) << m_epsilon << " "
				<< std::setw(12) << m_space->getOrNodes() << " "
				<< std::setw(12) << m_space->getAndNodes() << " "
				<< std::setw(12) << m_upperBound << " (" << 1/ELEM_DECODE(m_upperBound) << ")"
				<< std::endl;
		}

		if (m_epsilon == 1.0) { // proved optimality
			m_solved = true;
			break;
		}

		// update the weight
		m_epsilon = m_schedule->getNextWeight();

		// repair the search space
		repair();
	}

	return result;
}




