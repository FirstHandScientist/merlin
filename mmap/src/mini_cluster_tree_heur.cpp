/*
 * mini_bucket_tree_heur.cpp
 *
 *  Created on: Oct 23, 2013
 *      Author: radu
 */

#include "mini_cluster_tree_heur.h"

// returns the heuristic value of a variable (given a partial assignment)
double MiniClusterTreeHeur::getHeur(int var,
		const std::vector<int>& assignment) const {
	assert(false); // not used yet
	return ELEM_NAN;
}

// returns the heuristics associated to the domain values of an unassigned variable
double MiniClusterTreeHeur::getHeur(int var,
		std::vector<double>& bounds,
		const std::vector<int>& assignment) const {

	// safety checks
	assert( var >= 0 && var < m_problem->getN());
	assert( m_problem->isMap(var) );
	assert( assignment[var] == NONE );

	bounds.clear();
	JTreeNode* node = m_joinTree->getSearchCluster(var);

	mex::Factor M = node->marginal(var, m_ibound, assignment);
	for (size_t i = 0; i < M.numel(); ++i) {
		bounds.push_back(M[i]);
	}

	return M.max();
}

// update of the heuristic
void MiniClusterTreeHeur::update(const std::vector<int>& assignment, const bool full) {
	m_joinTree->propagate(m_ibound, assignment);
}

// builds the mini-bucket tree and run the initial message propagation
size_t MiniClusterTreeHeur::build(const std::vector<int>* assignment,
		bool computeTables) {

	// safety checks
	assert(m_problem->hasDummy());
	int N = ( m_problem->hasDummyFuns() ? m_problem->getN() : m_problem->getN()-1 );

	// create primal graph (without the dummy functions, if any)
	Graph g(N-1);
	std::vector<Function*> dummyFuns;
	const std::vector<Function*>& fns = m_problem->getFunctions();
	for (std::vector<Function*>::const_iterator it = fns.begin();
			it != fns.end(); ++it) {
		Function* f = (*it);
		if (f->getScope().find(m_problem->getDummyVar()) == f->getScope().end()) {
			g.addClique((*it)->getScope());
		} else {
			dummyFuns.push_back(f);
		}
	}

	// get the unconstrained elimination order (without the dummy variable)
	std::vector<int> elimOrder;
	int w = g.eliminate(elimOrder);
	cout << "Computed unconstrained elimination ordering "
			<< " (" << w << '/' << N << ")." << endl;

	// add the dummy functions scopes to the graph (if any)
	if (m_problem->hasDummyFuns()) {
		elimOrder.push_back(m_problem->getDummyVar());
		for (std::vector<Function*>::iterator it = dummyFuns.begin();
				it != dummyFuns.end(); ++it) {
			g.addClique((*it)->getScope());
		}
	}

	// build the bucket tree
	m_joinTree.reset(new JTree(m_problem, m_options));
	m_joinTree->build(g, elimOrder, true);

	// run the initial message propagation and report the initial upper bound
	std::vector<int> evidence;
	evidence.resize(N, NONE);
	m_joinTree->propagate(m_ibound, evidence);

	cout << "\tUpper Bound = " << m_joinTree->upperBound() << endl;

	return 0;
}
