/*
 * bucket_tree_heur.cpp
 *
 *  Created on: Sep 26, 2013
 *      Author: radu
 */


#include "join_tree_heur.h"

// computes the heuristic for variable var given a (partial) assignment
double JoinTreeHeur::getHeur(int var, std::vector<double>& bounds,
		const std::vector<int>& assignment) const {

	// safety checks
	assert( var >= 0 && var < m_problem->getN());
	assert( m_problem->isMap(var) );
	assert( assignment[var] == NONE );

	bounds.clear();
	const mex::Factor& M = m_joinTree->marginal(var);
	for (size_t i = 0; i < M.numel(); ++i) {
		bounds.push_back(M[i]);
	}

	return M.max();
}


// update the bucket tree
void JoinTreeHeur::update(const std::vector<int>& assignment, const bool full) {
	m_joinTree->propagate(assignment);
}

// reset the bucket tree heuristic
void JoinTreeHeur::reset() {

}

// build the bucket tree
size_t JoinTreeHeur::build(const std::vector<int> * assignment,
		bool computeTables) {

	// safety checks
	assert(m_problem->hasDummy());
	int N = m_problem->getN();

	// create primal graph (without the dummy variable)
	Graph g(N-1);
	std::vector<Function*> dummyFuns;
	const std::vector<Function*>& fns = m_problem->getFunctions();
	for (std::vector<Function*>::const_iterator it = fns.begin();
			it != fns.end(); ++it) {
		Function* f = (*it);
		if (f->getScope().find(m_problem->getDummyVar()) == f->getScope().end()) {
			g.addClique((*it)->getScope());
		} else {
			dummyFuns.push_back(f); // do not add a dummy function to the graph
		}
	}

	// get the unconstrained elimination order (without the dummy variable)
	std::vector<int> elimOrder;
	int w = g.eliminate(elimOrder);
	cout << "Computed unconstrained elimination ordering "
			<< " (" << w << '/' << N << ")." << endl;

	// if the model has dummy functions (ie, disconnected graph) then
	// add the dummy functions scopes to the graph and put the dummy variable
	// last in the ordering (will be the root of the bucket tree)
	if (m_problem->hasDummyFuns()) {
		elimOrder.push_back(m_problem->getDummyVar());
		for (std::vector<Function*>::iterator it = dummyFuns.begin();
				it != dummyFuns.end(); ++it) {
			g.addClique((*it)->getScope());
		}
	}

	// build the join tree
	m_joinTree.reset(new JTree(m_problem, m_options));
	m_joinTree->build(g, elimOrder);

	//m_joinTree->toDot(std::string("/home/radu/tmp/jointree.dot"));

	// run the initial message propagation and report the initial upper bound
	std::vector<int> assign;
	assign.resize(N, NONE);
	m_joinTree->propagate(assign);

	cout << "\tUpper Bound = " << m_joinTree->upperBound() << endl;

	return 0;
}
