/*
 * bucket_tree_incr_heur.cpp
 *
 *  Created on: Oct 1, 2013
 *      Author: radu
 */


#include "join_tree_incr_heur.h"

// builds the bucket tree and run the initial message propagation
size_t JoinTreeIncrHeur::build(const std::vector<int>* assignment,
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
	cout << "Computed unconstrained (minfill) elimination ordering "
			<< " (" << w << '/' << N << ")." << endl;

	// DEBUG only
//	m_problem->loadOrdering(m_options->orderingFile, elimOrder);

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

// incremental update of the heuristic (M should be the unassigned variables),
// which follows the assigning of the current variable (currVar)
void JoinTreeIncrHeur::updateIncr(int currVar, int nextVar,
		const std::vector<int>& assignment) {

	if (nextVar == NONE) {
		std::list<Change> noChanges;
		m_changes.push(noChanges);

		JTreeNode* from = m_joinTree->getSearchCluster( currVar );
		double upbo = from->upperBoundIncr(assignment);
		m_joinTree->setUpperBound(upbo);

		return; // no propagation; all MAP variables instantiated
	}

	// safety checks
	assert(assignment[ currVar ] != NONE );

	// partial message propagation is triggered only if the current search
	// cluster is fully assigned (ie, all its MAP variables have values)
	JTreeNode* from = m_joinTree->getSearchCluster( currVar );
	bool assigned = true;
	for (std::set<int>::const_iterator si = from->varsMax().begin();
			si != from->varsMax().end(); ++si) {
		int var = (*si);
		if (assignment[var] == NONE) {
			assigned = false;
			break;
		}
	}

	if (!assigned) { // current cluster still has unassigned MAP variables
		std::list<Change> noChanges;
		m_changes.push(noChanges);

		return; // no propagation
	}

	// target cluster found; attempt propagation on the pre-computed path
	assert( assignment[nextVar] == NONE );
	JTreeNode* to = m_joinTree->getSearchCluster(nextVar);
	assert( from != to );

	// save messages along the path between 'from' and 'to' clusters
	std::list<Change> changes;
	std::vector<PathEdge*>& path = m_joinTree->getPath( from->id() );
	for (std::vector<PathEdge*>::iterator it = path.begin();
			it != path.end(); ++it) {
		PathEdge* e = (*it);
		if (e->activeMessage == 1) {
			changes.push_back(Change(e->edge, 1, e->edge->getMessage1()));
		} else {
			changes.push_back(Change(e->edge, 2, e->edge->getMessage2()));
		}
	}

	// record changes; changes are associated with a variable-value assignment
	m_changes.push(changes);

	// update the messages
	for (std::vector<PathEdge*>::iterator it = path.begin();
			it != path.end(); ++it) {
		PathEdge* e = (*it);
		e->updateMessage(assignment);
	}

}

// rollbacks the last changes on the stack
void JoinTreeIncrHeur::rollback() {
	std::list<Change>& changes = m_changes.top();
	for (std::list<Change>::reverse_iterator it = changes.rbegin();
			it != changes.rend(); ++it) {
		Change& ch = (*it);
		JTreeEdge* e = ch.edge;
		if (ch.message == 1) e->setMessage1(ch.potential);
		else e->setMessage2(ch.potential);
	}
	m_changes.pop();
}

// returns the heuristics associated to the domain values of an unassigned variable
double JoinTreeIncrHeur::getHeur(int var,
		std::vector<double>& bounds,
		const std::vector<int>& assignment) const {

	// safety checks
	assert( var >= 0 && var < m_problem->getN());
	assert( m_problem->isMap(var) );
	assert( assignment[var] == NONE );

	bounds.clear();
	JTreeNode* node = m_joinTree->getSearchCluster(var);

	mex::Factor M = node->marginalIncr(var, assignment);
	for (size_t i = 0; i < M.numel(); ++i) {
		bounds.push_back(M[i]);
	}

	return M.max();
}

