/*
 * mini_bucket_tree_incr_heur.h
 *
 *  Created on: Oct 16, 2013
 *      Author: radu
 */

#ifndef MINI_CLUSTER_TREE_INCR_HEUR_H_
#define MINI_CLUSTER_TREE_INCR_HEUR_H_

#include "heuristic.h"
#include "join_tree.h"
#include "program_options.h"

class MiniClusterTreeIncrHeur : public Heuristic {
protected:

	// The ibound for this MCTE instance
	int m_ibound;

	// the bucket tree
	scoped_ptr<JTree> m_joinTree;

	// keeps track of the changes during depth-first search
	std::stack<std::list<Change> > m_changes;

public:

	// not used
	size_t limitSize(size_t memlimit, const std::vector<int> * assignment) {
		assert(false);
		return 0;
	}

	// not used here
	double getHeur(int var, const std::vector<int>& assignment) const {
		assert(false); // not used yet
		return ELEM_NAN;
	}

	// computes the heuristic for variable var given a (partial) assignment
	double getHeur(int var, std::vector<double>& bounds,
			const std::vector<int>& assignment) const;

	// update the heuristic (only for dynamic ones)
	void update(const std::vector<int>& assignment, const bool full = true) {
		assert(false); // not used
	}

	// incremental update of the heuristic
	void updateIncr(int currVar, int nextVar, const std::vector<int>& assignment);

	// builds the bucket tree and run the initial message propagation
	size_t build(const std::vector<int>* assignment = NULL, bool computeTables =
			true);

	// returns the global upper bound
	double getGlobalUB() const {
		return m_joinTree->upperBound();
	}

	// gets sum of tables sizes
	size_t getSize() const {
		assert(false);
		return 0;
	}

	// write to a file (not used)
	bool writeToFile(std::string fn) const {
		assert(false);
		return false;
	}

	// read from a file (not used)
	bool readFromFile(std::string fn) {
		assert(false);
		return false;
	}

	// returns the MAP search order (static)
	void getStaticSearchOrder(std::list<int>& M) {
		M = m_joinTree->getSearchOrder();
	}

	// returns the MAP search clusters (static)
	void getStaticSearchClusters(std::list<int>& C) {
		C.clear();
		const std::vector<JTreeNode*>& clusters = m_joinTree->getSearchClusters();
		for (std::vector<JTreeNode*>::const_iterator it = clusters.begin();
				it != clusters.end(); ++it) {
			C.push_back( (*it)->id() );
		}
	}

	// rollback the last changes
	void rollback();

	void dump() {};

	void check(int) {};

public:

	// constructor
	MiniClusterTreeIncrHeur(Problem* p, Pseudotree* pt, ProgramOptions* po);

	// destructor
	virtual ~MiniClusterTreeIncrHeur();


};

/* inline definitions */
inline MiniClusterTreeIncrHeur::MiniClusterTreeIncrHeur(Problem* p, Pseudotree* pt, ProgramOptions* po)
	: Heuristic(p, pt, po) {

	m_ibound = po->ibound;
}

inline MiniClusterTreeIncrHeur::~MiniClusterTreeIncrHeur() {

}



#endif /* MINI_BUCKET_TREE_INCR_HEUR_H_ */
