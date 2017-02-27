/*
 * mini_bucket_tree_heur.h
 *
 *  Created on: Oct 23, 2013
 *      Author: radu
 */

#ifndef MINI_CLUSTER_TREE_HEUR_H_
#define MINI_CLUSTER_TREE_HEUR_H_

#include "heuristic.h"
#include "join_tree.h"
#include "program_options.h"


/**
 * A dynamic mini-bucket tree heuristic
 */

class MiniClusterTreeHeur : public Heuristic {

protected:

	// The ibound for this MCTE instance
	int m_ibound;

	// the join tree
	scoped_ptr<JTree> m_joinTree;

public:

	// not used
	size_t limitSize(size_t memlimit, const std::vector<int> * assignment) {
		assert(false);
		return 0;
	}

	// used by MBE and MBTE
	double getHeur(int var, const std::vector<int>& assignment) const;

	// computes the heuristic for variable var given a (partial) assignment
	double getHeur(int var, std::vector<double>& bounds,
			const std::vector<int>& assignment) const;

	// update the heuristic (only for dynamic ones)
	void update(const std::vector<int>& assignment, const bool full = true);

	// incremental update of the heuristic
	void updateIncr(int currVar, int nextVar, const std::vector<int>& assignment) {
		assert(false); // not used
	}

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

	// not used
	void getStaticSearchOrder(std::list<int>& M) {
		M = m_joinTree->getSearchOrder();
	}

	// not used
	void getStaticSearchClusters(std::list<int>& C) {
		C.clear();
		const std::vector<JTreeNode*>& clusters = m_joinTree->getSearchClusters();
		for (std::vector<JTreeNode*>::const_iterator it = clusters.begin();
				it != clusters.end(); ++it) {
			C.push_back( (*it)->id() );
		}
	}

	// rollback the last changes
	void rollback() {
		// do nothing
	}
	void dump() {
		// do nothing
	}

	void check(int) {};

public:

	// constructor
	MiniClusterTreeHeur(Problem* p, Pseudotree* pt, ProgramOptions* po);

	// destructor
	virtual ~MiniClusterTreeHeur();
};

/* inline definitions */
inline MiniClusterTreeHeur::MiniClusterTreeHeur(Problem* p, Pseudotree* pt, ProgramOptions* po)
	: Heuristic(p, pt, po) {

	m_ibound = po->ibound;
}

inline MiniClusterTreeHeur::~MiniClusterTreeHeur() {

}


#endif /* MINI_BUCKET_TREE_HEUR_H_ */
