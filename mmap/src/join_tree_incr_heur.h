/*
 * bucket_tree_inc_heur.h
 *
 *  Created on: Oct 1, 2013
 *      Author: radu
 */

#ifndef IBM_MAP_JOIN_TREE_INCR_HEUR2_H_
#define IBM_MAP_JOIN_TREE_INCR_HEUR_H_

#include "heuristic.h"
#include "join_tree.h"

/**
 * Incremental join-tree heuristic
 *
 *
 */

class JoinTreeIncrHeur : public Heuristic {
protected:

	// the join tree
	scoped_ptr<JTree> m_joinTree;

	// the stack with changes
	std::stack<std::list<Change> > m_changes;

protected:

	// reset the bucket tree (internal use only)
	void reset();

public:

	// not used
	size_t limitSize(size_t memlimit, const std::vector<int> * assignment) {
		assert(false);
		return 0;
	}

	// used by MBE
	double getHeur(int, const std::vector<int>&) const {
		assert(false);
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
//		return m_joinTree->getMapAssignmentValue();
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

	/* Returns the static ordering of the MAP variables and clusters. This
	 * is used by the incremental bucket/join tree heuristics.
	 */
	void getStaticSearchOrder(std::list<int>& M) {
		M = m_joinTree->getSearchOrder();
	}

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

	void dump(){};

	void check(int){};

public:

	// constructor
	JoinTreeIncrHeur(Problem* p, Pseudotree* pt, ProgramOptions* po);

	// destructor
	virtual ~JoinTreeIncrHeur();

};

/* inline definitions */
inline JoinTreeIncrHeur::JoinTreeIncrHeur(Problem* p, Pseudotree* pt, ProgramOptions* po) :
		Heuristic(p, pt, po) {

}

inline JoinTreeIncrHeur::~JoinTreeIncrHeur() {

}

#endif /* BUCKET_TREE_INCR_HEUR_H_ */
