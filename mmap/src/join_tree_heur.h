/*
 * bucket_tree_heur.h
 *
 *  Created on: Sep 26, 2013
 *      Author: radu
 */

#ifndef JOIN_TREE_HEUR_H_
#define JOIN_TREE_HEUR_H_


#include "heuristic.h"
#include "join_tree.h"

class JoinTreeHeur : public Heuristic {

protected:

	// the join tree
	scoped_ptr<JTree> m_joinTree;

protected:

	// reset the bucket tree (internal use only)
	void reset();

public:

	// not used
	size_t limitSize(size_t memlimit, const std::vector<int> * assignment) {
		assert(false);
		return 0;
	}

	// update the heuristic (only for dynamic ones)
	void update(const std::vector<int>& assignment, const bool full = true);

	// incremental update of the heuristic
	void updateIncr(int currVar, int nextVar, const std::vector<int>& assignment) {
		assert(false);
	}

	// builds the bucket tree and run the initial message propagation
	size_t build(const std::vector<int>* assignment = NULL, bool computeTables =
			true);

	// returns the global upper bound
	double getGlobalUB() const {
		return m_joinTree->upperBound();
	}

	// used by MBE only
	double getHeur(int var, const std::vector<int>& assignment) const {
		assert(false);
		return ELEM_NAN;
	}

	// computes the heuristic for variable var given a (partial) assignment
	double getHeur(int var, std::vector<double>& bounds,
			const std::vector<int>& assignment) const;

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
	void getStaticSearchOrder(std::list<int>&) {
		assert(false); // not used
	}
	void getStaticSearchClusters(std::list<int>&) {
		assert(false); // not used
	}
	void rollback() {
		assert(false); // not used
	}
	void dump() {
		assert(false);
	}
	void check(int){}

public:

	// constructor
	JoinTreeHeur(Problem* p, Pseudotree* pt, ProgramOptions* po);

	// destructor
	virtual ~JoinTreeHeur();

};

/* inline definitions */

inline JoinTreeHeur::JoinTreeHeur(Problem* p, Pseudotree* pt,
		ProgramOptions* po) : Heuristic(p, pt, po) {

}

inline JoinTreeHeur::~JoinTreeHeur() {

}


#endif /* BUCKET_TREE_HEUR_H_ */
