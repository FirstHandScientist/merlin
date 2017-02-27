/*
 * mini_bucket_elim.h
 *
 *  Created on: Mar 28, 2013
 *      Author: radu
 */

#ifndef IBM_ANYTIME_MINI_BUCKET_ELIM_H_
#define IBM_ANYTIME_MINI_BUCKET_ELIM_H_

#include "heuristic.h"
#include "function.h"
#include "problem.h"
#include "program_options.h"
#include "pseudo_tree.h"
#include "mini_bucket.h"

/* The overall minibucket elimination */
class MiniBucketHeur: public Heuristic {

	friend class MiniBucket;

protected:
	// The ibound for this MB instance
	int m_ibound;

	// The global upper bound (MAP)
	double m_globalUB;

	// The augmented buckets that will store the minibucket functions (but not the original ones)
	std::vector<std::list<Function*> > m_augmented;
	// Precompute and store, for each variable v, the relevant intermediate functions that are
	// generated in a pseudotree descendant and passed to an ancestor of v
	// (points to the same function objects as m_augmented)
	std::vector<std::list<Function*> > m_intermediate;

protected:

	// Computes a dfs order of the pseudo tree, for building the bucket structure
	void findDfsOrder(std::vector<int>&) const;

	// Compares the size of the scope of two functions
    //bool scopeIsLarger(Function*, Function*) const;

	// reset the data structures
	void reset();

public:

	// update the heuristic give a partial assignment (not used)
	void update(const std::vector<int>&, const bool full = true) {
		assert(false);
	}

	// incremental update of the heuristic
	void updateIncr(int currVar, int nextVar,
			const std::vector<int>& assignment) {
		assert(false);
	}

	// checks if the given i-bound would exceed the memlimit and lowers
	// it accordingly.
	size_t limitSize(size_t memlimit, const std::vector<int> * assignment);

	// builds the heuristic, limited to the relevant subproblem, if applicable.
	// if computeTables=false, only returns size estimate (no tables computed)
	size_t build(const std::vector<int>* assignment = NULL, bool computeTables =
			true);

	// returns the global upper bound
	double getGlobalUB() const {
		return m_globalUB;
	}

	// computes the heuristic for variable var given a (partial) assignment
	double getHeur(int var, const std::vector<int>& assignment) const;

	// computes the heuristic for variable var given a (partial) assignment
	double getHeur(int var, std::vector<double>& bounds,
			const std::vector<int>& assignment) const {
		assert(false);
		return ELEM_NAN;
	}

	// reset the i-bound
	void setIbound(int ibound) {
		m_ibound = ibound;
	}
	// gets the i-bound
	int getIbound() const {
		return m_ibound;
	}

	// gets sum of tables sizes
	size_t getSize() const;

	bool writeToFile(std::string fn) const;
	bool readFromFile(std::string fn);

	bool isAccurate();

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
	void check(int) {}
public:
	MiniBucketHeur(Problem* p, Pseudotree* pt, ProgramOptions* po, int ib);
	virtual ~MiniBucketHeur();

};

/* Inline definitions */

inline bool MiniBucketHeur::isAccurate() {
	assert(m_pseudotree);
	return (m_pseudotree->getWidth() < m_ibound);
}

inline MiniBucketHeur::MiniBucketHeur(Problem* p, Pseudotree* pt,
		ProgramOptions* po, int ib) :
		Heuristic(p, pt, po), m_ibound(ib), m_globalUB(ELEM_ONE)
{
}

inline MiniBucketHeur::~MiniBucketHeur() {
	// make sure to delete each function only once
	for (std::vector<std::list<Function*> >::iterator it = m_augmented.begin();
			it != m_augmented.end(); ++it)
		for (std::list<Function*>::iterator it2 = it->begin(); it2 != it->end();
				++it2)
			delete (*it2);
}

#endif /* MINI_BUCKET_ELIM_H_ */
