/*
 * mini_bucket_mplp_heur.h
 *
 *  Created on: Oct 25, 2013
 *      Author: radu
 */

#ifndef MINI_BUCKET_WMB_JG_HEUR_H_
#define MINI_BUCKET_WMB_JG_HEUR_H_

#include "heuristic.h"
#include "function.h"
#include "problem.h"
#include "program_options.h"
#include "pseudo_tree.h"
#include "utils.h"

#include "mini_bucket.h"

#include "mex/jglp.h"

// Mini-bucket with iterative cost shifting heuristic (JGLP)

class MiniBucketJGHeur : public Heuristic {
public:

	// checks if the given i-bound would exceed the memlimit, in that
	// case compute highest possible i-bound that 'fits'
	size_t limitSize(size_t memlimit, const std::vector<val_t> *assignment) {
		return 0;
	}

	// update the heuristic give a partial assignment (not used)
	void update(const std::vector<int>&, const bool full = true) {
		assert(false);
	}

	// incremental update of the heuristic
	void updateIncr(int currVar, int nextVar,
			const std::vector<int>& assignment) {
		assert(false);
	}

	// Model pre-processing step, if any  (here, run mplp)
	bool preprocess(const std::vector<val_t>* assignment);

	// builds the heuristic, limited to the relevant subproblem, if applicable.
	// if computeTables=false, only returns size estimate (no tables computed)
	size_t build(const std::vector<val_t>* assignment = NULL, bool computeTables =
			true);

	// returns the global upper bound
	double getGlobalUB() const {
		return (ELEM_DECODE(_jglp.ub())) * _p->getGlobalConstant();
	}

	// computes the heuristic for variable var given a (partial) assignment
	double getHeur(int var, const vector<val_t>& assignment) const {
		double res = _jglp.logHeurToGo(_jglp.var(var), assignment);
		return ELEM_DECODE(res);
	}

	// computes the heuristic for variable var given a (partial) assignment
	double getHeur(int var, std::vector<double>& bounds,
			const std::vector<int>& assignment) const {
		assert(false);
		return ELEM_NAN;
	}


	// reset the i-bound
	void setIbound(int ibound) {
		_jglp.setIBound(ibound);
	}
	// gets the i-bound
	int getIbound() const {
		return _jglp.getIBound();
	}

	// gets sum of tables sizes
	size_t getSize() const {
		size_t sz = 0;
		const mex::vector<mex::Factor>& flist = _jglp.factors();
		for (size_t i = 0; i < flist.size(); ++i)
			sz += flist[i].numel();
		return sz;
	}

	bool writeToFile(string fn) const {
		return true;
	}
	bool readFromFile(string fn) {
		return false;
	}

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
	MiniBucketJGHeur(Problem* p, Pseudotree* pt, ProgramOptions* po, int ib);
	~MiniBucketJGHeur() {
	}

protected:
	Problem* _p;
	Pseudotree* _pt;
	mex::jglp _jglp;
	size_t _memlimit;
	ProgramOptions* _options;

	// Helper functions to move between DAO & mex factor formats
	mex::vector<mex::Factor> copyFactors(void);
	void rewriteFactors(const std::vector<mex::Factor>& factors);

	void findDfsOrder(vector<int>& order) const;
};

inline bool MiniBucketJGHeur::isAccurate() {
	assert(m_pseudotree);
	return (m_pseudotree->getWidth() < (int) _jglp.getIBound());
}

#endif /* MINI_BUCKET_MPLP_HEUR_H_ */
