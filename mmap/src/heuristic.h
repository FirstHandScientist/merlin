/*
 * heuristic.h
 *
 *  Created on: Mar 28, 2013
 *      Author: radu
 */

#ifndef IBM_ANYTIME_HEURISTIC_H_
#define IBM_ANYTIME_HEURISTIC_H_

#include "base.h"

/* forward class declarations */
class ProgramOptions;
class Problem;
class Pseudotree;

/*
 * Base class for all heuristic implementations.
 * The following are the three central functions that need to be implemented
 * in subclasses (other functions can have empty implementations):
 * - build(...) : to actually build the heuristic.
 * - getHeur(...) : to query heuristic values
 * - getGlobalUB() : to retrieve the problem upper bound
 */
class Heuristic {

protected:
	Problem* m_problem;            // The problem instance
	Pseudotree* m_pseudotree;      // The underlying pseudotree
	ProgramOptions* m_options;     // Program options instance

public:

	/**
	 * Updates the heuristic given a partial assignment (only for dynamic ones)
	 */
	virtual void update(const std::vector<int>& assignment,
			const bool full = true) = 0;

	// incremental update of the heuristic
	virtual void updateIncr(int currVar, int nextVar,
			const std::vector<int>& assignment) = 0;

	/* Limits the memory size of the heuristic (e.g. through lowering
	 * the mini bucket i-bound)
	 */
	virtual size_t limitSize(size_t limit,
			const std::vector<int> * assignment) = 0;

	/* Returns the memory size of the heuristic
	 */
	virtual size_t getSize() const = 0;

	/* Allows the heuristic to apply preprocessing to the problem instance, if
	 * applicable. Returns true if any changes were made to original problem,
	 * false otherwise. Optional argument used to specify partial assignment in
	 * case of conditioned subproblem.
	 */
	virtual bool preprocess(const std::vector<int>* = NULL) {
		return false;
	}

	/* Builds the heuristic. The first optional parameter can be used to specify
	 * a partial assignment (when solving a conditioned subproblem,for instance),
	 * the second optional parameter signals simulation-only mode (i.e. the
	 * heuristic compilation is just simulated). The basic heuristic will work
	 * without these two features, though.
	 * The function should return the memory size of the heuristic instance.
	 */
	virtual size_t build(const std::vector<int>* = NULL, bool computeTables =
			true) = 0;

	/* Reads and writes heuristic from/to specified file.
	 * Should return true/false on success or failure, respectively.
	 */
	virtual bool readFromFile(std::string filename) = 0;
	virtual bool writeToFile(std::string filename) const = 0;

	/* Returns the global upper bound on the problem solution (e.g., after
	 * marginalizing out the root bucket - assuming a minimization task)
	 */
	virtual double getGlobalUB() const = 0;

	/* Compute and return the heuristic estimate for the given variable.
	 * The second argument is an assignment of all problem variables (i.e.
	 * assignment.size() = number of problem variables incl. the dummy) from
	 * from which the relevant assignments (i.e. the context of the current
	 * variable) will be extracted.
	 */
	virtual double getHeur(int var,
			const std::vector<int>& assignment) const = 0;

	// computes the heuristic for variable var given a (partial) assignment
	virtual double getHeur(int var, std::vector<double>& bounds,
			const std::vector<int>& assignment) const = 0;

	/* Returns true if the heuristic is provably accurate. Default false,
	 * override in child class.
	 */
	virtual bool isAccurate() {
		return false;
	}

	/* Returns the static ordering of the MAP variables and clusters. This
	 * is used by the incremental bucket/join tree heuristics.
	 */
	virtual void getStaticSearchOrder(std::list<int>&) = 0;
	virtual void getStaticSearchClusters(std::list<int>&) = 0;

	// rollback the last changes
	virtual void rollback() = 0;

	virtual void dump() = 0;

	virtual void check(int) = 0;

protected:
	Heuristic(Problem* p, Pseudotree* pt, ProgramOptions* po) :
			m_problem(p), m_pseudotree(pt), m_options(po) { /* nothing */
	}
public:
	virtual ~Heuristic() { /* nothing */
	}
};

/* "Empty" heuristic, for testing and debugging */
class UnHeuristic: public Heuristic {
public:
	void update(const std::vector<int>&) {
		assert(false);
	}
	size_t limitSize(size_t, ProgramOptions*, const std::vector<int> *) {
		return 0;
	}
	size_t getSize() const {
		return 0;
	}
	size_t build(const std::vector<int>*, bool) {
		return 0;
	}
	bool readFromFile(std::string) {
		return true;
	}
	bool writeToFile(std::string) const {
		return true;
	}
	double getGlobalLB() const {
		assert(false);
		return 0;
	}
	double getHeur(int, const std::vector<int>&) const {
		assert(false);
		return 0;
	}
	bool isAccurate() {
		return false;
	}
	UnHeuristic() :
			Heuristic(NULL, NULL, NULL) {
	}
	virtual ~UnHeuristic() {
	}
};

#endif /* HEURISTIC_H_ */
