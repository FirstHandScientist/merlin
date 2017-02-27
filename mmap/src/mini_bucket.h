/*
 * mini_bucket.h
 *
 *  Created on: Mar 28, 2013
 *      Author: radu
 */

#ifndef IBM_ANYTIME_MINI_BUCKET_H_
#define IBM_ANYTIME_MINI_BUCKET_H_

#include "function.h"
#include "problem.h"

/* A single minibucket, i.e. a collection of functions */
class MiniBucket {

protected:
	int m_bucketVar;               // the bucket variable
	int m_ibound;                  // the ibound
	Problem* m_problem;      	   // pointer to the bucket elimination structure
	std::vector<Function*> m_functions; // the functions in the MB
	std::set<int> m_jointScope;    // keeps track of the joint scope if the functions

protected:
	// combines all MB functions into a new function over the joint scope (not assigned)
	Function* combine(const std::set<int>& jointScope, const std::vector<int>& assignment);

public:
	// checks whether the MB has space for a function
	bool allowsFunction(Function*);
	bool allowsFunction(Function*, const std::vector<int>&);
	// adds a function to the minibucket
	void addFunction(Function*);
	void addFunction(Function*, const std::vector<int>&);
	// returns the set of functions
	std::vector<Function*>& getFunctions();
	// Joins the MB functions, eliminate the bucket variable, and returns the resulting function
	// set buildTable==false to get only size estimate (table will not be computed)
	Function* eliminate(bool buildTable = true, bool first = false);
	Function* eliminate(const std::vector<int>& origElimVars,
			const std::vector<int>& evidence, const bool first = true);
public:
	MiniBucket(int var, int bound, Problem* p);
	MiniBucket(Problem* p, std::vector<Function*>& funs);

};

/* Inline definitions */

inline void MiniBucket::addFunction(Function* f) {
	assert(f);
	// insert function
	m_functions.push_back(f);
	// update joint scope
	m_jointScope.insert(f->getScope().begin(), f->getScope().end());
}

inline void MiniBucket::addFunction(Function* f, const std::vector<int>& assignment) {
	assert(f);
	// insert function
	m_functions.push_back(f);
	// update joint scope (ignore assigned variables)
	for (std::set<int>::const_iterator si = f->getScope().begin();
			si != f->getScope().end(); ++si) {
		int var = (*si);
		if (assignment[var] == NONE) {
			m_jointScope.insert(var);
		}
	}
}


inline std::vector<Function*>& MiniBucket::getFunctions() {
	return m_functions;
}

inline MiniBucket::MiniBucket(int v, int b, Problem* p) :
		m_bucketVar(v), m_ibound(b), m_problem(p) {
}

inline MiniBucket::MiniBucket(Problem* p, std::vector<Function*>& funs) :
		m_bucketVar(NONE), m_ibound(NONE), m_problem(p) {
	for (std::vector<Function*>::iterator it = funs.begin();
			it != funs.end(); ++it) {
		addFunction( *it );
	}
}

#endif /* MINI_BUCKET_H_ */
