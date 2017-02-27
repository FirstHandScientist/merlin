/*
 * function.h
 *
 *  Created on: Mar 28, 2013
 *      Author: radu
 */

#ifndef IBM_ANYTIME_FUNCTION_H_
#define IBM_ANYTIME_FUNCTION_H_

#include "base.h"
#include "utils.h"

#include "mex/Factor.h"

// Forward declaration
class Problem;

/**
 * A function in the problem as scope and table
 */

class Function {

protected:

    // Function id (number)
	int m_id;

    // Pointer to the problem
	Problem* m_problem;

    // The actual table of function values
	double* m_table;

    // Size of the table
	size_t m_tableSize;

    // Scope of the function
	std::set<int> m_scope;

	  // Precomputed offsets for value lookup
	std::vector<size_t> m_offsets;

	/* Tightness related information */
	size_t m_tightness;     // number of valid entries in table

	// Original input function
	bool m_original;

public:
	int getId() const {
		return m_id;
	}
	size_t getTableSize() const {
		return m_tableSize;
	}
	double* getTable() const {
		return m_table;
	}
	const std::set<int>& getScope() const {
		return m_scope;
	}
	int getArity() const {
		return m_scope.size();
	}
	bool isOriginal() const {
		return m_original;
	}
	/* returns true iff the function is constant */
	bool isConstant() const {
		return m_tableSize == 1;
	}

	/* true iff var. i is in scope */
	bool hasInScope(const int& i) const {
		return (m_scope.find(i) != m_scope.end());
	}

	/* true iff at least one var from S is in scope */
	bool hasInScope(const std::set<int>& S) const {
		return !intersectionEmpty(S, m_scope);
	}

	/* generates a new (smaller) function with reduced scope and the <var,val>
	 * pairs from the argument factored into the new table.
	 * ! HAS TO BE IMPLEMENTED IN SUBCLASSES !  */
	virtual Function* substitute(const std::map<int, int>& assignment) const = 0;

	/* translates the variables in the function scope */
	void translateScope(const std::map<int, int>& translate);

	/* checks if all variables in the function's scope are instantiated in assignment
	 * (note that assignment is a full problem assignment, not just the function scope) */
	bool isInstantiated(const std::vector<int>& assignment) const;

	/* returns the table entry for an assignment */
	double getValue(const std::vector<int>& assignment) const;

	/* returns the function value for the tuple (which is pointered function scope) */
	double getValuePtr(const std::vector<int*>& tuple) const;

	/* sets the table entry for an assignment */
	void setValue(double z, const std::vector<int>& assignment);

	// set the original flag to true (by default it is false)
	inline void setOriginal() {
		m_original = true;
	}

	/* ATI: convert between mex/Factor and daoopt/Function representations */
	mex::Factor asFactor();
	void fromFactor(const mex::Factor&);

protected:

	/* main work for substitution: computes new scope, new table and table size
	 * and stores them in the three non-const argument references */
	void substitute_main(const std::map<int, int>& assignment,
			std::set<int>&, double*&, size_t&) const;

public:

	/* generates a clone of this function object */
	virtual Function* clone() const = 0;

	/* gets static tightness */
	size_t getTightness() const {
		return m_tightness;
	}

public:
	virtual int getType() const = 0;

	virtual ~Function();

protected:

	Function(const int& id, Problem* p, const std::set<int>& scope, double* T,
			const size_t& size);

};

/***************************************************************************/

class FunctionBayes: public Function {

public:
	Function* substitute(const std::map<int, int>& assignment) const;
	Function* clone() const;
	inline int getType() const {
		return TYPE_BAYES;
	}

public:
	FunctionBayes(const int& id, Problem* p, const std::set<int>& scope, double* T,
			const size_t& size);
	virtual ~FunctionBayes() {
	}
};


/***************************************************************************/

/* Inline implementations */

inline Function::~Function() {
	if (m_table)
		delete[] m_table;
}

inline bool Function::isInstantiated(const std::vector<int>& assignment) const {
	for (std::set<int>::const_iterator it = m_scope.begin(); it != m_scope.end();
			++it) {
		if (assignment[*it] == UNKNOWN)
			return false;
	}
	return true;
}

inline FunctionBayes::FunctionBayes(const int& id, Problem* p,
		const set<int>& scope, double* T, const size_t& size) :
		Function(id, p, scope, T, size) {
	size_t t = 0;
	if (T) {
		for (size_t i = 0; i < size; ++i) {
			if (T[i] != ELEM_ZERO)
				++t;
		}
	}
	m_tightness = t;
}

/* cout operator */
inline std::ostream& operator <<(std::ostream& os, const std::set<int>& s) {
	os << "[ ";
	std::copy(s.begin(), s.end(), std::ostream_iterator<int>(os, " "));
	os << "]";
	return os;
}

/* cout operator */
inline std::ostream& operator <<(std::ostream& os, const Function& f) {
	os << 'f' << f.getId() << ':' << f.getScope() << "=( ";
	for (size_t i = 0; i < f.getTableSize(); ++i) {
		os << f.getTable()[i] << " ";
	}
	os << ")";
	return os;
}


inline bool scopeIsLarger(Function* p, Function* q) {
	assert(p && q);
	if (p->getArity() == q->getArity())
		return (p->getId() > q->getId());
	else
		return (p->getArity() > q->getArity());
}

inline bool scopeIsLargerSet(const std::pair<size_t, size_t>& p,
		const std::pair<size_t, size_t>& q) {

	return p.second > q.second;
}

#endif /* FUNCTION_H_ */
