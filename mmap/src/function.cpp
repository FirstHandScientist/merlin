/*
 * function.cpp
 *
 *  Created on: Mar 29, 2013
 *      Author: radu
 */

#include "function.h"
#include "problem.h"

/* Constructor */
Function::Function(const int& id, Problem* p, const set<int>& scope, double* T,
		const size_t& size) :
		m_id(id), m_problem(p), m_table(T), m_tableSize(size), m_scope(scope),
		m_tightness(0), m_original(false) {

	m_offsets.resize(scope.size());
	size_t offset = 1;
	set<int>::reverse_iterator rit;
	vector<size_t>::reverse_iterator ritOff;
	for (rit = m_scope.rbegin(), ritOff = m_offsets.rbegin();
			rit != m_scope.rend(); ++rit, ++ritOff) {
		*ritOff = offset;
		offset *= m_problem->getDomainSize(*rit);
	}
}

/* ATI: convert between mex/Factor representation and daoopt Function representation */
mex::Factor Function::asFactor() {
	mex::VarSet vs;
	for (set<int>::iterator it = m_scope.begin(); it != m_scope.end(); ++it)
		vs += mex::Var(*it, m_problem->getDomainSize(*it));
	mex::Factor F(vs, 0.0); //m_table);
	mex::vector<mex::Var> ord(vs.begin(), vs.end());
	mex::permuteIndex pi(ord, true);
	for (size_t j = 0; j < F.numel(); j++)
		F[j] = m_table[pi.convert(j)]; //F[pi.convertLinearIndex(j)]=m_table[j];
	return F;
}
void Function::fromFactor(const mex::Factor& F) {
	mex::VarSet vs;
	for (set<int>::iterator it = m_scope.begin(); it != m_scope.end(); ++it)
		vs += mex::Var(*it, m_problem->getDomainSize(*it));
	assert( vs == F.vars());

	mex::vector<mex::Var> ord(F.vars().begin(), F.vars().end());
	mex::permuteIndex pi(ord, true);
	for (size_t j = 0; j < F.numel(); j++)
		m_table[pi.convert(j)] = F[j]; //m_table[j] = F[pi.convertLinearIndex(j)];
}

/* returns the table entry for the assignment (input is vector of int) */
double Function::getValue(const std::vector<int>& assignment) const {
//  assert(isInstantiated(assignment)); // make sure scope is fully instantiated
#ifdef PRECOMP_OFFSETS
	size_t idx = 0;
	set<int>::const_iterator it=m_scope.begin();
	vector<size_t>::const_iterator itOff=m_offsets.begin();
	for (; it!=m_scope.end(); ++it, ++itOff)
	idx += assignment[*it] * (*itOff);
#else
	size_t idx = 0, offset = 1;
	for (set<int>::reverse_iterator rit = m_scope.rbegin();
			rit != m_scope.rend(); ++rit) {
		idx += assignment[*rit] * offset;
		offset *= m_problem->getDomainSize(*rit);
	}
#endif
	assert(idx < m_tableSize);
	return m_table[idx];
}

/* set the table entry for the assignment (input is vector of int) */
void Function::setValue(double z, const std::vector<int>& assignment) {
//  assert(isInstantiated(assignment)); // make sure scope is fully instantiated
#ifdef PRECOMP_OFFSETS
	size_t idx = 0;
	set<int>::const_iterator it=m_scope.begin();
	vector<size_t>::const_iterator itOff=m_offsets.begin();
	for (; it!=m_scope.end(); ++it, ++itOff)
	idx += assignment[*it] * (*itOff);
#else
	size_t idx = 0, offset = 1;
	for (set<int>::reverse_iterator rit = m_scope.rbegin();
			rit != m_scope.rend(); ++rit) {
		idx += assignment[*rit] * offset;
		offset *= m_problem->getDomainSize(*rit);
	}
#endif
	assert(idx < m_tableSize);
	m_table[idx] = z;
}

/* returns the table entry for the assignment (input is vector of pointers to int) */
double Function::getValuePtr(const vector<int*>& tuple) const {
	assert(tuple.size() == m_scope.size());
	// make sure tuple size matches scope size
#ifdef PRECOMP_OFFSETS
	size_t idx = 0;
	for (size_t i=0; i<m_scope.size(); ++i)
	idx += *(tuple[i]) * m_offsets[i];
#else
	size_t idx = 0, offset = 1;
	set<int>::reverse_iterator rit = m_scope.rbegin();
	vector<int*>::const_reverse_iterator ritTup = tuple.rbegin();
	for (; rit != m_scope.rend(); ++rit, ++ritTup) {
		idx += *(*ritTup) * offset;
		offset *= m_problem->getDomainSize(*rit);
	}
#endif
	assert(idx < m_tableSize);
	return m_table[idx];
}

void Function::translateScope(const map<int, int>& translate) {
	//cout << "Translating " << m_id << " :";
	set<int> newScope;
	for (set<int>::iterator it = m_scope.begin(); it != m_scope.end(); ++it) {
		newScope.insert(translate.find(*it)->second);
		//cout << ' ' << *it << "->" << translate.find(*it)->second;
	}
	//cout << endl;
	m_scope = newScope;
}

/* Main substitution function.
 * !! Changes the values of the newScope, newTable, newTableSize reference arguments !!
 */
void Function::substitute_main(const map<int, int>& assignment,
		set<int>& newScope, double*& newTable, size_t& newTableSize) const {

	map<int, int> localElim; // variables in scope to be substituted
	newTableSize = 1; // size of resulting table

	for (set<int>::const_iterator it = m_scope.begin(); it != m_scope.end();
			++it) {
		// assignment contains evidence and unary variables
		map<int, int>::const_iterator s = assignment.find(*it);
		if (s != assignment.end()) {
			localElim.insert(*s); // variable will be instantiated
		} else {
			newScope.insert(*it);
			newTableSize *= m_problem->getDomainSize(*it);
		}
	}

	if (newTableSize == m_tableSize) {
		// no reduction, just return the original table
		newTable = new double[newTableSize];
		for (size_t i = 0; i < newTableSize; ++i)
			newTable[i] = m_table[i];
		return;
	}

	// Collect domain sizes of variables in scope and
	// compute resulting table size
	vector<int> domains(m_scope.size(), UNKNOWN);
	int i = 0;
	for (set<int>::iterator it = m_scope.begin(); it != m_scope.end(); ++it) {
		domains[i] = m_problem->getDomainSize(*it);
		++i;
	}

	// initialize the new table
	newTable = new double[newTableSize];

	// keeps track of assignment while going over table
	vector<int> tuple(m_scope.size(), 0);
	size_t idx = 0, newIdx = 0;
	do {
		// does it match?
		bool match = true;
		set<int>::iterator itSco = m_scope.begin();
		vector<int>::iterator itTup = tuple.begin();
		map<int, int>::iterator itEli = localElim.begin();
		for (; match && itEli != localElim.end(); ++itEli) {
			while (itEli->first != *itSco) {
				++itSco;
				++itTup;
			}
			if (*itTup != itEli->second)
				match = false; // current tuple does not match evidence
		}

		// tuple is compatible with evidence, record in new table
		if (match) {
			newTable[newIdx] = m_table[idx];
			++newIdx;
		}

	} while (increaseTuple(idx, tuple, domains)); // increases idx and the tuple

}

Function* FunctionBayes::clone() const {
	double* newTable = new double[m_tableSize];
	for (size_t i = 0; i < m_tableSize; ++i)
		newTable[i] = m_table[i];
	Function* f = new FunctionBayes(m_id, m_problem, m_scope, newTable,
			m_tableSize);
	return f;
}

Function* FunctionBayes::substitute(const map<int, int>& assignment) const {

	set<int> newScope;
	double* newTable = NULL;
	size_t newTableSize;
	// compute the new scope, table, and table size
	substitute_main(assignment, newScope, newTable, newTableSize);

	// Create new, modified function object
	Function* newF = new FunctionBayes(m_id, m_problem, newScope, newTable,
			newTableSize);
	return newF;

}


