/*
 * problem.h
 *
 *  Created on: Mar 28, 2013
 *      Author: radu
 */

#ifndef IBM_ANYTIME_PROBLEM_H_
#define IBM_ANYTIME_PROBLEM_H_

#include "base.h"
#include "function.h"

/* holds a problem instance with variable domains and function tables */
class Problem {

protected:

    // is last variable a dummy variable?
	bool m_hasDummy;

	// checks if dummy functions added to the model
	bool m_hasDummyFuns;

    // No. of variables
	int m_n;

    // No. of variables (before evidence was removed)
	int m_nOrg;

    // Max. domain size
	int m_k;

	// No. of evidence
	int m_e;

    // No. of functions
	int m_c;

    // Max. function arity
	int m_r;

	// Global constant modifier for objective function
	double m_globalConstant;

     // Problem name
	std::string m_name;

    // Domain sizes of variables
	std::vector<int> m_domains;

	// Types of variables (includes dummy as well)
	std::vector<var_t> m_vartypes;

	// List of functions
	std::vector<Function*> m_functions;

    // List of evidence as <index,value>
	std::map<int, int> m_evidence;

	 // Translation of variable names after removing evidence
	std::map<int, int> m_old2new;

	// the dummy variable
	int m_dummyVar;

	// determinism
	double m_determinism;

	// Query variables (ie marginal MAP)
	std::set<int> m_query;

public:

	int getDomainSize(int i) const;
	double getGlobalConstant() const;

	int getN() const {
		return m_n;
	}
	int getNOrg() const {
		return m_nOrg;
	}
	int getK() const {
		return m_k;
	}
	int getE() const {
		return m_e;
	}
	int getC() const {
		return m_c;
	}
	int getR() const {
		return m_r;
	}
	const std::string& getName() const {
		return m_name;
	}
	const map<int, int>& getEvidence() const {
		return m_evidence;
	}
	const map<int, int>& getOld2New() const {
		return m_old2new;
	}
	const set<int>& getQuery() const {
		return m_query;
	}
	const std::vector<Function*>& getFunctions() const {
		return m_functions;
	}
	const std::vector<int>& getDomains() const {
		return m_domains;
	}
	const std::vector<var_t>& getVartypes() const {
		return m_vartypes;
	}
	bool hasDummy() const {
		return m_hasDummy;
	}
	bool hasDummyFuns() const {
		return m_hasDummyFuns;
	}
	double getDeterminism() const {
		return m_determinism;
	}

public:

	/* parses a UAI format input file */
	bool parseUAI(const std::string& prob, const std::string& evid, const bool positive = false);

	/* write to UAI format for input files */
	void writeUAI(const string& prob) const;

	/* parses a ERG format input file */
	bool parseERG(const std::string& prob, const std::string& evid, const bool positive = false);

	/* parses an ordering from file 'file' and stores it in 'elim' */
	bool parseOrdering(const std::string& file, std::vector<int>& elim) const;
	bool loadOrdering(const string& file, std::vector<int>& order) const;

	/* parses the MAP variables from file 'file' and stores it in 'vtypes' */
	bool parseMapVars(const std::string& file, std::vector<int>& map);
	bool loadMapVars(const string& file, std::vector<int>& mapVars);

	/* removes evidence and unary-domain variables */
	void removeEvidence();

	/* returns true iff the index variable from the full set has been eliminated
	 * as evidence or unary */
	bool isEliminated(int i) const;

	/* returns true iff i is a MAP variable */
	bool isMap(int i) const;

	/* returns true iff i is a SUM variable */
	bool isSum(int i) const;

	/* updates the variable types */
	void updateVartypes(const std::vector<int>& mapVars);

	/* sets the type of a variable */
	void setVartype(int var, var_t t);

	/* adds the dummy variable to connect disconnected pseudo tree components */
	void addDummy();

	/* indicates that dummy functions have been added to the model */
	void addDummyFuns();

	/* returns the dummy variable */
	int getDummyVar();

	/* add a new function to the problem */
	void addFunction(Function* f);

	void replaceFunctions(const vector<Function*>& newFunctions);

	/* clone the problem instance */
	Problem* clone();
	// create a new problem conditioned on a partial assignment
	Problem* conditioned(const std::vector<int>& assignment);
	// set evidence
	void setEvidence(const std::vector<int>& assignment);

public:

	Problem();
	virtual ~Problem();
};

/* Inline definitions */

inline void Problem::addFunction(Function* f) {
	m_functions.push_back(f);
	m_c = (int) m_functions.size();
}

inline int Problem::getDomainSize(int i) const {
	assert(i<m_n);
	return m_domains[i];
}

inline double Problem::getGlobalConstant() const {
	return m_globalConstant;
}

inline void Problem::addDummy() {
	m_n += 1;
	m_hasDummy = true;
	m_domains.push_back(1); // unary domain
	m_vartypes.push_back(VAR_MAX);
	m_dummyVar = m_n - 1; // the dummy variable
}

inline void Problem::addDummyFuns() {
	m_hasDummyFuns = true;
}

inline int Problem::getDummyVar() {
	return m_dummyVar;
}

inline Problem::Problem() :
		m_hasDummy(false), m_hasDummyFuns(false), m_n(UNKNOWN),
		m_nOrg(UNKNOWN), m_k(UNKNOWN),
		m_e(UNKNOWN), m_c(UNKNOWN), m_r(UNKNOWN), m_globalConstant(0),
		m_dummyVar(UNKNOWN), m_determinism(0.0) {
	/* empty*/
}

inline Problem::~Problem() {
	// delete functions
	for (std::vector<Function*>::iterator it = m_functions.begin();
			it != m_functions.end(); ++it) {
		if (*it)
			delete (*it);
	}
}

#endif /* PROBLEM_H_ */
