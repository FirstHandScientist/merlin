/*
 * problem.cpp
 *
 *  Created on: Mar 29, 2013
 *      Author: radu
 */

#include "problem.h"
#include <iostream>
#include <sstream>
#include <fstream>

// Remove the evidence variable in the problem (call right after load)
void Problem::removeEvidence() {

	assert( m_n != UNKNOWN );

	// record original no. of variables
	m_nOrg = m_n;

	//map<uint,uint> evidence(m_evidence);
	std::vector<int> new_domains;
	std::vector<Function*> new_funs;

	// eliminateVar[i]==TRUE iff var. is to be eliminated
	std::vector<bool> eliminateVar(m_n, false);
	for (map<int, int>::iterator it = m_evidence.begin();
			it != m_evidence.end(); ++it) {
		eliminateVar[it->first] = true;
	}

	//std::cout << "Removing " << m_evidence.size() << " variables." << std::endl;

	// Identify and tag unary domain variables
	for (int i = 0; i < m_n; ++i) {
		if (m_domains.at(i) == 1) { // regard unary domains as evidence
			m_evidence.insert(make_pair(i, 0));
			++m_e;
			eliminateVar.at(i) = true;
		}
	}

	// Project functions to account for evidence
	m_globalConstant = ELEM_ONE;
	int new_r = 0; // max. arity
	std::vector<Function*>::iterator fi = m_functions.begin();
	for (; fi != m_functions.end(); ++fi) {
		Function *fn = (*fi);
		// Substitute evidence variables
		Function* new_fn = fn->substitute(m_evidence);
		if (new_fn->isConstant()) { // now empty scope
			m_globalConstant *= new_fn->getTable()[0];
			delete new_fn;
		} else {
			new_funs.push_back(new_fn); // record new function
			new_r = max(new_r, (int) new_fn->getScope().size());
		}
		delete fn; // delete old function

	}
	m_functions = new_funs;

	// eliminate tagged variables and reorder remaining ones
	int idx = 0, new_n = 0, new_k = 0;
	for (int i = 0; i < m_n; ++i) {
		if (!eliminateVar.at(i)) {

			m_old2new.insert(make_pair(i, idx));

			int k = m_domains.at(i);
			new_domains.push_back(k);
			new_k = max(new_k, k);

			++idx;
			++new_n;
		}
	}

	// update variable information
	m_domains = new_domains;
	m_n = new_n;
	m_k = new_k;

	// translate scopes of the new functions
	for (fi = m_functions.begin(); fi != m_functions.end(); ++fi)
		(*fi)->translateScope(m_old2new);

	// update function information
	m_c = m_functions.size();
}

// Parse and load the ordering file (required argument)
bool Problem::parseOrdering(const string& file, vector<int>& elim) const {

	// Safety checks
	assert( m_n != UNKNOWN );

	// Open the file
	std::ifstream in(file.c_str());
	if (in.fail()) { // file not existent yet
		std::cerr << "Error while opening the ordering file (file not found)." << std::endl;
		exit(EXIT_FAILURE);
	}

	// ignore first line if there's a pound '#' sign (comment)
	if (in.peek() == '#') {
		in.ignore(8192, '\n');
	}

	int nIn;
	in >> nIn; // length of ordering
	if (nIn != m_n && nIn != m_nOrg) {
		std::cerr << "Problem reading ordering, number of variables doesn't match"
				<< std::endl;
		in.close();
		exit(EXIT_FAILURE);
	}

	// read into buffer first
	std::list<int> buffer;
	int x = UNKNOWN;
	while (nIn-- && !in.eof()) {
		in >> x;
		buffer.push_back(x);
	}

	bool fullOrdering = false; // default
	if (buffer.size() == (size_t) m_nOrg) {
		fullOrdering = true;
	}

	int n = 0;
	vector<bool> check(m_n, false);
	for (std::list<int>::iterator it = buffer.begin(); it != buffer.end(); ++it) {
		if (!fullOrdering) {

			if (*it < 0 || *it >= m_n) {
				std::cerr << "Problem reading ordering, variable index " << *it
						<< " out of range" << std::endl;
				in.close();
				exit(EXIT_FAILURE);
			}

			if (check[*it]) {
				cerr << "Problem reading ordering, variable " << *it
						<< " appears more than once." << std::endl;
				in.close();
				exit(EXIT_FAILURE);
			} else
				check[*it] = true;

			elim.push_back(*it);
			++n;

		} else { // full order, needs filtering

			if (*it < 0 || *it >= m_nOrg) {
				std::cerr << "Problem reading ordering, variable index " << *it
						<< " out of range" << std::endl;
				in.close();
				exit(EXIT_FAILURE);
			}

			std::map<int, int>::const_iterator it2 = m_old2new.find(*it);
			if (it2 != m_old2new.end()) {
				x = it2->second;
				if (check[x]) {
					std::cerr << "Problem reading ordering, variable " << *it
							<< " appears more than once." << endl;
					in.close();
					exit(EXIT_FAILURE);
				} else
					check[x] = true;

				elim.push_back(x);
				++n;
			} else { /* evidence */
			}

		}

	}

	if (n != m_n) {
		std::cerr << "Problem reading ordering, number of variables doesn't match."
				<< std::endl;
		in.close();
		exit(EXIT_FAILURE);
	}

	in.close();
	return true;
}

bool Problem::parseUAI(const string& prob, const string& evid, const bool positive) {

	// Open the problem file
	std::ifstream in(prob.c_str());
	if (in.fail()) {
		std::cerr << "Error opening the problem file: " << prob << std::endl;
		in.close();
		return false;
	}

	size_t found = prob.find_last_of("/");
	m_name = (found != std::string::npos) ? prob.substr(found + 1) : prob;
	std::cout << "Reading problem " << m_name << " ... ";

	std::vector<int> arity;
	std::vector<std::vector<int> > scopes;
	std::string s;
	int x, y;
	int xs;
	size_t z;

	in >> s; // Problem type
	bool forcePosDist = positive;
	assert (s == "BAYES" || s == "MARKOV");

	in >> x; // No. of variables
	m_n = m_nOrg = x;
	m_domains.resize(m_n, UNKNOWN);
	m_k = -1;
	for (int i = 0; i < m_n; ++i) { // Domain sizes
		in >> x; // read into int first
		if (x > numeric_limits<int>::max()) {
			std::cerr << "Domain size " << x
					<< " out of range for internal representation.\n"
					<< "(Recompile with different type for variable values.)"
					<< std::endl;
			in.close();
			return false;
		}
		xs = (int) x;
		m_domains[i] = xs;
		m_k = max(m_k, xs);
	}

	in >> x; // No. of functions
	m_c = x;
	scopes.reserve(m_c);

	// Scope information for functions
	m_r = -1;
	for (int i = 0; i < m_c; ++i) {
		std::vector<int> scope;
		in >> x; // arity

		m_r = max(m_r, x);
		for (int j = 0; j < x; ++j) {
			in >> y; // the actual variables in the scope
			if (y >= m_n) {
				std::cerr << "Variable index " << y << " out of range." << endl;
				in.close();
				return false;
			}
			scope.push_back(y); // preserve order from file
		}
		scopes.push_back(scope);
	}

	// Read functions
	size_t numTuples = 0, numZeroTuples = 0;
	for (int i = 0; i < m_c; ++i) {
		in >> z; // No. of entries
		size_t tab_size = 1;

		for (vector<int>::iterator it = scopes[i].begin();
				it != scopes[i].end(); ++it) {
			tab_size *= m_domains[*it];
		}

		assert(tab_size==z);
		// product of domain sizes matches no. of entries

		// create set version of the scope (ordered)
		std::set<int> scopeSet(scopes[i].begin(), scopes[i].end());
		z = scopeSet.size();

		// compute reindexing map from specified scope to ordered, internal one
		std::map<int, int> mapping;
		int k = 0;
		for (std::vector<int>::const_iterator it = scopes[i].begin();
				it != scopes[i].end(); ++it)
			mapping[*it] = k++;
		std::vector<int> reidx(z);
		std::vector<int>::iterator itr = reidx.begin();
		for (std::set<int>::iterator it = scopeSet.begin(); itr != reidx.end();
				++it, ++itr) {
			*itr = mapping[*it];
		}

		// read the full table into an temp. array (to allow reordering)
		std::vector<double> temp(tab_size);
		for (size_t j = 0; j < tab_size; ++j) {
			in >> temp[j];
			numTuples++;
			if (temp[j] == 0.0) numZeroTuples++;
		}

		// force a positive distribution (for BAYES and MARKOV only)
		if (forcePosDist) {
			for (size_t j = 0; j < tab_size; ++j) {
				if (temp[j] == 0.0) {
					temp[j] += EPSILON;
				} else if (temp[j] == 1.0) {
					//temp[j] -= EPSILON;
				}
			}
		}

		// get the variable domain sizes
		std::vector<int> limit;
		limit.reserve(z);
		for (std::vector<int>::const_iterator it = scopes[i].begin();
				it != scopes[i].end(); ++it)
			limit.push_back(m_domains[*it]);
		std::vector<int> tuple(z, 0);

		// create the new table (with reordering)
		double* table = new double[tab_size];
		for (size_t j = 0; j < tab_size;) {
			size_t pos = 0, offset = 1;
			// j is the index in the temp. table
			for (k = z - 1; k >= 0; --k) { // k goes backwards through the ordered scope
				pos += tuple[reidx[k]] * offset;
				offset *= m_domains[scopes[i][reidx[k]]];
			}
			//table[pos] = (neglog) ? -ELEM_ENCODE(temp[j]) : temp[j];
			table[pos] = temp[j];
			increaseTuple(j, tuple, limit);
		}

		Function* f = new FunctionBayes(i, this, scopeSet, table, tab_size);
		f->setOriginal(); // mark it as original function
		m_functions.push_back(f);

	} // All function tables read
	in.close();
	std::cout << "done." << std::endl;
	std::cout << "Number of tuples:\t " << numTuples << std::endl;
	std::cout << "Number of zero tuples:\t " << numZeroTuples << std::endl;
	m_determinism = ((double)numZeroTuples/(double)numTuples) * 100;
	std::cout << "Determinism ratio:\t " << std::ceil(m_determinism) << " %" << std::endl;

	// Read evidence?
	if (evid.empty()) {
		m_e = 0;
		return true; // No evidence, return
	}

	std::cout << "Reading evidence ... ";

	in.open(evid.c_str());
	if (in.fail()) {
		std::cerr << "Error while opening the evidence file: " << evid << std::endl;
		in.close();
		return false;
	}

	in >> x;
	m_e = x; // Number of evidence

	for (int i = 0; i < m_e; ++i) {
		in >> x; // Variable index
		in >> y; // Variable value
		xs = (int) y;
		if (xs >= m_domains[x]) {
			std::cout << "Variable " << x << " has domain size "
					<< (int) m_domains[x] << ", evidence value " << y
					<< " out of range." << std::endl;
			in.close();
			return false;
		}
		m_evidence.insert(make_pair(x, xs));
	}

	std::cout << "done." << std::endl;

	in.close();
	return true;
}

bool Problem::parseERG(const string& prob, const string& evid, const bool positive) {

	// Open the problem file
	std::ifstream in(prob.c_str());
	if (in.fail()) {
		std::cerr << "Error opening the problem file: " << prob << std::endl;
		in.close();
		return false;
	}

	size_t found = prob.find_last_of("/");
	m_name = (found != std::string::npos) ? prob.substr(found + 1) : prob;
	std::cout << "Reading problem " << m_name << " ... ";

	std::vector<int> arity;
	std::vector<std::vector<int> > scopes;
	std::string s;
	int x, y;
	int xs;
	size_t z;

	bool forcePosDist = positive;

	in >> x; // No. of variables
	m_n = m_nOrg = x;
	m_domains.resize(m_n, UNKNOWN);
	m_k = -1;
	for (int i = 0; i < m_n; ++i) { // Domain sizes
		in >> x; // read into int first
		if (x > numeric_limits<int>::max()) {
			std::cerr << "Domain size " << x
					<< " out of range for internal representation.\n"
					<< "(Recompile with different type for variable values.)"
					<< std::endl;
			in.close();
			return false;
		}
		xs = (int) x;
		m_domains[i] = xs;
		m_k = max(m_k, xs);
	}

	m_c = m_n; // No. of functions (BAYES)
	scopes.reserve(m_c);

	// Scope information for functions
	m_r = -1;
	for (int i = 0; i < m_n; ++i) {
		std::vector<int> scope;
		in >> x; // number of parents

		m_r = max(m_r, x);
		for (int j = 0; j < x; ++j) {
			in >> y; // the actual variables in the scope
			if (y >= m_n) {
				std::cerr << "Variable index " << y << " out of range." << endl;
				in.close();
				return false;
			}
			scope.push_back(y); // preserve order from file
		}
		scope.push_back(i); // the child variable
		scopes.push_back(scope);
	}

	// Read functions
	size_t numTuples = 0, numZeroTuples = 0;
	for (int i = 0; i < m_c; ++i) {
		in >> z; // No. of entries
		size_t tab_size = 1;

		for (vector<int>::iterator it = scopes[i].begin();
				it != scopes[i].end(); ++it) {
			tab_size *= m_domains[*it];
		}

		assert(tab_size==z);
		// product of domain sizes matches no. of entries

		// create set version of the scope (ordered)
		std::set<int> scopeSet(scopes[i].begin(), scopes[i].end());
		z = scopeSet.size();

		// compute reindexing map from specified scope to ordered, internal one
		std::map<int, int> mapping;
		int k = 0;
		for (std::vector<int>::const_iterator it = scopes[i].begin();
				it != scopes[i].end(); ++it)
			mapping[*it] = k++;
		std::vector<int> reidx(z);
		std::vector<int>::iterator itr = reidx.begin();
		for (std::set<int>::iterator it = scopeSet.begin(); itr != reidx.end();
				++it, ++itr) {
			*itr = mapping[*it];
		}

		// read the full table into an temp. array (to allow reordering)
		std::vector<double> temp(tab_size);
		for (size_t j = 0; j < tab_size; ++j) {
			in >> temp[j];
			numTuples++;
			if (temp[j] == 0.0) numZeroTuples++;
		}

		// force a positive distribution (for BAYES and MARKOV only)
		if (forcePosDist) {
			for (size_t j = 0; j < tab_size; ++j) {
				if (temp[j] == 0.0) {
					temp[j] += EPSILON;
				} else if (temp[j] == 1.0) {
					temp[j] -= EPSILON;
				}
			}
		}

		// get the variable domain sizes
		std::vector<int> limit;
		limit.reserve(z);
		for (std::vector<int>::const_iterator it = scopes[i].begin();
				it != scopes[i].end(); ++it)
			limit.push_back(m_domains[*it]);
		std::vector<int> tuple(z, 0);

		// create the new table (with reordering)
		double* table = new double[tab_size];
		for (size_t j = 0; j < tab_size;) {
			size_t pos = 0, offset = 1;
			// j is the index in the temp. table
			for (k = z - 1; k >= 0; --k) { // k goes backwards through the ordered scope
				pos += tuple[reidx[k]] * offset;
				offset *= m_domains[scopes[i][reidx[k]]];
			}
			//table[pos] = (neglog) ? -ELEM_ENCODE(temp[j]) : temp[j];
			table[pos] = temp[j];
			increaseTuple(j, tuple, limit);
		}

		Function* f = new FunctionBayes(i, this, scopeSet, table, tab_size);
		f->setOriginal(); // mark it as original function
		m_functions.push_back(f);

	} // All function tables read
	in.close();
	std::cout << "done." << std::endl;
	std::cout << "Number of tuples:\t " << numTuples << std::endl;
	std::cout << "Number of zero tuples:\t " << numZeroTuples << std::endl;
	double det = ((double)numZeroTuples/(double)numTuples) * 100;
	std::cout << "Determinism ratio:\t " << std::ceil(det) << " %" << std::endl;

	// Read evidence?
	if (evid.empty()) {
		m_e = 0;
		return true; // No evidence, return
	}

	std::cout << "Reading evidence ... ";

	in.open(evid.c_str());
	if (in.fail()) {
		std::cerr << "Error while opening the evidence file: " << evid << std::endl;
		in.close();
		return false;
	}

	in >> x;
	m_e = x; // Number of evidence

	for (int i = 0; i < m_e; ++i) {
		in >> x; // Variable index
		in >> y; // Variable value
		xs = (int) y;
		if (xs >= m_domains[x]) {
			std::cout << "Variable " << x << " has domain size "
					<< (int) m_domains[x] << ", evidence value " << y
					<< " out of range." << std::endl;
			in.close();
			return false;
		}
		m_evidence.insert(make_pair(x, xs));
	}

	std::cout << "done." << std::endl;

	in.close();
	return true;
}

// Assumes that (pseudo) evidence was NOT removed from the problem
void Problem::writeUAI(const string& prob) const {
	assert(prob.size());

	std::ofstream out;
	out.open(prob.c_str(), ios::out | ios::trunc);

	if (!out) {
		cerr << "Error writing network to file " << prob << endl;
		exit(1);
	}

	out << "MARKOV" << endl; // TODO hard-coding is not optimal

	// variable info
	out << (m_n) << endl;
	for (vector<val_t>::const_iterator it = m_domains.begin();
			it != m_domains.end(); ++it)
		out << ' ' << ((int) *it);

	// function information
	out << endl << m_functions.size() << endl; // no. of functions
	for (vector<Function*>::const_iterator it = m_functions.begin();
			it != m_functions.end(); ++it) {
		const set<int>& scope = (*it)->getScope();
		out << scope.size() << '\t'; // scope size
		for (set<int>::const_iterator itS = scope.begin(); itS != scope.end();
				++itS)
			out << *itS << ' '; // variables in scope
		out << endl;
	}
	out << endl;

	// write the function tables
	for (vector<Function*>::const_iterator it = m_functions.begin();
			it != m_functions.end(); ++it) {
		double * T = (*it)->getTable();
		out << (*it)->getTableSize() << endl; // table size
		for (size_t i = 0; i < (*it)->getTableSize(); ++i)
			out << ' ' << T[i]; // table entries
		out << endl;
	}

	// done
	out << endl;
	out.close();

}

// this is called after removeEvidence
bool Problem::parseMapVars(const string& file, std::vector<int>& mapVars) {

	// Safety checks
	assert( m_n != UNKNOWN );

	// Open the file
	std::ifstream in(file.c_str());
	if (in.fail()) { // file not existent yet
		std::cerr << "Error while opening the MAP variables file." << std::endl;
		exit(EXIT_FAILURE);
	}

	// ignore first line if there's a pound '#' sign (comment)
	if (in.peek() == '#') {
		in.ignore(8192, '\n');
	}

	int nIn;
	in >> nIn; // length of the set of variables

	// read into buffer first
	m_query.clear();
	std::list<int> buffer;
	int x = UNKNOWN;
	while (nIn-- && !in.eof()) {
		in >> x;
		buffer.push_back(x);
		m_query.insert(x); // original var indexes
	}

	int n = 0;
	vector<bool> check(m_n, false);
	for (std::list<int>::iterator it = buffer.begin(); it != buffer.end(); ++it) {

		if (*it < 0 || *it >= m_nOrg) {
			std::cerr << "Problem reading ordering, variable index " << *it
					<< " out of range (0, " << m_nOrg << "]" << std::endl;
			in.close();
			exit(EXIT_FAILURE);
		}

		std::map<int, int>::const_iterator it2 = m_old2new.find(*it);
		if (it2 != m_old2new.end()) {
			x = it2->second;
			if (check[x]) {
				std::cerr << "Problem reading map vars, variable " << *it
						<< " appears more than once." << endl;
				in.close();
				exit(EXIT_FAILURE);
			} else
				check[x] = true;

			mapVars.push_back(x);
			++n;
		}
	}

	in.close();
	return true;
}

// this can be called before removeEvidence
bool Problem::loadMapVars(const string& file, std::vector<int>& mapVars) {

	// Safety checks
	assert( m_n != UNKNOWN );

	// Open the file
	std::ifstream in(file.c_str());
	if (in.fail()) { // file not existent yet
		std::cerr << "Error while opening the MAP variables file." << std::endl;
		exit(EXIT_FAILURE);
	}

	// ignore first line if there's a pound '#' sign (comment)
	if (in.peek() == '#') {
		in.ignore(8192, '\n');
	}

	int nIn;
	in >> nIn; // length of the set of variables

	// read into buffer first
	m_query.clear();
	std::list<int> buffer;
	int x = UNKNOWN;
	while (nIn-- && !in.eof()) {
		in >> x;
		buffer.push_back(x);
		m_query.insert(x); // original var indexes
	}

	for (std::list<int>::iterator it = buffer.begin(); it != buffer.end(); ++it) {
		x = (*it);
		mapVars.push_back(x);
	}

	in.close();
	return true;
}

// this can be called before removeEvidence
bool Problem::loadOrdering(const string& file, std::vector<int>& order) const {

	// Safety checks
	assert( m_n != UNKNOWN );

	// Open the file
	std::ifstream in(file.c_str());
	if (in.fail()) { // file not existent yet
		std::cerr << "Error while opening the ordering file." << std::endl;
		exit(EXIT_FAILURE);
	}

	// ignore first line if there's a pound '#' sign (comment)
	if (in.peek() == '#') {
		in.ignore(8192, '\n');
	}

	int nIn;
	in >> nIn; // length of the set of variables

	// read into buffer first
	std::list<int> buffer;
	int x = UNKNOWN;
	while (nIn-- && !in.eof()) {
		in >> x;
		buffer.push_back(x);
	}

	for (std::list<int>::iterator it = buffer.begin(); it != buffer.end(); ++it) {
		x = (*it);
		order.push_back(x);
	}

	in.close();
	return true;
}


void Problem::updateVartypes(const std::vector<int>& mapVars) {
	m_vartypes.resize(m_n);
	for (size_t i = 0; i < m_vartypes.size(); ++i) {
		m_vartypes[i] = VAR_SUM;
	}
	for (size_t i = 0; i < mapVars.size(); ++i) {
		m_vartypes[mapVars[i]] = VAR_MAX;
	}
}

void Problem::setVartype(int var, var_t t) {
	assert(var >= 0 && var < (int)m_vartypes.size());
	m_vartypes[var] = t;
}

bool Problem::isEliminated(int i) const {
	std::map<int, int>::const_iterator itRen = m_old2new.find(i);
	return (itRen == m_old2new.end());
}

bool Problem::isMap(int i) const {
	assert(i >= 0 && i < m_n);
	return (m_vartypes[i] == VAR_MAX ? true : false);
}

bool Problem::isSum(int i) const {
	assert(i >= 0 && i < m_n);
	return (m_vartypes[i] == VAR_SUM ? true : false);
}

void Problem::replaceFunctions(const vector<Function*>& newFunctions) {
	// delete current functions
	for (vector<Function*>::iterator it = m_functions.begin();
			it != m_functions.end(); ++it) {
		if (*it)
			delete (*it);
	}
	// store new functions
	m_functions = newFunctions;
	m_c = m_functions.size();
	// update function scopes???
}

// create a copy of the current problem instance
Problem* Problem::clone() {

	Problem* cl = new Problem();

	// copy the functions
	for (vector<Function*>::iterator it = m_functions.begin();
			it != m_functions.end(); ++it) {
		Function* fn = (*it);
		cl->addFunction(fn->clone());
	}

	cl->m_hasDummy = m_hasDummy;
	cl->m_hasDummyFuns = m_hasDummyFuns;
	cl->m_n = m_n;
	cl->m_nOrg = m_nOrg;
	cl->m_k = m_k;
	cl->m_e = m_e;
	cl->m_c = m_c;
	cl->m_r = m_r;
	cl->m_globalConstant = m_globalConstant;
	cl->m_name = m_name;
	cl->m_domains = m_domains;
	cl->m_vartypes = m_vartypes;
	cl->m_evidence = m_evidence;
	cl->m_old2new = m_old2new;
	cl->m_dummyVar = m_dummyVar;
	cl->m_determinism = m_determinism;
	cl->m_query = m_query;

	return cl;
}

// create a new problem conditioned on a partial assignment
Problem* Problem::conditioned(const std::vector<int>& assignment) {

	Problem* cond = this->clone();
	cond->setEvidence(assignment);
	cond->removeEvidence();
	return cond;
}

void Problem::setEvidence(const std::vector<int>& assignment) {
	assert( (int)assignment.size() == m_n );
	m_evidence.clear();
	for (size_t var = 0; var < assignment.size(); ++var) {
		if (assignment[var] != UNKNOWN) {
			m_evidence[var] = assignment[var];
		}
	}
}

