/*
 * mini_bucket_elim.cpp
 *
 *  Created on: Apr 2, 2013
 *      Author: radu
 */

#include "mini_bucket_heur.h"

//#define DEBUG

#ifdef DEBUG
/* ostream operator for debugging */
std::ostream& operator <<(std::ostream& os, const std::list<Function*>& l) {
	std::list<Function*>::const_iterator it = l.begin();
	os << '[';
	while (it!=l.end()) {
		os << (**it);
		if (++it != l.end()) os << ',';
	}
	os << ']';
	return os;
}
#endif


/* computes the augmented part of the heuristic estimate */
double MiniBucketHeur::getHeur(int var,
		const std::vector<int>& assignment) const {

	// safety checks
	assert( var >= 0 && var < m_problem->getN());

	// variable 'var' is assumed to be already assigned (in 'assignment')
	double h = ELEM_ONE;

	// go over augmented and intermediate lists and combine all values
	std::list<Function*>::const_iterator itF = m_augmented[var].begin();
	for (; itF != m_augmented[var].end(); ++itF) {
		h *= (*itF)->getValue(assignment);
	}

	itF = m_intermediate[var].begin();
	for (; itF != m_intermediate[var].end(); ++itF) {
		h *= (*itF)->getValue(assignment);
	}

	return h;
}

// computes the heuristic for variable var given a (partial) assignment
//double MiniBucketHeur::getHeur(int var, std::vector<double>& bounds,
//		const std::vector<int>& assignment) const {
//
//	// safety checks
//	assert( var >= 0 && var < m_problem->getN());
//
//	// variable 'var' is unassigned
//	assert( assignment[var] == NONE );
//
//	double upbo = -INFINITY;
//	std::vector<int> temp = assignment;
//	for (val_t i = 0; i < m_problem->getDomainSize(var); ++i) {
//
//		temp[var] = i;
//		double h = ELEM_ONE;
//
//		// go over augmented and intermediate lists and combine all values
//		std::list<Function*>::const_iterator itF = m_augmented[var].begin();
//		for (; itF != m_augmented[var].end(); ++itF) {
//			h *= (*itF)->getValue(temp);
//		}
//
//		itF = m_intermediate[var].begin();
//		for (; itF != m_intermediate[var].end(); ++itF) {
//			h *= (*itF)->getValue(temp);
//		}
//
//		bounds.push_back(h);
//		upbo = std::max(upbo, h);
//	}
//
//	return upbo;
//}

void MiniBucketHeur::reset() {

  std::vector<std::list<Function*> > empty;
  m_augmented.swap(empty);

  std::vector<std::list<Function*> > empty2;
  m_intermediate.swap(empty2);

}


size_t MiniBucketHeur::build(const std::vector<int> * assignment,
		bool computeTables) {

	// includes the dummy variable as well (last in the ordering)

#ifdef DEBUG
	cout << "$ Building MBE(" << m_ibound << ")" << endl;
#endif

	this->reset();

	std::vector<int> elimOrder; // will hold dfs order
	findDfsOrder(elimOrder); // computes dfs ordering of relevant subtree

	m_augmented.resize(m_problem->getN());
	m_intermediate.resize(m_problem->getN());

	// keep track of total memory consumption
	size_t memSize = 0;

	// ITERATES OVER BUCKETS, FROM LEAVES TO ROOT
	for (std::vector<int>::reverse_iterator itV = elimOrder.rbegin();
			itV != elimOrder.rend(); ++itV) {

#ifdef DEBUG
		cout << "$ Bucket for variable " << *itV << endl;
#endif

		// collect relevant functions in funs
		std::vector<Function*> funs;
		const std::list<Function*>& fnlist =
				m_pseudotree->getNode(*itV)->getFunctions();
		funs.insert(funs.end(), fnlist.begin(), fnlist.end());
		funs.insert(funs.end(), m_augmented[*itV].begin(),
				m_augmented[*itV].end());

#ifdef DEBUG
		for (vector<Function*>::iterator itF=funs.begin(); itF!=funs.end(); ++itF)
		cout << ' ' << (**itF);
		cout << endl;
#endif

		// compute global upper bound for root (dummy) bucket
		if (*itV == elimOrder[0]) { // variable is dummy root variable
			if (computeTables && assignment) { // compute lower bound if assignment is given
				m_globalUB = ELEM_ONE;
				for (std::vector<Function*>::iterator itF = funs.begin();
						itF != funs.end(); ++itF) {
					m_globalUB *= (*itF)->getValue(*assignment); // all constant functions, tablesize==1
				}

				cout << "    Upper Bound = " << m_globalUB << std::endl;
				m_globalUB *= m_problem->getGlobalConstant();
				cout << "    MBE-ALL  = " << m_globalUB << std::endl;
			}
			continue; // skip the dummy variable's bucket
		}

		// sort functions by decreasing scope size
		std::sort(funs.begin(), funs.end(), scopeIsLarger);

		// partition functions into minibuckets
		std::vector<MiniBucket> minibuckets;
		for (std::vector<Function*>::iterator itF = funs.begin();
				itF != funs.end(); ++itF) {
			bool placed = false;
			for (vector<MiniBucket>::iterator itB = minibuckets.begin();
					!placed && itB != minibuckets.end(); ++itB) {
				if (itB->allowsFunction(*itF)) { // checks if function fits into bucket
					itB->addFunction(*itF);
					placed = true;
				}
			}
			if (!placed) { // no fit, need to create new bucket
				MiniBucket mb(*itV, m_ibound, m_problem);
				mb.addFunction(*itF);
				minibuckets.push_back(mb);
			}
		}

#ifdef DEBUG
		cout << " There are " << minibuckets.size() << " mini-buckets" << endl;
#endif

		// minibuckets for current bucket are now ready, process each
		// and place resulting function
		bool first = true; // flag indicating the first mini-bucket in the list
		for (std::vector<MiniBucket>::iterator itB = minibuckets.begin();
				itB != minibuckets.end(); ++itB) {

			Function* newf = itB->eliminate(computeTables, first); // process the minibucket
			const set<int>& newscope = newf->getScope();
			memSize += newf->getTableSize();
			// go up in tree to find target bucket
			PseudotreeNode* n = m_pseudotree->getNode(*itV)->getParent();
			while (newscope.find(n->getVar()) == newscope.end()
					&& n != m_pseudotree->getRoot()) {
				m_intermediate[n->getVar()].push_back(newf);
				n = n->getParent();
			}
			// matching bucket found OR root of pseudo tree reached
			m_augmented[n->getVar()].push_back(newf);
			if (first) first = false;

#ifdef DEBUG
			cout << "  Generated function: " << *newf << std::endl;
#endif
		}
		// all minibuckets processed and resulting functions placed
	}

#ifdef DEBUG
	// output augmented and intermediate buckets
	if (computeTables) {
		for (int i=0; i<m_problem->getN(); ++i) {
			cout << "$ AUG" << i << ": " << m_augmented[i] << " + " << m_intermediate[i] << endl;
		}
	}
#endif

	// clean up for estimation mode
	if (!computeTables) {
		for (std::vector<std::list<Function*> >::iterator itA =
				m_augmented.begin(); itA != m_augmented.end(); ++itA)
			for (std::list<Function*>::iterator itB = itA->begin();
					itB != itA->end(); ++itB)
				delete *itB;
		m_augmented.clear();
		m_augmented.clear();
	}

	return memSize;
}


/* finds a dfs order of the pseudotree (or the locally restricted subtree)
 * and writes it into the argument vector */
void MiniBucketHeur::findDfsOrder(vector<int>& order) const {
	order.clear();
	std::stack<PseudotreeNode*> dfs;
	dfs.push(m_pseudotree->getRoot());
	PseudotreeNode* n = NULL;
	while (!dfs.empty()) {
		n = dfs.top();
		dfs.pop();
		order.push_back(n->getVar());
		for (std::vector<PseudotreeNode*>::const_iterator it =
				n->getChildren().begin(); it != n->getChildren().end(); ++it) {
			dfs.push(*it);
		}
	}
}


size_t MiniBucketHeur::limitSize(size_t memlimit, const std::vector<int> * assignment) {

  // convert to bits
  memlimit *= 1024 *1024 / sizeof(double);

  int ibound = m_options->ibound;

  cout << "Adjusting mini bucket i-bound..." << endl;
  this->setIbound(ibound);
  size_t mem = this->build(assignment, false);
  cout << " i=" << ibound << " -> " << ((mem / (1024*1024.0)) * sizeof(double) )
       << " MBytes" << endl;

  while (mem > memlimit && ibound > 1) {
    this->setIbound(--ibound);
    mem = this->build(assignment, false);
    cout << " i=" << ibound << " -> " << ((mem / (1024*1024.0)) * sizeof(double) )
         << " MBytes" << endl;
  }

  m_options->ibound = ibound;
  return mem;
}


size_t MiniBucketHeur::getSize() const {
  size_t S = 0;
  for (std::vector<std::list<Function*> >::const_iterator it=m_augmented.begin(); it!= m_augmented.end(); ++it) {
    for (std::list<Function*>::const_iterator itF=it->begin(); itF!=it->end(); ++itF)
      S += (*itF)->getTableSize();
  }
  return S;
}


/*
 * mini bucket file format (all data in binary):
 * - size_t: no. of variables
 * - int: i-bound
 * - double: global upper bound
 * for every variable:
 *   - size_t: number of functions in bucket structure
 *   for every such function:
 *     - int: function ID
 *     - size_t: scope size
 *     for every scope variable:
 *       - int: variable index
 *     - size_t: table size
 *     for every table entry:
 *       - double: CPT entry
 * for every variable:
 *   - size_t: number of intermediate function pointers
 *   for every function pointer:
 *     - size_t: function index (implicit order from above)
 */

bool MiniBucketHeur::writeToFile(std::string fn) const {

	std::ofstream out(fn.c_str());
	if (!out) {
		cerr << "Error writing mini buckets to file " << fn << endl;
		return false;
	}

	// used later
	int x = NONE;
	size_t y = NONE;

	// number of variables
	size_t sz = m_augmented.size();
	out.write((char*) &(sz), sizeof(sz));

	// i-bound
	out.write((char*) &(m_ibound), sizeof(m_ibound));

	// global UB
	out.write((char*) &(m_globalUB), sizeof(m_globalUB));

	std::map<const Function*, size_t> funcMap;

	// over m_augmented
	for (size_t i = 0; i < sz; ++i) {
		size_t sz2 = m_augmented[i].size();
		out.write((char*) &(sz2), sizeof(sz2));

		std::list<Function*>::const_iterator itF = m_augmented[i].begin();
		for (size_t j = 0; j < sz2; ++j, ++itF) {
			const Function* f = *itF;
			funcMap.insert(std::make_pair(f, funcMap.size()));

			// function ID
			int id = f->getId();
			out.write((char*) &(id), sizeof(id));

			// scope
			size_t sz3 = f->getScope().size();
			out.write((char*) &(sz3), sizeof(sz3));
			// scope size
			for (set<int>::const_iterator it = f->getScope().begin();
					it != f->getScope().end(); ++it) {
				x = *it;
				out.write((char*) &(x), sizeof(x));
				// vars from scope
			}

			// table size
			sz3 = f->getTableSize();
			out.write((char*) &(sz3), sizeof(sz3));

			// table
			out.write((char*) (f->getTable()), sizeof(double) * sz3);

		}

	}

	// over m_intermediate
	for (size_t i = 0; i < sz; ++i) {
		size_t sz2 = m_intermediate[i].size();
		out.write((char*) &(sz2), sizeof(sz2));

		std::list<Function*>::const_iterator itF = m_intermediate[i].begin();
		for (size_t j = 0; j < sz2; ++j, ++itF) {
			y = funcMap.find(*itF)->second;
			out.write((char*) &(y), sizeof(y));
		}

	}

	out.close();

	return true;

}

bool MiniBucketHeur::readFromFile(string fn) {

	std::ifstream inTemp(fn.c_str());
	inTemp.close();
	if (inTemp.fail()) { // file not existent yet
		return false;
	}

	std::ifstream in(fn.c_str());

	this->reset();

	// used later
	int x = NONE;
	size_t y = NONE;
	std::vector<Function*> allFuncs;

	// no. of variables
	size_t sz;
	in.read((char*) &(sz), sizeof(sz));

	if (sz != (size_t) m_problem->getN()) {
		cerr << "Number of variables in mini bucket file doesn't match" << endl;
		return false;
	}

	m_augmented.resize(sz);
	m_intermediate.resize(sz);

	// i-bound
	int ibound;
	in.read((char*) &(ibound), sizeof(ibound));
	m_ibound = ibound;

	// global UB
	double ub;
	in.read((char*) &(ub), sizeof(ub));
	m_globalUB = ub;

	// over variables for m_augmented
	for (size_t i = 0; i < sz; ++i) {
		size_t sz2;
		in.read((char*) &(sz2), sizeof(sz2));

		// over functions
		for (size_t j = 0; j < sz2; ++j) {
			int id;
			in.read((char*) &(id), sizeof(id));

			// scope
			size_t sz3;
			in.read((char*) &(sz3), sizeof(sz3));
			std::set<int> scope;
			for (size_t k = 0; k < sz3; ++k) {
				in.read((char*) &(x), sizeof(x));
				scope.insert(x);
			}

			// table size and table
			in.read((char*) &(sz3), sizeof(sz3));
			double* T = new double[sz3];
			in.read((char*) (T), sizeof(double) * sz3);

			// create function and store it
			Function* f = new FunctionBayes(id, m_problem, scope, T, sz3);
			m_augmented[i].push_back(f);
			allFuncs.push_back(f);
		}
	}

	for (size_t i = 0; i < sz; ++i) {
		// no. of function pointers
		size_t sz2;
		in.read((char*) &(sz2), sizeof(sz2));

		for (size_t j = 0; j < sz2; ++j) {
			// function index
			in.read((char*) &(y), sizeof(y));
			m_intermediate[i].push_back(allFuncs.at(y));
		}
	}

	in.close();
	cout << "Read mini bucket with i-bound " << ibound << " from file " << fn
			<< std::endl;
	return true;
}

