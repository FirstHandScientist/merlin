/*
 * FactorSet.h
 *
 *  Created on: 14 Aug 2016
 *      Author: radu
 */

#ifndef LIB_MEX_FACTORSET_H_
#define LIB_MEX_FACTORSET_H_

#include "Factor.h"

namespace mex {

class FactorSet : public virtual mxObject {
public:
	/// Useful typedefs
	typedef double value;
	typedef mex::midx<uint32_t>    				findex;		// factor index
	typedef mex::midx<uint32_t>    				vindex;		// variable index
	typedef mex::set<mex::midx<uint32_t> >  	flist; 		// collection of factor indices

	// Constructors ////////////////////////////////////////////////////////////

	// copy ctor
	FactorSet(FactorSet const& fs) : _factors(fs._factors) {}

	// default constructor
	FactorSet() : _factors() {};

	// constant factor over given vars
	FactorSet(VarSet const& vs, value s = 1.0) {
		Factor f(vs, s);
		_factors.push_back(f);
	}

	FactorSet(VarSet const& vs, value* T) {
		Factor f(vs, T);
		_factors.push_back(f);
	}

	FactorSet(const vector<Factor>& fs) : _factors(fs) {};

    // destructor
	~FactorSet() {};

	void add(const Factor& f) {
		_factors.push_back(f);
	}

	// number of factors
	size_t numel() const {
		return _factors.size();
	}

	size_t size() const {
		return _factors.size();
	}

	bool empty() const {
		return _factors.empty();
	}

	Factor& operator[](size_t v) {
		return _factors[v];
	}

	Factor operator[](size_t v) const {
		return _factors[v];
	}

	// get the factors
	vector<Factor>& factors() {
		return _factors;
	}

	// return the scope of the current factor set
	VarSet variables() {
		if (empty()) {
			return VarSet(); // empty set
		} else {
			return _factors[0].variables();
		}
	}

	// output operator
	friend std::ostream& operator<<(std::ostream& out, const mex::FactorSet& FS) {
		if (FS.empty()) {
			out << "FactorSet over []\n";
		} else {
			out << "FactorSet over " << FS[0].variables() << "\n";
			for (size_t j = 0; j < FS.size(); j++)
				out << " " << FS[j] << std::endl;
		}
		return out;
	}

	// Compute an upper bounding factor for the current factor set
	Factor upbo() {
		if (empty()) {
			return Factor(); // constant 1 factor
		} else {
			VarSet vs = this->variables();
			Factor ub(vs);
			for (size_t i = 0; i < ub.numel(); ++i) {
				value m = -infty();
				for (size_t j = 0; j < _factors.size(); ++j) {
					m = std::max(m, _factors[j][i]);
				}
				ub[i] = m;
			}

			return ub;
		}
	}

	// prune the current factor set by removing the pareto dominated ones
	FactorSet& pareto() {

		vector<Factor> V;

		// Iterate through the factor set
		for (size_t i = 0; i < _factors.size(); ++i) {
			const Factor& u = _factors[i];

			// Check if the factor 'u' is Pareto dominated by V
			bool found = false;
			for (size_t j = 0; j < V.size(); ++j) {
				const Factor& v = V[j];

				// Check for Pareto dominance
				if (v.dominates(u)) {
					found = true; // ignore the tuple 'u'
					break;
				}
			}

			if (!found) {
				// Remove from V all factors that are Pareto dominated by u
				vector<Factor>::iterator it = V.begin();
				for (; it != V.end(); ) {
					const Factor& v = (*it);
					if (u.dominates(v)) { // Check Pareto dominance
						V.erase(it);
					} else {
						++it;
					}
				}

				V.push_back(u);
			}
		}

		// Now V contains all the non-dominated factors
		_factors.swap(V);

		return (*this);
	}

	// Measure how much factor f is a "representative" of another factor g
	double repr(const Factor& f, const Factor& g) const {
		Factor tmp = f/g;
		return tmp.max();
	}


	// Given the current factor set, generate k clusters such that each cluster
	// minmizes the error measure from the paper. Basically, select k factors,
	// which are called the "representative" factors, then assign the remaining
	// ones to the cluster that minimizes the repr() measure wrt the "representative"
	// of that cluster. Should overall minimize the error.
	void prune_mu(size_t k, FactorSet& result, double& error) {

		this->pareto(); // remove pareto dominated factors from this set

		if (_factors.size() <= k || k == 0) {
			//std::cout << "*No pruning mu\n";
			result = (*this); // nothing to do
			error = 0;
			return;
		}

		assert(_factors.size() > k);

		std::list<Factor> remaining;
		vector<vector<Factor> > clusters;
		for (size_t i = 0; i < k; ++i) {
			vector<Factor> cl;
			cl.push_back(_factors[i]);
			clusters.push_back(cl);
		}

		for (size_t j = k; j < _factors.size(); ++j) {
			remaining.push_back(_factors[j]);
		}

		// Assign the remaining factors to the appropriate cluster
		while (!remaining.empty()) {
			Factor f = remaining.front();
			remaining.pop_front();
			int cl = -1;
			double err = infty();
			for (size_t i = 0; i < k; ++i) {
				const Factor& g = clusters[i][0];
				double e = repr(g, f);
				if (e < err) {
					cl = i;
					err = e;
				}
			}

			assert(cl >= 0);
			clusters[cl].push_back(f);
		}

		// Create the representative factor set
		error = 0; // epsilon(Vi)
		FactorSet Vi;
		for (size_t i = 0; i < k; ++i) {
			Vi.add(clusters[i][0]);
			double err = 0;
			for (size_t j = 1; j < clusters[i].size(); ++j) {
				const Factor& mu = clusters[i][j];
				err = std::max(err, repr(mu, clusters[i][0]));
			}

			error = std::max(error, err);
		}

		result = Vi;
	}

	// Given the current factor set, generate k clusters such that each cluster
	// minmizes the error measure from the paper. Basically, select k factors,
	// which are called the "representative" factors, then assign the remaining
	// ones to the cluster that minimizes the repr() measure wrt the "representative"
	// of that cluster. Should overall minimize the error.
	void prune_ub(size_t k, FactorSet& upbo) {

		//this->pareto();

		if (_factors.size() <= k || k == 0) {
			//std::cout << "*No pruning ub\n";
			upbo = (*this); // nothing to do
			return;
		}

		assert(_factors.size() > k);

		std::list<Factor> remaining;
		vector<vector<Factor> > clusters;
		for (size_t i = 0; i < k; ++i) {
			vector<Factor> cl;
			cl.push_back(_factors[i]);
			clusters.push_back(cl);
		}

		for (size_t j = k; j < _factors.size(); ++j) {
			remaining.push_back(_factors[j]);
		}

		// Assign the remaining factors to the appropriate cluster
		while (!remaining.empty()) {
			Factor f = remaining.front();
			remaining.pop_front();
			int cl = -1;
			double err = infty();
			for (size_t i = 0; i < k; ++i) {
				const Factor& g = clusters[i][0];
				double e = repr(g, f);
				if (e < err) {
					cl = i;
					err = e;
				}
			}

			assert(cl >= 0);
			clusters[cl].push_back(f);
		}

		// Create the representative factor set
		FactorSet UB;
		for (size_t i = 0; i < k; ++i) {
			UB.add(FactorSet(clusters[i]).upbo());
		}

		upbo = UB;
	}

	void sum(Var VX) {
		for (size_t i = 0; i < _factors.size(); ++i) {
			_factors[i] = _factors[i].sum(VX);
		}
	}

	value max() {
		value v = -infty();
		for (size_t i = 0; i < _factors.size(); ++i) {
			v = std::max(v, _factors[i].max());
		}
		return v;
	}

protected:

	vector<Factor> _factors;		// factors in the set (defined over the same scope)
};


}; // end namespace



#endif /* LIB_MEX_FACTORSET_H_ */
