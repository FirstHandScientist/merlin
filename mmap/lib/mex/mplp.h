#ifndef __MEX_MPLP_H
#define __MEX_MPLP_H

#include <assert.h>
#include <stdexcept>
#include <stdlib.h>
#include <stdint.h>
#include <cstdlib>
#include <cstring>

#include "factorgraph.h"
#include "alg.h"

namespace mex {

// Factor graph algorithm specialization to MPLP-like optimization
// 
// derived factorGraph contains original factors
// beliefs form a reparameterization of the same distribution

class mplp: public factorGraph, public gmAlg, virtual public mxObject {
public:
	typedef factorGraph::findex findex;				// factor index
	typedef factorGraph::vindex vindex;				// variable index
	typedef factorGraph::flist flist; 			// collection of factor indices

public:
	mplp() : factorGraph() {
		setProperties();
	}
	mplp(const factorGraph& fg) : factorGraph(fg) {
		setProperties();
	}
	mplp(vector<Factor> fs) : factorGraph(fs) {
		setProperties();
	}
	template<class InputIterator>
	mplp(InputIterator f, InputIterator l) : factorGraph(f, l) {
		setProperties();
	}

	virtual mplp* clone() const {
		mplp* fg = new mplp(*this);
		return fg;
	}

	double lb() const {
		return 0;
	}	// !!!
	double ub() const {
		return _UB;
	}
	vector<index> best() const {
		return vector<index>();
	} // !!!

	// Not a summation algorithm
	double logZ() const {
		throw std::runtime_error("Not available");
	}
	double logZub() const {
		throw std::runtime_error("Not available");
	}
	double logZlb() const {
		throw std::runtime_error("Not available");
	}

	Factor& belief(findex f) {
		return _beliefs[f];
	}	//!!! const
	const Factor& belief(size_t f) const {
		return _beliefs[f];
	}
	const Factor& belief(Var v) const {
		return belief(localFactor(v));
	}
	const Factor& belief(VarSet vs) const {
		throw std::runtime_error("Not implemented");
	}
	const vector<Factor>& beliefs() const {
		return _beliefs;
	}
	Factor& belief(Var v) {
		return belief(localFactor(v));
	}

	MEX_ENUM( Property , StopTime, StopIter)
	;

	virtual void setProperties(std::string opt = std::string()) {
		if (opt.length() == 0) {
			setProperties("StopTime=-1,StopIter=-1");
			return;
		}
		std::vector<std::string> strs = mex::split(opt, ',');
		for (size_t i = 0; i < strs.size(); ++i) {
			std::vector<std::string> asgn = mex::split(strs[i], '=');
			switch (Property(asgn[0].c_str())) {
			case Property::StopTime:
				_stopTime = strtod(asgn[1].c_str(), NULL);
				break;
			case Property::StopIter:
				_stopIter = strtod(asgn[1].c_str(), NULL);
				break;
			default:
				break;
			}
		}
	}

	void init(const VarSet& vs) {
		init();
	}	// !!! inefficient

	void init() {
		_UB = 0.0; // not really an upper bound (approximate)
		_beliefs.resize(this->nvar());
	}

	void run() {
		double startTime = timeSystem();

		std::cout << "Mixed BP preprocessing ... " << std::endl;
		std::cout << "  iter: " << _stopIter << ", time: " << _stopTime << std::endl;

		int iter = 0;
		double tol = 1e-6;
		double Err = infty();
		vector<double> dists;
		_factors2 = this->factors();

		while (Err > tol) {
			if (_stopIter > 0 && iter >= _stopIter)
				break;		// iterations check
			if (_stopTime > 0 && _stopTime <= (timeSystem() - startTime))
				break;       // time-out check

			_dists.resize(nvar(), 0.0);
			for (size_t i = 0; i < nvar(); ++i) {
				Var v = var(i);
				updateVar(v);
			}

			Err = std::accumulate(_dists.begin(), _dists.end(), 0.0, std::plus<double>());

			// print some information for convergence
			if (iter%10 == 0) {
				std::cout << "  UB = " << _UB << "; err: " << Err << "; iter: " << iter << std::endl;
			}

			iter++;

			// find a configuration of the max variables from their local beliefs
			_solution.resize(nvar(), -1);
			for (VarSet::const_iterator it = _maxVars.begin();
					it != _maxVars.end(); ++it) {
				Var vv = (*it);
				size_t xtmp = _beliefs[vv].argmax();
				_solution[vv] = xtmp;
			}
		}
	}

	void setMaxVars(const std::list<int>& vars) {
		for (std::list<int>::const_iterator li = vars.begin(); li != vars.end(); ++li) {
			_maxVars.insert( var(*li) );
		}
	}

	std::vector<int>& getSolution() {
		return _solution;
	}

protected:
	// Contained objects
	vector<Factor> _beliefs;
	vector<Factor> _factors2;
	double _UB;
	double _stopTime;
	double _stopIter;
	VarSet _maxVars;
	vector<double> _dists;
	std::vector<int> _solution;

	// Update a single variable and all surrounding factors
	// Fixed point is all factors have equal max-marginals for all variables
	void updateVar(const Var& v) {

		// Damping factor / step size (helps convergence)
		double DAMP = 0.5;

		findex vf = localFactor(v);				// collect factors: var node +
		const mex::set<EdgeID>& nbrs = neighbors(vf); //   its neighbors
		bool isMax = ( _maxVars.contains(v) == true );

		Factor MM;

//		std::cout << "   updating var " << v << std::endl;
//		std::cout << "     # neighbors: " << nbrs.size() << std::endl;

		// compute each max-marginal and variable belief
		for (set<EdgeID>::const_iterator i = nbrs.begin(); i != nbrs.end(); ++i) {

			findex j = i->second;
			Factor balpha = _factors2[j];

//			std::cout << "    neighbor("<<j<<"): " << balpha << std::endl;
			VarSet vs = balpha.vars();
			for (VarSet::const_iterator it = vs.begin(); it != vs.end(); ++it) {
				Var vv = (*it);
				balpha *= _beliefs[vv];
			}

			Factor mm, bbeta = balpha.marginal(_maxVars);
			if (isMax) {
				mm = bbeta.maxmarginal(v) / _beliefs[v];
			} else {
				size_t xbeta = bbeta.argmax();
				Factor fcond = balpha.condition(bbeta.vars(), xbeta);
				mm = fcond.marginal(v) / _beliefs[v];
			}

			_factors2[j] = _factors2[j] / (mm^DAMP);
			MM = MM * mm^DAMP;
		}

		Factor tmp = _beliefs[v] * MM;
		if (isMax) {
			_UB = _UB + std::log(tmp.max());
			tmp /= tmp.max();
		} else {
			_UB = _UB + std::log(tmp.sum());
			tmp /= tmp.sum();
		}

		if (tmp.vars() == _beliefs[v].vars()) {
			_dists[v] = tmp.distance(_beliefs[v], Factor::Distance::L1);
		} else {
			_dists[v] = infty();
		}
		_beliefs[v] = tmp;
	}
};

//////////////////////////////////////////////////////////////////////////////////////////////
}
       // namespace mex
#endif  // re-include
