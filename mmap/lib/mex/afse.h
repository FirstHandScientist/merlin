/*
 * afse.h
 *
 *  Created on: Dec 19, 2013
 *      Author: radu
 */

#ifndef __MEX_MINIBUCKET_AFSE_H
#define __MEX_MINIBUCKET_AFSE_H

#include <assert.h>
#include <stdexcept>
#include <stdlib.h>
#include <stdint.h>
#include <cstdlib>
#include <cstring>

#include "factorgraph.h"
#include "alg.h"
#include "FactorSet.h"

namespace mex {

// Anytime Factor Set Elimination (AFSE) algorithm.
//
// This is an anytime Marginal MAP algorithm (by Maua and de Campos, 2012)

class afse: public graphModel, public gmAlg, virtual public mxObject {
public:
	typedef graphModel::findex findex;        // factor index
	typedef graphModel::vindex vindex;        // variable index
	typedef graphModel::flist flist;         // collection of factor indices

public:
	afse() : graphModel() {
		setProperties();
	}
	afse(const graphModel& gm) : graphModel(gm), _gmo(gm) {
		clearFactors();
		setProperties();
	}
	virtual afse* clone() const {
		afse* gm = new afse(*this);
		return gm;
	}

	graphModel _gmo;	// Original graphical model

	// Can be an optimization algorithm or a summation algorithm....
	double ub() const {
		return _logZ;
	}
	double lb() const {
		throw std::runtime_error("Not implemented");
	}
	vector<index> best() const {
		throw std::runtime_error("Not implemented");
	}

	double logZ() const {
		return _logZ;
	}
	double logZub() const {
		return _logZub;
	}
	double logZlb() const {
		return _logZlb;
	}

	// No beliefs defined currently
	const Factor& belief(size_t f) const {
		throw std::runtime_error("Not implemented");
	}
	const Factor& belief(Var v) const {
		throw std::runtime_error("Not implemented");
	}
	const Factor& belief(VarSet vs) const {
		throw std::runtime_error("Not implemented");
	}
	const vector<Factor>& beliefs() const {
		throw std::runtime_error("Not implemented");
	}

	const graphModel& gmOrig() const {
		return _gmo;
	}

	void build(const graphModel& gmo, size_t iBound, const VarOrder& elimOrder);
	virtual void run() {
	}

	MEX_ENUM( Property , Order,Distance )
	;

	bool _byScope;
	Factor::Distance distMethod;
	graphModel::OrderMethod ordMethod;
	size_t _iBound;
	double _logZ;
	double _logZub;
	double _logZlb;
	VarOrder _order;
	vector<flist> atElim;
	vector<double> atElimNorm;
	vector<vindex> _parents;
	vector<bool> _varTypes; // _varTypes[i] true if MAX, false if SUM

	///////////////// JGLP local structures ///////////////////////////////////
	vector<bool> _types;		// the type of each cluster (SUM or MAX)
	vector<flist> _clusters;	// clusters for each variable
	std::map<size_t, size_t> _cl2var; // cluster to var mapping
	size_t _numClusters;

	vector<flist> _originals;	// original factors (index) for each cluster
	vector<flist> _deltas;		// delta/indicator factors for each of the MAP variables
	vector<VarSet> _scopes;		// the scope (vars) for each cluster

	vector<vector<VarSet> > _separators;

	vector<flist> _in;			// incoming to each cluster
	vector<flist> _out; 		// outgoing from each cluster
	flist _roots;				// root cluster(s)

	vector<Factor> _forward; 	// forward messages (by edge)
	vector<FactorSet> _forward2; // forward factorset messages (by edge)
	vector<double> _errors;		// errors corresp. to factorset messages
	vector<size_t> _card;	// set cardinalities for clustering
	vector<FactorSet> _sigmas;	// upper bounds propagated

	vector<std::pair<findex, findex> > _schedule;
	vector<vector<size_t> > _edgeIndex;

	double _startTime;

	/////////////////////////////////////////////////////////////////
	// Setting properties (directly or through property string)
	/////////////////////////////////////////////////////////////////

	void setVarTypes(const vector<bool>& types) {
		_varTypes = types;
	}
	const vector<bool>& getVarTypes() const {
		return _varTypes;
	}

	void setOrder(const VarOrder& ord) {
		_order = ord;
	}
	void setOrder(graphModel::OrderMethod method) {
		_order.clear();
		ordMethod = method;
	}
	const VarOrder& getOrder() {
		return _order;
	}

	const vector<vindex>& getPseudotree() {
		return _parents;
	}
	void setPseudotree(const vector<vindex>& p) {
		_parents = p;
	}

	void setModel(const graphModel& gm) {
		_gmo = gm;
	}
	void setModel(const vector<Factor>& fs) {
		_gmo = graphModel(fs);
	}

	virtual void setProperties(std::string opt = std::string()) {
		if (opt.length() == 0) {
			setProperties("Order=MinFill");
			_byScope = true;
			return;
		}
		std::vector<std::string> strs = mex::split(opt, ',');
		for (size_t i = 0; i < strs.size(); ++i) {
			std::vector<std::string> asgn = mex::split(strs[i], '=');
			switch (Property(asgn[0].c_str())) {
			case Property::Order:
				_order.clear();
				_parents.clear();
				ordMethod = graphModel::OrderMethod(asgn[1].c_str());
				break;
			case Property::Distance:
				distMethod = Factor::Distance(asgn[1].c_str());
				_byScope = false;
				break;
			default:
				break;
			}
		}
	}

	// elimination (sum)
	Factor elim(const Factor& F, const VarSet& vs) {
		return F.sum(vs);
	}

	// marginals (sum)
	Factor marg(const Factor& F, const VarSet& vs) {
		return F.marginal(vs);
	}

	// Scoring function for bucket aggregation
	//   Unable to combine => -3; Scope-only => 1.0; otherwise a positive double score
	double score(const vector<Factor>& fin, const Var& VX, size_t i, size_t j) {
		double err;
		const Factor& F1 = fin[i], &F2 = fin[j];           // (useful shorthand)
		size_t iBound = std::max(std::max(_iBound, F1.nvar() - 1),
				F2.nvar() - 1);      // always OK to keep same size
		VarSet both = F1.vars() + F2.vars();
		if (both.nvar() > iBound)
			err = -3;  // too large => -3
		else
			err = 1.0 / (F1.nvar() + F2.nvar()); // greedy scope-based 2 (check if useful???)
		//else if (_byScope) err = 1;            // scope-based => constant score
		return err;
	}

	// Scoring function for bucket aggregation
	//   Unable to combine => -3; Scope-only => 1.0; otherwise a positive double score
	double score(const vector<VarSet>& fin, const Var& VX, size_t i, size_t j) {
		double err;
		const VarSet& F1 = fin[i], &F2 = fin[j];           // (useful shorthand)
		size_t iBound = std::max(std::max(_iBound, F1.nvar() - 1),
				F2.nvar() - 1);      // always OK to keep same size
		VarSet both = F1 + F2;
		if (both.nvar() > iBound)
			err = -3;  // too large => -3
		else
			err = 1.0 / (F1.nvar() + F2.nvar()); // greedy scope-based 2 (check if useful???)
		//else if (_byScope) err = 1;            // scope-based => constant score
		return err;
	}

	// helper class for pairs of sorted indices
	struct sPair: public std::pair<size_t, size_t> {
		sPair(size_t ii, size_t jj) {
			if (ii < jj) {
				first = jj;
				second = ii;
			} else {
				first = ii;
				second = jj;
			}
		}
	};

	void init(const VarSet& vs) {
		init();
	}            // !!! inefficient

	// create the mini-bucket tree (join graph): symbols only
	void init() {
		_startTime = timeSystem();
		if (_order.size() == 0) { // if we need to construct an elimination ordering
			double tic = timeSystem();
			_order = _gmo.order(ordMethod);
			_parents.clear(); // (new elim order => need new pseudotree) !!! should do together
			std::cout << "Order in " << timeSystem() - tic << " sec\n";
			std::cout << "Order: ";
			std::copy(_order.begin(), _order.end(),
					std::ostream_iterator<size_t>(std::cout, " "));
			std::cout << std::endl;
		}
		if (_parents.size() == 0) {     // if we need to construct a pseudo-tree
			double tic = timeSystem();
			_parents = _gmo.pseudoTree(_order);
			std::cout << "Pseudo in " << timeSystem() - tic << " sec\n";
		}

		// Get the factors scopes
		vector<VarSet> fin;
		for (vector<Factor>::const_iterator i = _gmo.factors().begin();
				i != _gmo.factors().end(); ++i) {
			fin.push_back((*i).vars());
		}

		// Mark factors depending on variable i
		vector<flist> vin;
		for (size_t i = 0; i < _gmo.nvar(); ++i) {
			vin.push_back(_gmo.withVariable(var(i)));
		}

		vector<flist> Orig(_gmo.nFactors()); 	// origination info: which original factors are
		for (size_t i = 0; i < Orig.size(); ++i)
			Orig[i] |= i;    					// included for the first time, and which newly
		vector<flist> New(_gmo.nFactors()); 	// created clusters feed into this cluster

		// First downward pass to initialize the bucket tree
		_clusters.resize(_order.size());
		for (VarOrder::const_iterator x = _order.begin(); x != _order.end(); ++x) {

			//std::cout << "Eliminating " << *x << (_varTypes[*x] ? "(MAP)\n" : "(SUM)\n");

			Var VX = var(*x);
			if (*x >= vin.size() || vin[*x].size() == 0)
				continue;  // check that we have some factors over this variable

			flist ids = vin[*x];  // list of factor IDs contained in this bucket

			// Merge all factor scopes in this bucket into a single scope
			vector<size_t> temp;
			typedef flist::const_iterator flistIt;
			size_t jj = *(ids.begin());
			for (flistIt i = ids.begin(); i != ids.end(); ++i) {
				size_t ii = *i;
				if (ii == jj) continue;
				fin[jj] |= fin[ii];     // combine into j
				erase(vin, ii, fin[ii]);
				fin[ii] = VarSet();  //   & remove i

				Orig[jj] |= Orig[ii];
				Orig[ii].clear(); // keep track of list of original factors in this cluster
				New[jj] |= New[ii];
				New[ii].clear(); //  list of new "message" clusters incoming to this cluster

				temp.push_back(ii);
			}

			for (size_t i = 0; i < temp.size(); ++i) {
				size_t ii = temp[i];
				ids /= ii;
			}

			// Sanity checks
			assert(ids.size() == 1);
//			std::cout << "  After merging: " << ids.size() << std::endl;
//			for (flist::const_iterator i = ids.begin(); i != ids.end(); ++i) {
//				std::cout << "  Factor id " << *i << std::endl;
//				std::cout << "  Scope: " << fin[*i] << std::endl;
//			}

			//// Eliminate bucket variable
			vector<findex> alphas;
			for (flistIt i = ids.begin(); i != ids.end(); ++i) {
				//
				// Create new cluster alpha over this set of variables; save function parameters also
				findex alpha = findex(-1);
				alpha = addFactor(Factor(fin[*i]));
				alphas.push_back(alpha);
				_clusters[*x] |= alpha;
				_cl2var[alpha] = *x;

				fin[*i] = fin[*i] - VX;

				// add inter clusters edges
				for (flistIt j = New[*i].begin(); j != New[*i].end(); ++j) {
					addEdge(*j, alpha);
					_schedule.push_back(std::make_pair(*j, alpha));
				}

				// update cluster types
				_types.push_back(_varTypes[*x]);

				// keep track of original factors
				_originals.push_back(flist());
				_originals[alpha] |= Orig[*i];

				// now incoming nodes to *i is just alpha
				Orig[*i].clear();
				New[*i].clear();
				New[*i] |= alpha;

				// recompute and update adjacency
				insert(vin, *i, fin[*i]);
			}
		}
		/// end for: variable elim order ///////////////////////////////////////

		// separators and cluster scopes
		size_t C = _factors.size();
		_numClusters = C;
		_separators.resize(C);
		for (size_t i = 0; i < C; ++i) _separators[i].resize(C);
		_scopes.resize(C);
		for (size_t i = 0; i < C; ++i) _scopes[i] = _factors[i].vars();
		const mex::vector<EdgeID>& elist = edges();
		for (size_t i = 0; i < elist.size(); ++i) {
			findex a,b;
			a = elist[i].first;
			b = elist[i].second;
			if (a > b) continue;
			VarSet sep = _factors[a].vars() & _factors[b].vars();
			_separators[a][b] = sep;
			_separators[b][a] = sep;
		}

		// incoming and outgoing
		_in.resize(C);
		_out.resize(C);
		for (vector<std::pair<findex, findex> >::const_iterator i = _schedule.begin();
				i != _schedule.end(); ++i) {
			findex from = (*i).first;
			findex to = (*i).second;
			_in[to] |= from;
			_out[from] |= to;
		}

		// init the root cluster(s)
		for (size_t i = 0; i < _out.size(); ++i) {
			if ( _out[i].empty() )
				_roots |= i;
		}

		// init forward and backward messages
		size_t N = _schedule.size();
		_forward.resize(N);
		_forward2.resize(N);
		_errors.resize(N, 0.0);
		_card.resize(N, 1);
		_sigmas.resize(N);
		_edgeIndex.resize(C);
		for (size_t i = 0; i < C; ++i) _edgeIndex[i].resize(C);
		for (size_t i = 0; i < N; ++i) {
			size_t from = _schedule[i].first;
			size_t to = _schedule[i].second;
			_edgeIndex[from][to] = i;
		}

		// init clique potentials
		_logZ = 0;
		for (size_t i = 0; i < _factors.size(); ++i) {
			_factors[i] = Factor(1.0); // init

			// clique potential
			for (flist::const_iterator j = _originals[i].begin();
					j != _originals[i].end(); ++j) {
				_factors[i] *= _gmo.factor(*j);
			}
		}

		// creating the indicator factors for the MAP variables
		size_t num_ind = 0;
		_deltas.resize(_gmo.nvar());
		for (size_t v = 0; v < _gmo.nvar(); ++v) {
			if (_varTypes[v] == true) {
				Var VX = var(v);
				for (size_t k = 0; k < VX.states(); ++k) {
					Factor f(VX, 0.0);
					for (size_t l = 0; l < f.numel(); ++l) {
						if (l == k) {
							f[l] = 1.0;
							break;
						}
					}

					findex fid = addFactor(f); // add this factor to _factors
					_deltas[v] |= fid;
					++num_ind;
				}
			}
		}

		std::cout << "Added " << num_ind << " indicator factors" << std::endl;

		/////// DEBUG purpose //////////////////////////////////////////////////
//		std::cout << "------DEBUG-------------------------------------------\n";
//		std::cout << "Bucket tree with " << _numClusters << " clusters and "
//				<< elist.size() << " edges" << std::endl;
//		for (size_t i = 0; i < elist.size(); ++i) {
//			findex a,b;
//			a = elist[i].first;
//			b = elist[i].second;
//			if (a > b) continue;
//			std::cout << "  edge from "
//					<< _scopes[a] << " to "
//					<< _scopes[b] << " (a=" << a << ", b=" << b << ")"
//					<< " sep: " << _separators[a][b]
//					<< std::endl;
//		}
//
//		std::cout << "Forward propagation schedule:" << std::endl;
//		for (size_t i = 0; i < _schedule.size(); ++i) {
//			std::cout << " msg " << _schedule[i].first << " --> "
//					<< _schedule[i].second << std::endl;
//		}
//
//		std::cout << "Original factors per cluster:" << std::endl;
//		for (size_t i = 0; i < _originals.size(); ++i) {
//			std::cout << " cl " << i << " : ";
//			std::copy(_originals[i].begin(), _originals[i].end(),
//					std::ostream_iterator<int>(std::cout, " "));
//			std::cout << std::endl;
//		}
//
//		std::cout << "Delta factors per MAP variable:" << std::endl;
//		for (size_t i = 0; i < _deltas.size(); ++i) {
//			if (_varTypes[i] == true) {
//				std::cout << " var " << i << " : ";
//				std::copy(_deltas[i].begin(), _deltas[i].end(),
//						std::ostream_iterator<int>(std::cout, " "));
//				std::cout << std::endl;
//			}
//		}
//
//		// _in and _out lists
//		std::cout << "_IN list:" << std::endl;
//		for (size_t i = 0; i < _in.size(); ++i) {
//			std::cout << "  _in[" << i << "] = ";
//			std::copy(_in[i].begin(), _in[i].end(),
//					std::ostream_iterator<int>(std::cout, " "));
//			std::cout << std::endl;
//		}
//		std::cout << "_OUT list:" << std::endl;
//		for (size_t i = 0; i < _out.size(); ++i) {
//			std::cout << "  _out[" << i << "] = ";
//			std::copy(_out[i].begin(), _out[i].end(),
//					std::ostream_iterator<int>(std::cout, " "));
//			std::cout << std::endl;
//		}
//		std::cout << "_ROOTS: ";
//		std::copy(_roots.begin(), _roots.end(),
//				std::ostream_iterator<int>(std::cout, " "));
//		std::cout << std::endl;
//
//		// factors, forward and backward
//		std::cout << "_factors:" << std::endl;
//		for (size_t i = 0; i < _factors.size(); ++i) {
//			std::cout << "[" << i << "]: " << _factors[i] << std::endl;
//		}
//		std::cout << "_forward2:" << std::endl;
//		for (size_t i = 0; i < _forward2.size(); ++i) {
//			std::cout << "(" << i << "): " << _forward2[i] << std::endl;
//		}
//
		////////////////////////////////////////////////////////////////////////
	}

	// compute the incoming messages to a cluster 'a'
	FactorSet incoming_mu(findex a) {

		FactorSet fs;
		const Factor& phi = _factors[a]; // clique factor
		size_t v = _cl2var[a]; // get the clique var
		if (_deltas[v].empty()) {
			fs.add(phi);
		} else {
			for (flist::const_iterator j = _deltas[v].begin();
					j != _deltas[v].end(); ++j) {
				Factor c = _factors[*j] * phi;
				fs.add(c);
			}
		}

		// incoming factorset messages to 'a'
		for (flist::const_iterator ci = _in[a].begin();
				ci != _in[a].end(); ++ci) {
			findex p = (*ci);
			size_t j = _edgeIndex[p][a];
			FactorSet& inc = _forward2[j];

			vector<Factor> tmp;
			vector<Factor>::const_iterator it1, it2;
			for (it1 = inc.factors().begin(); it1 != inc.factors().end(); ++it1) {
				const Factor& mu = (*it1);
				for (it2 = fs.factors().begin(); it2 != fs.factors().end(); ++it2) {
					Factor a = (*it2)*mu;
					tmp.push_back(a);
				}
			}

			fs.factors().swap(tmp);
		}

		return fs;
	}

	// compute the incoming sigma messages to a cluster 'a'
	FactorSet incoming_ub(findex a) {

		FactorSet fs;
		const Factor& phi = _factors[a]; // clique factor
		size_t v = _cl2var[a]; // get the clique var
		if (_deltas[v].empty()) {
			fs.add(phi);
		} else {
			for (flist::const_iterator j = _deltas[v].begin();
					j != _deltas[v].end(); ++j) {
				Factor c = _factors[*j] * phi;
				fs.add(c);
			}
		}

		// incoming factorset messages to 'a'
		for (flist::const_iterator ci = _in[a].begin();
				ci != _in[a].end(); ++ci) {
			findex p = (*ci);
			size_t j = _edgeIndex[p][a];
			FactorSet& inc = _sigmas[j];

			vector<Factor> tmp;
			vector<Factor>::const_iterator it1, it2;
			for (it1 = inc.factors().begin(); it1 != inc.factors().end(); ++it1) {
				const Factor& mu = (*it1);
				for (it2 = fs.factors().begin(); it2 != fs.factors().end(); ++it2) {
					Factor a = (*it2)*mu;
					tmp.push_back(a);
				}
			}

			fs.factors().swap(tmp);
		}

		return fs;
	}

	// Variable elimination: return lower and upper bounds
	std::pair<double, double> forward(int iter) {

		size_t max_el = 0;
		for (VarOrder::const_iterator x = _order.begin(); x != _order.end(); ++x) {

			//std::cout << "Eliminating " << *x << (_varTypes[*x] ? "(MAP)\n" : "(SUM)\n");

			// generate forward messages
			Var VX = var(*x);
			findex a = _clusters[*x][0]; // get source bucket of the variable
			if ( _out[a].size() > 0 ) {
				findex b = *(_out[a].begin());
				size_t i = _edgeIndex[a][b];

				FactorSet tmp = incoming_mu(a); max_el = std::max(max_el, tmp.size());
//				std::cout << "Incoming mu set: " << tmp.size() << " elements" << std::endl;
//				std::cout << tmp;
				tmp.sum(VX); // eliminate the bucket variable
//				std::cout << "After sum: " << tmp.size() << " elements" << std::endl;
//				std::cout << tmp;

				FactorSet sig = incoming_ub(a);
//				std::cout << "Incoming sig set: " << sig.size() << " elements" << std::endl;
//				std::cout << sig;
				sig.sum(VX); // eliminate the bucket variable
//				std::cout << "After sum: " << sig.size() << " elements" << std::endl;
//				std::cout << sig;

				// reduce the cardinality
//				std::cout << "*Pruning with card set: " << _card[i] << std::endl;
				tmp.prune_mu(_card[i], _forward2[i], _errors[i]);
				sig.prune_ub(_card[i], _sigmas[i]);

//				std::cout << "forward msg (" << a << "," << b << "): elim = " << VX << " -> " << std::endl;
//				std::cout << _forward2[i];
//				std::cout << "sigmas  msg (" << a << "," << b << "): elim = " << VX << " -> " << std::endl;
//				std::cout << _sigmas[i];
//				std::cout << "error found: " << _errors[i] << std::endl;
			}
		}

		// compute the lower and upper bounds
		double L = -infty(), U = infty();
		assert(_roots.size() == 1); // assume connected graph with 1 root only
		for (flist::const_iterator ci = _roots.begin();
				ci != _roots.end(); ++ci) {

			double err;
			size_t r = *ci; // root cluster
			FactorSet mu, sig, fs;
			FactorSet Lr = incoming_mu(r); // not pruned
			FactorSet Ur = incoming_ub(r); // not pruned
//			std::cout << "Lr set: " << Lr.size() << " elements" << std::endl;
//			std::cout << Lr;
//			std::cout << "Ur set: " << Ur.size() << " elements" << std::endl;
//			std::cout << Ur;

			size_t k = 0;
			for (flist::const_iterator it = _in[r].begin();
					it != _in[r].end(); ++it) {
				findex p = (*it);
				size_t j = _edgeIndex[p][r];
				k = std::max(k, _card[j]);
			}

			Lr.prune_mu(k, mu, err);
			Ur.prune_ub(k, sig);

			size_t v = _cl2var[r];
			Var VX = var(v);
			mu.sum(VX);

			double lb = mu.max(); // max over constants
			double ub = -infty(); // max over unary factors over V
			for (size_t i = 0; i < sig.factors().size(); ++i) {
				Factor& f = sig.factors()[i];
				ub = std::max(ub, f.max());
			}

			L = lb;
			U = ub;
		}

//		std::cout << "-------------------------------------------------------" << std::endl;
//		std::cout << "  AFSE Bounds: " << L.max() << " (" << std::log(L.max()) << ") "
//				<< U.max() << " (" << std::log(U.max()) << ") t="
//				<< (timeSystem() - _startTime) << " i=" << iter << " [" << max_el << "]"
//				<< std::endl;

		return std::make_pair(L, U);
	}

	// see if they're equal (in log space, up to 15 dec places)
	bool equals(double a, double b) {
		double d1 = std::log(a), d2 = std::log(b);
		double eps = fabs(d1 - d2);
		if (std::floor(std::log(eps)) <= -6.0) return true;
		else return false;
	}

	// tighten the MAP upper bound
	size_t solve(double stopTime, size_t step = 1) {

		size_t nIter = 0;
		double upbo = infty(), lowbo = -infty();
		_logZlb = -infty(); _logZub = infty();
		do {

			std::pair<double, double> bounds = forward(nIter);
			lowbo = std::max(lowbo, bounds.first);
			upbo = std::min(upbo, bounds.second);

			_logZlb = std::log(lowbo);
			_logZub = std::log(upbo);

			++nIter;

			// select the message with highest error and increase the cardinality
			int cand = -1;
			double err = -1.0;
			for (size_t i = 0; i < _errors.size(); ++i) {
				if (_errors[i] > err) {
					err = _errors[i];
					cand = i;
				}
			}

			assert(cand >= 0);
			_card[cand] += step;

//			std::cout << "**Updated card set to " << cand << " as " << _card[cand] << std::endl;
			std::cout << "  AFSE Bounds: " << lowbo << " (" << std::log(lowbo) << ") "
					<< upbo << " (" << std::log(upbo) << ") t="
					<< (timeSystem() - _startTime) << " i=" << nIter
					<< std::endl;

		} while ( lowbo < upbo );

		std::cout << "  AFSE Solution: " << lowbo << " (" << std::log(lowbo) << ") "
				<< upbo << " (" << std::log(upbo) << ") t="
				<< (timeSystem() - _startTime)
				<< std::endl;
		_logZ = std::log(lowbo);
		_logZlb = _logZub = _logZ; // found optimal solution

		return nIter;
	}
};

//////////////////////////////////////////////////////////////////////////////////////////////
}// namespace mex

#endif /* JGLP_H_ */
