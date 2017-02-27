/*
 * jglp.h
 *
 *  Created on: Dec 19, 2013
 *      Author: radu
 */

#ifndef __MEX_MINIBUCKET_JGLP_H
#define __MEX_MINIBUCKET_JGLP_H

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

class jglp: public graphModel, public gmAlg, virtual public mxObject {
public:
	typedef graphModel::findex findex;        // factor index
	typedef graphModel::vindex vindex;        // variable index
	typedef graphModel::flist flist;         // collection of factor indices

public:
	jglp() : graphModel() {
		setProperties();
	}
	jglp(const graphModel& gm) : graphModel(gm), _gmo(gm) {
		clearFactors();
		setProperties();
	}
	virtual jglp* clone() const {
		jglp* gm = new jglp(*this);
		return gm;
	}

	graphModel _gmo;

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
		return _logZ;
	}
	double logZlb() const {
		return _logZ;
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

	MEX_ENUM( Property , iBound,Order,Distance,DoMatch,DoJG )
	;

	bool _byScope;
	bool _doMatch;
	bool _doJG;
	Factor::Distance distMethod;
	graphModel::OrderMethod ordMethod;
	size_t _iBound;
	double _logZ;
	VarOrder _order;
	vector<flist> atElim;
	vector<double> atElimNorm;
	vector<vindex> _parents;
	vector<bool> _varTypes; // _varTypes[i] true if MAX, false if SUM

	///////////////// JGLP local structures ///////////////////////////////////
	vector<bool> _types;			// the type of each cluster (SUM or MAX)
	vector<double> _weights;	// the weight of each cluster
	vector<flist> _match;		// clusters to be matched for each variable

	vector<flist> _originals;	// original factors (index) for each cluster
	vector<VarSet> _scopes;		// the scope (vars) for each cluster

	vector<vector<VarSet> > _separators;

	vector<flist> _in;			// incoming to each cluster
	vector<flist> _out; 		// outgoing from each cluster
	flist _roots;				// root cluster(s)

	vector<Factor> _forward; // forward messages (by edge)
	vector<Factor> _backward; // backward messages (by edge)
	vector<Factor> _reparam; // reparameterization function (by cluster)
	double _normR;   		// normalization constants (by cluster)
	vector<double> _norm;   // normalization constants (by edge)

	vector<pair<findex, findex> > _schedule;
	vector<vector<size_t> > _edgeIndex;

	double _startTime;

	/////////////////////////////////////////////////////////////////
	// Setting properties (directly or through property string)
	/////////////////////////////////////////////////////////////////

	void setIBound(size_t i) {
		_iBound = i ? i : std::numeric_limits<size_t>::max();
	}
	size_t getIBound() const {
		return _iBound;
	}
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
			setProperties("iBound=4,Order=MinWidth,DoMatch=1,DoJG=0");
			_byScope = true;
			return;
		}
		std::vector<std::string> strs = mex::split(opt, ',');
		for (size_t i = 0; i < strs.size(); ++i) {
			std::vector<std::string> asgn = mex::split(strs[i], '=');
			switch (Property(asgn[0].c_str())) {
			case Property::iBound:
				setIBound(atol(asgn[1].c_str()));
				break;
			case Property::Order:
				_order.clear();
				_parents.clear();
				ordMethod = graphModel::OrderMethod(asgn[1].c_str());
				break;
			case Property::Distance:
				distMethod = Factor::Distance(asgn[1].c_str());
				_byScope = false;
				break;
			case Property::DoMatch:
				_doMatch = atol(asgn[1].c_str());
				break;
			case Property::DoJG:
				_doJG = atol(asgn[1].c_str());
				break;
			default:
				break;
			}
		}
	}

	// weighted elimination
	Factor elim(const Factor& F, const VarSet& vs, const double w) {
		return F.sumPower(vs, w);
	}

	// weighted marginals
	Factor marg(const Factor& F, const VarSet& vs, const double w) {
		return F.marginal(vs, w);
	}

	// !!!! for Lars : lookup remaining cost given context
	template<class MapType>
	double logHeurToGo(Var v, MapType vals) const {
		double s = 0.0;
		for (size_t i = 0; i < atElim[_vindex(v)].size(); ++i) {
			findex ii = atElim[_vindex(v)][i];
			const VarSet& vs = _forward[ii].vars();
			s += std::log(_forward[ii][sub2ind(vs, vals)]);
		}

//		std::cout << "var " << v << " s = " << s + atElimNorm[_vindex(v)] << std::endl;
		return s + atElimNorm[_vindex(v)];
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

		// First downward pass to initialize the mini-bucket tree and backward messages
		_match.resize(_order.size());
		for (VarOrder::const_iterator x = _order.begin(); x != _order.end(); ++x) {

			//std::cout << "Eliminating "<<*x << (_varTypes[*x] ? "(MAP)\n" : "(SUM)\n");

			Var VX = var(*x);
			if (*x >= vin.size() || vin[*x].size() == 0)
				continue;  // check that we have some factors over this variable

			flist ids = vin[*x];  // list of factor IDs contained in this bucket

			//// Select allocation into mini-buckets ///////////////////////////
			typedef flist::const_iterator flistIt;
			typedef std::pair<double, sPair> _INS;
			std::multimap<double, sPair> scores;
			std::map<sPair, std::multimap<double, sPair>::iterator> reverseScore;

			//// Populate list of pairwise scores for aggregation //////////////
			for (flistIt i = ids.begin(); i != ids.end(); ++i) {
				for (flistIt j = ids.begin(); j != i; ++j) {
					double err = score(fin, VX, *i, *j);
					sPair sp(*i, *j);
					reverseScore[sp] = scores.insert(_INS(err, sp)); // save score
				}
				reverseScore[sPair(*i, *i)] = scores.insert(
						_INS(-1, sPair(*i, *i)));       // mark self index at -1
			}

			//// Run through until no more pairs can be aggregated: ////////////
			//   Find the best pair (ii,jj) according to the scoring heuristic and join
			//   them as jj; then remove ii and re-score all pairs with jj
			for (;;) {
				std::multimap<double, sPair>::reverse_iterator top =
						scores.rbegin();
				if (top->first < 0)
					break;                         // if can't do any more, quit
				else {
					size_t ii = top->second.first, jj = top->second.second;
					//std::cout<<"Joining "<<ii<<","<<jj<<"; size "<<(fin[ii].vars()+fin[jj].vars()).nrStates()<<"\n";
					fin[jj] |= fin[ii];                        // combine into j
					erase(vin, ii, fin[ii]);
					fin[ii] = VarSet();  //   & remove i

					Orig[jj] |= Orig[ii];
					Orig[ii].clear(); // keep track of list of original factors in this cluster
					New[jj] |= New[ii];
					New[ii].clear(); //  list of new "message" clusters incoming to this cluster

					for (flistIt k = ids.begin(); k != ids.end(); ++k) { // removing entry i => remove (i,k) for all k
						scores.erase(reverseScore[sPair(ii, *k)]);
					}
					ids /= ii;

					for (flistIt k = ids.begin(); k != ids.end(); ++k) { // updated j; rescore all pairs (j,k)
						if (*k == jj)
							continue;
						double err = score(fin, VX, jj, *k);
						sPair sp(jj, *k);
						scores.erase(reverseScore[sp]);    // change score (i,j)
						reverseScore[sp] = scores.insert(_INS(err, sp));  //
					}
				}
			}

			//// Weight for mini-buckets /////////
			double weight = 1.0/((double)ids.size()); // uniform weights
			// check if bucket variable is a MAP variable
			if (_varTypes[*x] == true) weight = infty();

			//// Eliminate individually each mini-bucket ///////////////////////
			vector<findex> alphas;
			for (flistIt i = ids.begin(); i != ids.end(); ++i) {
				//
				// Create new cluster alpha over this set of variables; save function parameters also
				findex alpha = findex(-1);
				alpha = addFactor(Factor(fin[*i]));
				alphas.push_back(alpha);
				_match[*x] |= alpha;

				fin[*i] = fin[*i] - VX;

				// add inter clusters edges
				for (flistIt j = New[*i].begin(); j != New[*i].end(); ++j) {
					addEdge(*j, alpha);
					_schedule.push_back(make_pair(*j, alpha));
				}

				// update cluster types
				_types.push_back(_varTypes[*x]);
				_weights.push_back(weight);

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

			//std::cout<<"\n";

		}
		/// end for: variable elim order ///////////////////////////////////////

		// separators and cluster scopes
		size_t C = _factors.size();
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
		for (vector<pair<findex, findex> >::const_iterator i = _schedule.begin();
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
		_backward.resize(N);
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

		// init reparam function
		_reparam.resize( _factors.size(), Factor(1.0) );
		//_normR.resize( _factors.size(), 0.0 );
		_normR = 0.0;

		/////// DEBUG purpose //////////////////////////////////////////////////
//		std::cout << "------DEBUG-------------------------------------------\n";
//		std::cout << "Join graph (mini-bucket-tree) with " << _factors.size() << " clusters and "
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
//		std::cout << "Backward propagation schedule:" << std::endl;
//		vector<pair<findex, findex> >::reverse_iterator ri = _schedule.rbegin();
//		for (; ri != _schedule.rend(); ++ri) {
//			std::cout << " msg " << ri->second << " --> "
//					<< ri->first << std::endl;
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
//		// _match list
//		std::cout << "_MATCH list:" << std::endl;
//		for (size_t i = 0; i < _match.size(); ++i) {
//			std::cout << "  var " << i << ": ";
//			std::copy(_match[i].begin(), _match[i].end(),
//					std::ostream_iterator<int>(std::cout, " "));
//			std::cout << std::endl;
//		}
//
//		// factors, forward and backward
//		std::cout << "_factors:" << std::endl;
//		for (size_t i = 0; i < _factors.size(); ++i) {
//			std::cout << "[" << i << "]: " << _factors[i] << std::endl;
//		}
//		std::cout << "_forward:" << std::endl;
//		for (size_t i = 0; i < _forward.size(); ++i) {
//			std::cout << "(" << i << "): " << _forward[i] << std::endl;
//		}
//		std::cout << "_backward:" << std::endl;
//		for (size_t i = 0; i < _backward.size(); ++i) {
//			std::cout << "(" << i << "): " << _backward[i] << std::endl;
//		}

		////////////////////////////////////////////////////////////////////////
	}

	// compute the belief of a cluster 'a'
	Factor belief(findex a) {

		Factor bel = _factors[a] * _reparam[a];

		// forward messages to 'a'
		for (flist::const_iterator ci = _in[a].begin();
				ci != _in[a].end(); ++ci) {
			findex p = (*ci);
			size_t j = _edgeIndex[p][a];
			bel *= _forward[j];
		}

		// backward message to 'a'
		for (flist::const_iterator ci = _out[a].begin();
				ci != _out[a].end(); ++ci) {
			findex p = (*ci);
			size_t j = _edgeIndex[a][p];
			bel *= _backward[j];
		}

		return bel;
	}

	// compute the belief of a cluster 'a', excluding 'b' from incoming msgs
	Factor incoming(findex a, size_t i) {

		Factor bel = _factors[a] * _reparam[a];

		// forward messages to 'a'
		for (flist::const_iterator ci = _in[a].begin();
				ci != _in[a].end(); ++ci) {
			findex p = (*ci);
			size_t j = _edgeIndex[p][a];
			bel *= _forward[j];
			_norm[i] += _norm[j];
		}

		return bel;
	}

	// compute the belief of a cluster 'a'
	Factor incoming(findex a) {

		Factor bel = _factors[a] * _reparam[a];

		// forward messages to 'a'
		for (flist::const_iterator ci = _in[a].begin();
				ci != _in[a].end(); ++ci) {
			findex p = (*ci);
			size_t j = _edgeIndex[p][a];
			bel *= _forward[j];
		}

		return bel;
	}

	// forward pass with matching
	void forward(double step) {

		double prevZ = _logZ;
		_logZ = 0;
		_normR = 0;
		_norm.resize(_forward.size(), 0.0);

		for (VarOrder::const_iterator x = _order.begin(); x != _order.end(); ++x) {

//			std::cout << "Eliminating "<<*x << (_varTypes[*x] ? "(MAP)\n" : "(SUM)\n");

			// moment-match the clusters of this bucket
			match(*x, step);

			// generate forward messages
			Var VX = var(*x);
			for (flist::const_iterator it = _match[*x].begin();
					it != _match[*x].end(); ++it) {
				findex a = (*it);
				if ( _out[a].size() > 0 ) {
					findex b = *(_out[a].begin());
					size_t i = _edgeIndex[a][b];

//					std::cout << "forward msg (" << a << "," << b << "): elim = " << VX << " -> ";

					Factor tmp = incoming(a, i);
					if (_varTypes[*x] == false) {
						_forward[i] = tmp.sumPowerEx(VX, 1.0/_weights[a]);
					} else {
						_forward[i] = tmp.max(VX);
					}

					// normalize for numerical stability
					double maxf = _forward[i].max();
					_forward[i] /= maxf;
					double lnmaxf = std::log(maxf);
					_logZ += lnmaxf;
//					std::cout << _forward[i] << std::endl;
				}
			}
		}

		// compute the upper bound
		Factor F(0.0);
		for (flist::const_iterator ci = _roots.begin();
				ci != _roots.end(); ++ci) {

			Factor bel = belief(*ci);
			F += log( bel.max() );
		}
		_logZ += F.max();

		std::cout << "JG Bound: " << _logZ << " (" << std::exp(_logZ) << ") "
				<< (timeSystem() - _startTime) << " d=" << (fabs(_logZ - prevZ))
				<< std::endl;

		// adding back the normR constant
//		double Rdist = std::exp(_normR/_reparam.size());
//		for (size_t a = 0; a < _reparam.size(); ++a) {
//			_reparam[a] *= Rdist;
//		}
	}

	// backward pass
	void backward(size_t iter) {

		// update backward messages
		vector<pair<findex, findex> >::reverse_iterator ri = _schedule.rbegin();
		for (; ri != _schedule.rend(); ++ri ) {

			// compute backward message m(b->a)
			findex a = (*ri).first;
			findex b = (*ri).second;
			size_t i = _edgeIndex[a][b]; // edge index

			VarSet VX = _scopes[b] - _separators[a][b];

//			std::cout << "backward msg (" << b << "," << a << "): elim = " << VX << " -> ";

			// compute the belief at b
			Factor bel = belief(b);

			if (_types[b] == false && _types[a] == false) { // SUM-SUM

				bel ^= 1.0/_weights[b];
				bel /= (_forward[i]^(1.0/_weights[a])); // divide out m(a->b)

				//_backward[i] = elim(bel, VX, 1);
				_backward[i] = bel.sum(VX);
				_backward[i] ^= (_weights[b]);

			} else if (_types[b] == true && _types[a] == true) { // MAX-MAX

				bel /= _forward[i]; // divide out m(a->b)
				_backward[i] = bel.max(VX);

			} else if (_types[b] == true && _types[a] == false) { // MAX-SUM

				bel = bel.sigma(iter); // the weird sigma operator that focuses on max
				bel /= (_forward[i]^(1.0/_weights[a])); // divide out m(a->b)

				_backward[i] = bel.sum(VX);
				_backward[i] ^= (_weights[a]);

			} else {
				assert(false); // cannot reach this case!!
			}

//			std::cout << _backward[i] << std::endl;

			// normalize
			double mx = _backward[i].max();
			_backward[i] /= mx;
		}

	}

	// match the clusters of a bucket
	void match(size_t x, double step) {

		if (_match[x].size() <= 1)
			return; // no matching

		Var VX = var(x);
		if (_varTypes[x] == true) { // max marginals matching

//			std::cout << "matching max marginals" << std::endl;

			size_t R = _match[x].size();
			vector<Factor> ftmp(R); // compute geometric mean

			VarSet var;
			var |= VX; // on mutual variable (bucket variable)
			Factor fmatch(var,1.0);

			size_t i = 0;
			for (flist::const_iterator it = _match[x].begin();
					it != _match[x].end(); ++it, i++) {

				findex a = (*it);
				Factor bel = belief(a);
				ftmp[i] = bel.maxmarginal(var); // max-marginal
				fmatch *= ftmp[i];
			}

			i = 0;
			fmatch ^= (1.0/R); // and match each bucket to it
			for (flist::const_iterator it = _match[x].begin();
					it != _match[x].end(); ++it, ++i) {
				findex a = (*it);
				_reparam[a] *= (fmatch/ftmp[i]);

//				std::cout << " reparam      : " << _reparam[a] << std::endl;
//				double mx = _reparam[a].max();
//				_reparam[a] /= mx;
//				double lnmax = std::log(mx);
//				_normR += lnmax;
			}

		} else { // weighted marginals matching

//			std::cout << "matching weighted marginals on cliques: ";
//			std::copy(_match[x].begin(), _match[x].end(),
//					std::ostream_iterator<int>(std::cout, " "));
//			std::cout << std::endl;

			size_t R = _match[x].size();
			vector<Factor> ftmp(R);   // compute geometric mean

			VarSet var;
			var |= VX; // on mutual variable (bucket variable)
			Factor fmatch(var,1.0);
//			double W = 0.0;
			size_t i = 0;

			for (flist::const_iterator it = _match[x].begin();
					it != _match[x].end(); ++it, i++) {

				findex a = (*it);
				Factor bel = belief(a);
				bel ^= (1.0/_weights[a]);
				ftmp[i] = bel.marginal(var);
				fmatch *= (ftmp[i] ^ _weights[a]);
//				W += _weights[a];

//				std::cout << " clique belief: " << bel << std::endl;
//				std::cout << " marginal     : " << ftmp[i] << std::endl;
			}

			i = 0;
//			fmatch ^= (1.0/W);         // and match each bucket to it

//			std::cout << " geom mean    : " << fmatch << std::endl;
			for (flist::const_iterator it = _match[x].begin();
					it != _match[x].end(); ++it, ++i) {
				findex a = (*it);
				_reparam[a] *= ((fmatch/ftmp[i])^(step*_weights[a]));

//				std::cout << " reparam      : " << _reparam[a] << std::endl;

//				double mx = _reparam[a].max();
//				_reparam[a] /= mx;
//				double lnmax = std::log(mx);
//				_normR += lnmax;
			}
		}
	}

	// tighten the MAP upper bound
	void tighten(size_t nIter, double stopTime = -1, double stopObj = -1) {

		for (size_t iter = 1; iter <= nIter; ++iter) {
			double step = 1.0/(double)iter;

			double prevZ = _logZ;
			forward(step);
			backward(iter);

			double d = fabs(_logZ - prevZ);
			if (d < stopObj) break;

			// do at least one iterations
			if (stopTime > 0 && stopTime <= (timeSystem() - _startTime))
				break;

		}
	}

	// tighten the MAP upper bound and keep track of the argmax in MAX buckets
	std::vector<int> propagate(size_t iter) {

		double step = 1.0/(double)iter;

		// propagate messages
		forward(step);
		backward(iter);

		// recover assignment to MAX variables
		return recoverConfig();
	}

	// recover the argmax of the MAP variables
	std::vector<int> recoverConfig() {
		std::vector<int> sol;
		sol.resize(nvar(), -1);

		for (VarOrder::reverse_iterator x = _order.rbegin();
				x != _order.rend(); ++x) {

			if (_varTypes[*x] == false) {
				break; // summation variable so we're done.
			}

			// moment-match the clusters of this bucket
			Var VX = var(*x);

			VarSet var;
			var |= VX; // on mutual variable (bucket variable)
			Factor fmatch(var,1.0);
			for (flist::const_iterator it = _match[*x].begin();
					it != _match[*x].end(); ++it) {

				findex a = (*it);
				Factor bel = incoming(a);

				// project the current assignment to max variables
				VarSet vs = bel.vars();
				for (VarSet::const_iterator vi = vs.begin(); vi != vs.end(); ++vi) {
					Var y = (*vi);
					if (sol[y] != -1) {
						bel = bel.condition(y, sol[y]);
					}
				}

				fmatch *= bel.maxmarginal(var); // max-marginal
			}

			// get the argmax
			size_t xtmp = fmatch.argmax();
			sol[*x] = xtmp;
		}

		return sol;
	}

	// final forward pass to build the search heuristic
	void finalize(int iter) {

		_logZ = 0;

		atElim.clear();
		atElim.resize(_gmo.nvar());
		atElimNorm.clear();
		atElimNorm.resize(_gmo.nvar(), 0.0);
		_norm.resize(_forward.size(), 0.0);

		for (VarOrder::const_iterator x = _order.begin(); x != _order.end(); ++x) {

//			std::cout << "Eliminating "<<*x << (_varTypes[*x] ? "(MAP)\n" : "(SUM)\n");

			// moment-match the clusters of this bucket
			match(*x, 1.0/(double)iter);

			// generate forward messages
			Var VX = var(*x);
			for (flist::const_iterator it = _match[*x].begin();
					it != _match[*x].end(); ++it) {

				findex a = (*it);
				if ( _out[a].size() > 0 ) {
					findex b = *(_out[a].begin());
					size_t i = _edgeIndex[a][b];

					Factor tmp = incoming(a, i);
					if (_varTypes[*x] == false) {
						_forward[i] = tmp.sumPowerEx(VX, 1.0/_weights[a]);
					} else {
						_forward[i] = tmp.max(VX);
					}

					// normalize for numerical stability
					double maxf = _forward[i].max();
					_forward[i] /= maxf;
					double lnmaxf = std::log(maxf);
					_logZ += lnmaxf;
					_norm[i] += lnmaxf;

					//  mark next bucket and intermediates with msg for heuristic calc
					size_t k = _parents[*x];
					for (; k != vindex(-1) && !_forward[i].vars().contains(var(k)); k =
							_parents[k]) {
						atElim[k] |= i;
						atElimNorm[k] += _norm[i];
					}
					if (k != vindex(-1)) {
						atElim[k] |= i;
						atElimNorm[k] += _norm[i];
					}  // need check?

				}
			}
		}

		// compute the upper bound
		Factor F(0.0);
		for (flist::const_iterator ci = _roots.begin();
				ci != _roots.end(); ++ci) {

			Factor bel = belief(*ci);
			F += log( bel.max() );
		}
		_logZ += F.max();

		std::cout << "JG Bound: " << _logZ << " (" << std::exp(_logZ) << ") " << (timeSystem() - _startTime) << std::endl;
		std::cout << "JG Time:  " << (timeSystem() - _startTime) << " seconds" << std::endl;

		// reparameterize factors
		for (size_t a = 0; a < _factors.size(); ++a) {
			_factors[a] *= _reparam[a];
		}
	}

};

//////////////////////////////////////////////////////////////////////////////////////////////
}// namespace mex

#endif /* JGLP_H_ */
