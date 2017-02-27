/*
 * wmb.h
 *
 *  Created on: Oct 16, 2014
 *      Author: radu
 */

#ifndef WMB_H_
#define WMB_H_

#include <assert.h>
#include <stdexcept>
#include <stdlib.h>
#include <stdint.h>
#include <cstdlib>
#include <cstring>

#include "factorgraph.h"
#include "alg.h"

namespace mex {

// Weighted mini-buckets (WMB)
//   elimination is done by sumPower
//   the moment matching is done via weighted marginals
//

class wmb: public graphModel, public gmAlg, virtual public mxObject {
public:
	typedef graphModel::findex findex;        // factor index
	typedef graphModel::vindex vindex;        // variable index
	typedef graphModel::flist flist;         // collection of factor indices

public:
	wmb() : graphModel() {
		setProperties();
	}
	wmb(const graphModel& gm) : graphModel(gm), _gmo(gm) {
		clearFactors();
		setProperties();
	}
	virtual wmb* clone() const {
		wmb* gm = new wmb(*this);
		return gm;
	}

	graphModel _gmo;

	// Can be an optimization algorithm or a summation algorithm....
	double ub() const {
		return _logZ;
	}
	double lb() const {
		return _logZ;
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
	}  // !!! init? or run?

	bool _byScope;
	Factor::Distance distMethod;
	graphModel::OrderMethod ordMethod;
	size_t _iBound;
	double _logZ;
	VarOrder _order;
	vector<vindex> _parents;
	size_t _inducedWidth;
	bool _negativeWeights;	// set negative weights for computing a lower bound

	///////////////// JG local structures ///////////////////////////////////
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
	size_t getInducedWidth() const {
		return _inducedWidth;
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
		ordMethod = graphModel::OrderMethod::MinFill;
	}

	// weighted elimination
	Factor elim(const Factor& F, const VarSet& vs, const double w) {
		if (w == infty()) {
			return F.max(vs);
		} else if (w == -infty()) {
			return F.min(vs);
		} else {
			return F.sumPower(vs, w);
		}
	}

	// weighted marginals
	Factor marg(const Factor& F, const VarSet& vs, const double w) {
		return F.marginal(vs, w);
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

	void preprocess() {
		_order = _gmo.order(ordMethod);
		_parents = _gmo.pseudoTree(_order);
		_inducedWidth = _gmo.inducedWidth(_order);
	}

	void init(const VarSet& vs) {
		init();
	}            // !!! inefficient

	void init() {
		_logZ = 0.0;
		_inducedWidth = 0;
		_negativeWeights = true;
		if (_order.size() == 0) { // if we need to construct an elimination ordering
//			double tic = timeSystem();
			_order = _gmo.order(ordMethod);
			_parents.clear(); // (new elim order => need new pseudotree) !!! should do together
//			std::cout << "Order in " << timeSystem() - tic << " sec\n";
		}
		if (_parents.size() == 0) {     // if we need to construct a pseudo-tree
//			double tic = timeSystem();
			_parents = _gmo.pseudoTree(_order);
//			std::cout << "Pseudo in " << timeSystem() - tic << " sec\n";
		}

		// Get the factors and normalize them
		vector<Factor> fin(_gmo.factors());
		vector<double> Norm(_gmo.nFactors(), 0.0);
		for (size_t i = 0; i < _gmo.nFactors(); ++i) {
			double mx = fin[i].max();
			fin[i] /= mx;
			Norm[i] = std::log(mx);
			_logZ += Norm[i];
		}
		vector<flist> vin;

		for (size_t i = 0; i < _gmo.nvar(); ++i) {
			vin.push_back(_gmo.withVariable(var(i)));
		}

		//// Eliminate each variable in the sequence given: ////////////////////
		for (VarOrder::const_iterator x = _order.begin(); x != _order.end(); ++x) {

//			std::cout << "Eliminating "<<*x << (_varTypes[*x] ? "(MAP)\n" : "(SUM)\n");

			Var VX = var(*x);
			if (*x >= vin.size() || vin[*x].size() == 0)
				continue;  // check that we have some factors over this variable

			flist ids = vin[*x];  // list of factor IDs contained in this bucket

			//// Select allocation into buckets ///////////////////////////////////////
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

			//// Run through until no more pairs can be aggregated: ////////////////////
			//   Find the best pair (ii,jj) according to the scoring heuristic and join
			//   them as jj; then remove ii and re-score all pairs with jj
			for (;;) {
				std::multimap<double, sPair>::reverse_iterator top =
						scores.rbegin();
				//multimap<double,_IDX>::reverse_iterator  last=scores.lower_bound(top->first);  // break ties randomly !!!
				//std::advance(last, randi(std::distance(top,last)));
				//std::cout<<top->first<<" "<<top->second.first<<" "<<top->second.second<<"\n";

				if (top->first < 0)
					break;                         // if can't do any more, quit
				else {
					size_t ii = top->second.first, jj = top->second.second;
					//std::cout<<"Joining "<<ii<<","<<jj<<"; size "<<(fin[ii].vars()+fin[jj].vars()).nrStates()<<"\n";

					fin[jj] *= fin[ii];                      // combine into j
					Norm[jj] += Norm[ii];
					double mx = fin[jj].max();
					fin[jj] /= mx;
					mx = std::log(mx);
					_logZ += mx;
					Norm[jj] += mx;
					erase(vin, ii, fin[ii].vars());
					fin[ii] = Factor();  //   & remove i

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
			double R = (double)ids.size();
			double w0 = (2*R-1)/R;
			double w = -1/R;
			double weight;

			//// Eliminate individually within buckets /////////////////////////////////
			//   currently does not use weights except 0/1; !!! add sumPower alternatives from matlab code
//			vector<findex> alphas;
//			std::cout << "  # mini-buckets: " << ids.size() << std::endl;
//			std::cout << "    mini-buckets: "; for (flistIt i=ids.begin();i!=ids.end();++i) std::cout<<fin[*i]<<" "; std::cout<<"\n";
			int pos = 0;
			for (flistIt i = ids.begin(); i != ids.end(); ++i) {
				//
				// Create new cluster alpha over this set of variables; save function parameters also
//				findex alpha = findex(-1), alpha2 = findex(-1);
				findex alpha2 = findex(-1);

				// weighted elimination
				if (pos == 0) weight = w0;
				else weight = w;
				fin[*i] = elim(fin[*i], VX, weight);

//				std::cout<< "  -> mb after:  " << fin[*i] << std::endl;

				// normalize for numerical stability
				double maxf = fin[*i].max();
				fin[*i] /= maxf;
				maxf = std::log(maxf);
				_logZ += maxf;
				Norm[*i] += maxf; // save normalization for overall bound

//				std::cout<< "  -> mb norm:   " << fin[*i] << std::endl;
				alpha2 = addFactor(fin[*i]);

				insert(vin, *i, fin[*i].vars()); // recompute and update adjacency
				pos++;
			}
//			std::cout<<"\n";

		}
		/// end for: variable elim order /////////////////////////////////////////////////////////

		Factor F(0.0);
		for (size_t i = 0; i < fin.size(); ++i) {
			F += log(fin[i]);
		}
		assert( F.nvar() == 0);
		_logZ += F.sum();

//		std::cout<<"SUM log Bound "<<_logZ<<"\n";

	}

	// tighten the lower bound by iterative message passing along the join-graph
	// create the mini-bucket tree (join graph): symbols only
	void initJG(bool negWeights = true) {
		_startTime = timeSystem();
		_negativeWeights = negWeights;
		if (_order.size() == 0) { // if we need to construct an elimination ordering
			//double tic = timeSystem();
			_order = _gmo.order(ordMethod);
			_parents.clear(); // (new elim order => need new pseudotree) !!! should do together
			//std::cout << "Order in " << timeSystem() - tic << " sec\n";
		}
		if (_parents.size() == 0) {     // if we need to construct a pseudo-tree
			//double tic = timeSystem();
			_parents = _gmo.pseudoTree(_order);
			//std::cout << "Pseudo in " << timeSystem() - tic << " sec\n";
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
			double R = (double)ids.size();
			double w0 = (2*R-1)/R;
			double w = -1/R;
			double weight;

			//// Eliminate individually each mini-bucket ///////////////////////
			vector<findex> alphas;
			int pos = 0;
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
				if (_negativeWeights) weight = (pos == 0) ? w0 : w; // negative weights
				else weight = fabs(w); // positive weights
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
				pos++;
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
	Factor _belief(findex a) {

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

//		double prevZ = _logZ; belief
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
					_forward[i] = tmp.sumPower(VX, 1.0/_weights[a]);

					// normalize for numerical stability
					double maxf = _forward[i].max();
					_forward[i] /= maxf;
					double lnmaxf = std::log(maxf);
					_logZ += lnmaxf;
//					std::cout << _forward[i] << std::endl;
				}
			}
		}

		// compute the lower bound
		Factor F(0.0);
		for (flist::const_iterator ci = _roots.begin();
				ci != _roots.end(); ++ci) {

			Factor bel = _belief(*ci);
			F += log(bel.sum());
		}
		_logZ += F.sum();

//		std::cout << "JG Bound: " << _logZ << " (" << std::exp(_logZ) << ") "
//				<< (timeSystem() - _startTime) << " d=" << (fabs(_logZ - prevZ))
//				<< std::endl;

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
			Factor bel = _belief(b);
			bel ^= 1.0/_weights[b];
			bel /= (_forward[i]^(1.0/_weights[a])); // divide out m(a->b)

			//_backward[i] = elim(bel, VX, 1);
			_backward[i] = bel.sum(VX);
			_backward[i] ^= (_weights[b]);

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

//			std::cout << "matching weighted marginals on cliques: ";
//			std::copy(_match[x].begin(), _match[x].end(),
//					std::ostream_iterator<int>(std::cout, " "));
//			std::cout << std::endl;

		size_t R = _match[x].size();
		vector<Factor> ftmp(R);   // compute geometric mean

		VarSet var;
		var |= VX; // on mutual variable (bucket variable)
		Factor fmatch(var,1.0);
		size_t i = 0;

		for (flist::const_iterator it = _match[x].begin();
				it != _match[x].end(); ++it, i++) {

			findex a = (*it);
			Factor bel = _belief(a);
			bel ^= (1.0/_weights[a]);
			ftmp[i] = bel.marginal(var);
			fmatch *= (ftmp[i] ^ _weights[a]);
		}

		i = 0;

//		std::cout << " geom mean    : " << fmatch << std::endl;
		for (flist::const_iterator it = _match[x].begin();
				it != _match[x].end(); ++it, ++i) {
			findex a = (*it);
			_reparam[a] *= ((fmatch/ftmp[i])^(step*_weights[a]));

//			std::cout << " reparam      : " << _reparam[a] << std::endl;
		}
	}

	// tighten the lower bound
	void tighten(size_t nIter, double stopTime = -1, double stopObj = -1) {

		double bestZ = (_negativeWeights) ? -infty() : infty();
		std::cout << "WMB: tightening the bound:" << std::endl;
		for (size_t iter = 1; iter <= nIter; ++iter) {
			double step = 1.0/(double)iter; //1.0;//

//			double prevZ = _logZ;
			forward(step);
			backward(iter);
			if (_negativeWeights) bestZ = max(bestZ, _logZ);
			else bestZ = min(bestZ, _logZ);

//			double d = fabs(_logZ - prevZ);
//			if (d < stopObj) break;

			std::cout << "  logZ is " << _logZ << std::endl;

			// do at least one iterations
			if (stopTime > 0 && stopTime <= (timeSystem() - _startTime))
				break;

		}

		_logZ = bestZ;
	}

	// For local debugging of the bounds
	void test() {
		std::cout << " .. generating initial factors" << std::endl;

		vector<VarSet> fin;
		for (vector<Factor>::const_iterator i = _gmo.factors().begin();
				i != _gmo.factors().end(); ++i) {
			fin.push_back((*i).vars());
		}

		for (vector<VarSet>::const_iterator i = fin.begin(); i != fin.end(); ++i) {
			std::cout << (*i) << std::endl;
		}

		Factor Fab = _gmo.factor(0);
		Factor Fac = _gmo.factor(1);
		Factor Fad = _gmo.factor(2);
		Factor Fae = _gmo.factor(3);
		Factor Fbc = _gmo.factor(4);
		Factor Fbd = _gmo.factor(5);
		Factor Fbe = _gmo.factor(6);
		Factor Fbf = _gmo.factor(7);
		Factor Fce = _gmo.factor(8);
		Factor Fcf = _gmo.factor(9);
		Factor Fde = _gmo.factor(10);
		Factor Fef = _gmo.factor(11);

		std::cout << "Initial factors:" << std::endl;
		std::cout << Fab << "\n" << Fac << "\n" << Fad << "\n" << Fae
				<< "\n" << Fbc << "\n" << Fbd << "\n" << Fbe << "\n" << Fbf
				<< "\n" << Fce << "\n" << Fcf << "\n" << Fde << "\n" << Fef << std::endl;

		Factor F = Fab*Fac*Fad*Fae*Fbc*Fbd*Fbe*Fbf*Fce*Fcf*Fde*Fef;
		double exact = F.sum();
		std::cout << "Exact value is: " << exact << std::endl;

		// Upward messages
		Factor Uace_ad, Uceb_bed, Uceb_bf, Uace_ab, Uce_ceb, Uce_ace, Uce_cef;

		// Reparameterization
		Factor Rad, Rbed, Rbf, Rcef, Rab, Rceb;

		double wd1 = .5, wd2 = .5;
		double wf1 = .5, wf2 = .5;
		double wb1 = .5, wb2 = .5;

		std::cout << "Bounds:" << std::endl;
		double bestBound = infty();
		for (size_t iter = 1; iter <= 30; ++iter) {
			//Minibucket
			// Eliminate D:
			Var D = var(3);
			Factor Cad = Fad*Rad*Uace_ad, Cbed = Fbd*Fde*Rbed*Uceb_bed; // compute clique belief
			Factor Md1 = Cad.marginal(VarSet(D)), Md2 = Cbed.marginal(VarSet(D)); //compute reparameterization functions R
			Rad = (Rad*(Md1*Md2)^.5)/Md1; Rbed = (Rbed*(Md1*Md2)^.5)/Md2;
			Factor Mad_ace = (Fad*Rad).sumPower(VarSet(D), 1/wd1);     //compute downward messages
			Factor Mbed_ceb = (Fbd*Fde*Rbed).sumPower(VarSet(D), 1/wd2);

			// Eliminate F:
			Var F = var(5);
			Factor Cbf = Fbf*Rbf*Uceb_bf,  Ccef = Fcf*Fef*Rcef*Uce_cef; // compute clique beliefs
			Factor Mf1 = Cbf.marginal(VarSet(F)), Mf2 = Ccef.marginal(VarSet(F));       //compute reparam functions R
			Rbf = (Rbf*(Mf1*Mf2)^.5)/Mf1;  Rcef = (Rcef*(Mf1*Mf2)^.5)/Mf2;
			Factor Mbf_ceb  = (Fbf*Rbf).sumPower(VarSet(F), 1/wf1 );      //compute downward messages
			Factor Mcef_ce  = (Fcf*Fef*Rcef).sumPower(VarSet(F), 1/wf2);
			// Eliminate B:
			Var B = var(1);
			Factor Cab = Fab*Rab*Uace_ab,  Cceb = Fbc*Fbe*Rceb*Mbf_ceb*Mbed_ceb*Uce_ceb;  // compute clique beliefs
			Factor Mb1 = Cab.marginal(VarSet(B)), Mb2 = Cceb.marginal(VarSet(B));       //compute reparam functions R
			Rab = (Rab*(Mb1*Mb2)^.5)/Mb1; Rceb = (Rceb*(Mb1*Mb2)^.5)/Mb2;
			Factor Mab_ace  = (Fab*Rab).sumPower(VarSet(B), 1/wb1 );          // compute down msgs
			Factor Mceb_ce  = (Fbc*Fbe*Mbf_ceb*Mbed_ceb*Rceb).sumPower(VarSet(B), 1/wb2);

			// Eliminate A,C,E:
			Var A = var(0), C = var(2), E = var(4);
			Factor Mace_ce  = (Fac*Fae*Mab_ace*Mad_ace).sum(VarSet(A));         // compute down msg
			Factor Mce_e    = (Fce*Mace_ce*Mceb_ce*Mcef_ce).sum(VarSet(C));    // ""
			Factor BoundMB  = (Mce_e).sum(VarSet(E));                          // ""
			std::cout << iter << ":    bound is   " << BoundMB.max() << std::endl;
			bestBound = min(bestBound, BoundMB.max());


			// Now upward messages

			// Ue_ce : nothing (no factors or info below clique CE)
			VarSet CE(C,E);

			Uce_ace = (Fce*Mace_ce*Mceb_ce*Mcef_ce) / Mace_ce;      // sum-to-sum
			Uce_ceb = (((Fce*Mace_ce*Mceb_ce*Mcef_ce)/Mceb_ce)^(1/wb2))^(wb2);    //sum-to-sum
			Uce_cef = (((Fce*Mace_ce*Mceb_ce*Mcef_ce)/Mcef_ce)^(1/wf2))^(wf2);    //

			Uace_ab = (((( (Fac*Fae*Mab_ace*Mad_ace*Uce_ace).sum(CE)) / Mab_ace)^(1/wb1)).sum(CE))^(wb1);
			Uace_ad = (((( (Fac*Fae*Mab_ace*Mad_ace*Uce_ace).sum(CE)) / Mad_ace)^(1/wd1)).sum(CE))^(wd1);

			Uceb_bf = (((((Fbc*Fbe*Mbf_ceb*Mbed_ceb*Rceb*Uce_ceb)^(1/wb2)) / Mbf_ceb)^(1/wf1)).sum(CE))^wf1;
			Uceb_bed= (((((Fbc*Fbe*Mbf_ceb*Mbed_ceb*Rceb*Uce_ceb)^(1/wb2)) / Mbf_ceb)^(1/wd2)).sum(CE))^wd2;

		}

		std::cout << "Found best bound: " << bestBound << std::endl;
		std::cout << " .. done." << std::endl;
	} // done

};

//////////////////////////////////////////////////////////////////////////////////////////////
}// namespace mex

#endif /* WMB_H_ */
