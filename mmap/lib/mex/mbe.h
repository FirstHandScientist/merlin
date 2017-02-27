#ifndef __MEX_WEIGHTED_MINIBUCKET_H
#define __MEX_WEIGHTED_MINIBUCKET_H

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

class mbe: public graphModel, public gmAlg, virtual public mxObject {
public:
	typedef graphModel::findex findex;        // factor index
	typedef graphModel::vindex vindex;        // variable index
	typedef graphModel::flist flist;         // collection of factor indices

public:
	mbe() : graphModel() {
		setProperties();
	}
	mbe(const graphModel& gm) : graphModel(gm), _gmo(gm) {
		clearFactors();
		setProperties();
	}
	virtual mbe* clone() const {
		mbe* gm = new mbe(*this);
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
	}  // !!! init? or run?

	MEX_ENUM( Property , iBound,Order,Distance,DoMatch,DoWeight )
	;

	bool _byScope;
	bool _doMatch;
	bool _doWeight;
	Factor::Distance distMethod;
	graphModel::OrderMethod ordMethod;
	size_t _iBound;
	double _logZ;
	VarOrder _order;
	vector<flist> atElim;
	vector<double> atElimNorm;
	vector<vindex> _parents;
	vector<bool> _varTypes; // true if map, false if sum

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
	const vector<bool>& getVarTypes() {
		return _varTypes;
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
			setProperties("iBound=4,Order=MinWidth,DoMatch=0,DoWeight=1");
			_byScope = true;
			return;
		}
		std::vector<std::string> strs = mex::split(opt, ',');
		for (size_t i = 0; i < strs.size(); ++i) {
			std::vector<std::string> asgn = mex::split(strs[i], '=');
			switch (Property(asgn[0].c_str())) {
			case Property::iBound:
				setIBound(atol(asgn[1].c_str()));
				_byScope = true;
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
			case Property::DoWeight:
				_doWeight = atol(asgn[1].c_str());
				break;

			default:
				break;
			}
		}
	}

	// weighted elimination
	Factor elim(const Factor& F, const VarSet& vs, const double w) {
		if (w == infty()) {
			return F.max(vs);
		} else {
			return F.sumPower(vs, w);
		}
	}

	// weighted marginals
	Factor marg(const Factor& F, const VarSet& vs, const double w) {
		return F.marginal(vs, w);
	}

	// for the mini-bucket heuristic: lookup remaining cost given context
	template<class MapType>
	double logHeurToGo(Var v, MapType vals) const {
		double s = 0.0;
		for (size_t i = 0; i < atElim[_vindex(v)].size(); ++i) {
			findex ii = atElim[_vindex(v)][i];
			const VarSet& vs = factor(ii).vars();
			s += std::log(factor(ii)[sub2ind(vs, vals)]);
		}
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

	void init() {
		_logZ = 0.0;
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

		std::cout << "Use scope only? " << _byScope << "\n";
		assert(_byScope == true);
		_doMatch = true;
		_doWeight = true;
		std::cout << "Weighted Mini-Buckets? " << _doWeight << "\n";
		std::cout << "Moment Matching? " << _doMatch << "\n";

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

		atElim.clear();
		atElim.resize(_gmo.nvar());
		atElimNorm.clear();
		atElimNorm.resize(_gmo.nvar(), 0.0);

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

			//// Perform any matching? /////////////////////////////////////////////////
			//    "Matching" here is: compute the largest overlap of all buckets,
			//    and ensure that the moments on that subset of variables
			//    (ie, weighted marginals) are identical in all buckets.
			//
			//    Also, add extra variables if we can afford them?
			//size_t beta=0;

			//// Weight for mini-buckets /////////
			size_t select = ids[0];
			double weight = 1.0/((double)ids.size()); // uniform weights
			// check if bucket variable is a MAP variable
			if (_varTypes[*x]) weight = infty();

			if (_doMatch && ids.size() > 1) {

				if (_varTypes[*x]) { // match max marginals
					vector<Factor> ftmp(ids.size());       // compute geometric mean
					VarSet var = fin[ids[0]].vars();      // on all mutual variables
					for (size_t i = 1; i < ids.size(); i++)
						var &= fin[ids[i]].vars();
					Factor fmatch(var,1.0);
					for (size_t i = 0; i < ids.size(); i++) {
						ftmp[i] = marg(fin[ids[i]],var,weight);
						fmatch *= ftmp[i];
						//ftmp[i] = marg(fin[ids[i]], var).log();
						//fmatch += ftmp[i];
					}
					fmatch ^= (1.0/ids.size());         // and match each bucket to it
					//fmatch *= (1.0 / ids.size());     // and match each bucket to it
					for (size_t i=0;i<ids.size();i++)
						fin[ids[i]] *= (fmatch/ftmp[i]);

				} else { // match weighted marginals
					double W = 0.0;
					vector<Factor> ftmp(ids.size());       // compute geometric mean
					VarSet var = fin[ids[0]].vars();      // on all mutual variables
					for (size_t i = 1; i < ids.size(); i++)
						var &= fin[ids[i]].vars();
					Factor fmatch(var,1.0);
					//Factor fmatch(var, 0.0);
					for (size_t i = 0; i < ids.size(); i++) {
						ftmp[i] = marg(fin[ids[i]],var,weight);
						fmatch *= (ftmp[i]^weight);
						//fmatch *= ftmp[i];
						//ftmp[i] = marg(fin[ids[i]], var).log();
						//fmatch += ftmp[i];
						W += weight;
					}
					fmatch ^= (1.0/W);                  // and match each bucket to it
					//fmatch *= (1.0 / ids.size());     // and match each bucket to it
					for (size_t i=0;i<ids.size();i++)
						fin[ids[i]] *= ((fmatch/ftmp[i])^weight);

					//for (size_t i = 0; i < ids.size(); i++)
					//	fin[ids[i]] *= (fmatch - ftmp[i]).exp();

					//beta = addFactor( Factor(fmatch.vars(),1.0) );  // add node to new cluster graph
					//atElim[*x] |= beta;
				}
			}

			//// Eliminate individually within buckets /////////////////////////////////
			//   currently does not use weights except 0/1; !!! add sumPower alternatives from matlab code
//			vector<findex> alphas;
//			std::cout << "  # mini-buckets: " << ids.size() << std::endl;
//			std::cout << "    mini-buckets: "; for (flistIt i=ids.begin();i!=ids.end();++i) std::cout<<fin[*i]<<" "; std::cout<<"\n";
			for (flistIt i = ids.begin(); i != ids.end(); ++i) {
				//
				// Create new cluster alpha over this set of variables; save function parameters also
//				findex alpha = findex(-1), alpha2 = findex(-1);
				findex alpha2 = findex(-1);

//				std::cout << "    weight is " << weight << std::endl;
//				std::cout<< "  -> mb before: " << fin[*i] << std::endl;

				// weighted elimination
				if (_doWeight) {
					fin[*i] = elim(fin[*i], VX, weight);
				} else {
					if (_varTypes[*x] == false) { // SUM bucket
						if (*i == select) fin[*i] = elim(fin[*i], VX, 1.0);
						else fin[*i] = elim(fin[*i], VX, infty());
					} else { // MAX bucket
						fin[*i] = elim(fin[*i], VX, infty());
					}
				}
//				std::cout<< "  -> mb after:  " << fin[*i] << std::endl;

				// normalize for numerical stability
				double maxf = fin[*i].max();
				fin[*i] /= maxf;
				maxf = std::log(maxf);
				_logZ += maxf;
				Norm[*i] += maxf; // save normalization for overall bound

//				std::cout<< "  -> mb norm:   " << fin[*i] << std::endl;
				alpha2 = addFactor(fin[*i]);

//				Orig[*i].clear();
//				New[*i].clear();
//				New[*i] |= alpha;  // now incoming nodes to *i is just alpha

				size_t k = _parents[*x]; //  mark next bucket and intermediates with msg for heuristic calc
				for (; k != vindex(-1) && !fin[*i].vars().contains(var(k)); k =
						_parents[k]) {
					atElim[k] |= alpha2;
					atElimNorm[k] += Norm[*i];
				}
				if (k != vindex(-1)) {
					atElim[k] |= alpha2;
					atElimNorm[k] += Norm[*i];
				}  // need check?

				insert(vin, *i, fin[*i].vars()); // recompute and update adjacency
			}
//			std::cout<<"\n";

		}
		/// end for: variable elim order /////////////////////////////////////////////////////////

		Factor F(0.0);
		for (size_t i = 0; i < fin.size(); ++i) {
			F += log(fin[i]);
		}
		assert( F.nvar() == 0);
		_logZ += F.max();

		std::cout<<"MAP log Bound "<<_logZ<<"\n";

	}

	/*
	 *
	 WMBE: (1) forward pass MBE
	 (2) constructing cluster graph and connectivity

	 MPLP (1) "by variable" for any graphical model
	 (2) "by cluster set" as a generic function
	 - add "reverse schedule" for MBE tightening
	 (2) "by cluster" subgradient descent
	 - for each cluster, choose random maximizer and convert to vector form
	 - save by cluster id and push solution by variable id into "map" of vectors
	 - compute averages by variable id
	 - for each cluster, update parameters
	 (3) "by edge" in a factor graph, scheduled
	 (4) by tree block
	 (5) by monotone chains


	 *
	 */

};

//////////////////////////////////////////////////////////////////////////////////////////////
}// namespace mex

#endif  // re-include

