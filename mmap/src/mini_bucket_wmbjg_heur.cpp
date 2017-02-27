/*
 * mini_bucket_mplp_heur.cpp
 *
 *  Created on: Oct 25, 2013
 *      Author: radu
 */



#include "mini_bucket_wmbjg_heur.h"
#include "timer.h"

// Constructor: allow for blank entries?  if so, how to respecify later?
// MPLP function (with or without elimination ordering?)
//    either static function that looks like constructor, or enable blank constructor?
//    construct mplp factorgraph, pass messages (stop crit specified), write back to problem
// JGLP function

// Copy DaoOpt Function class into mex::Factor class structures
mex::vector<mex::Factor> MiniBucketJGHeur::copyFactors(void) {
	mex::vector<mex::Factor> fs(_p->getC());
	for (int i = 0; i < _p->getC(); ++i)
		fs[i] = _p->getFunctions()[i]->asFactor();
	return fs;
}

// Mini-bucket may have re-parameterized the original functions; if so, replace them
void MiniBucketJGHeur::rewriteFactors(const vector<mex::Factor>& factors) {
	vector<Function*> newFunctions(factors.size()); // to hold replacement, reparameterized functions
	for (size_t f = 0; f < factors.size(); ++f) { // allocate memory, copy variables into std::set
		double* tablePtr = new double[factors[f].nrStates()];
		std::set<int> scope;
		for (mex::VarSet::const_iterator v = factors[f].vars().begin();
				v != factors[f].vars().end(); ++v)
			scope.insert(v->label());
		newFunctions[f] = new FunctionBayes(f, _p, scope, tablePtr,
				factors[f].nrStates());

		//newFunctions[f]->fromFactor(log(factors[f])); // write in log factor functions
		newFunctions[f]->fromFactor(factors[f]);
	}
	_p->replaceFunctions(newFunctions); // replace them in the problem definition
}

bool MiniBucketJGHeur::preprocess(const vector<val_t>* assignment) {
	return false;
}

size_t MiniBucketJGHeur::build(const vector<val_t>* assignment,
		bool computeTables) {
	assert(_pt);
	// we need a pseudo tree

	if (_options == NULL)
		std::cout << "Warning (MBE-ATI): ProgramOptions not available!\n";

	Timer tm(_options->jglps);
	tm.start();
	_jglp.init();

	size_t numIter = (_options->jglp > 0 ? (size_t)_options->jglp : 100);

	// message passing
	for (size_t iter = 1; iter <= numIter; ++iter) {

		std::vector<int> sol = _jglp.propagate(iter);

		// reconstruct full solution (original variables)
		int nOrg = m_problem->getNOrg();
		vector<val_t> curSolution;
		curSolution.resize(nOrg, UNKNOWN);

		const map<int, int>& old2new = m_problem->getOld2New();
		const map<int, int>& evidence = m_problem->getEvidence();
		for (int i = 0; i < nOrg; ++i) {
			map<int, int>::const_iterator itRen = old2new.find(i);
			if (itRen != old2new.end()) { // var part of solution
				curSolution.at(i) = sol.at(itRen->second);
			} else {
				map<int, val_t>::const_iterator itEvid = evidence.find(i);
				if (itEvid != evidence.end())  // var part of evidence
					curSolution.at(i) = itEvid->second;
				else
					// var had unary domain
					curSolution.at(i) = 0;
			}
		}

		// save the assignment to the output file
		std::string output = m_options->outputFile;
		if (output.empty() == false) {
			std::ofstream out;
			out.open(output.c_str(), std::ofstream::out | std::ofstream::app);
			if (out.fail()) {
				std::cerr << "Cannot append to output file " << output << std::endl;
				exit(EXIT_FAILURE);
			}

			out << "-BEGIN-" << std::endl;
			const std::set<int>& query = m_problem->getQuery();
			std::set<int>::const_iterator li = query.begin();
			out << query.size();
			for (; li != query.end(); ++li) {
				out << " " << (*li) << " " << curSolution[*li];
			}
			out << std::endl;
			out.close();
		}

		// check for timeout
		if (tm.timeout()) break;
	}

	return this->getSize();
}

/* finds a dfs order of the pseudotree (or the locally restricted subtree)
 * and writes it into the argument vector */
void MiniBucketJGHeur::findDfsOrder(vector<int>& order) const {
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

	std::reverse(order.begin(), order.end());
}

MiniBucketJGHeur::MiniBucketJGHeur(Problem* p, Pseudotree* pt,
		ProgramOptions* po, int ib) :
		Heuristic(p, pt, po), _p(p), _pt(pt), _jglp(), _memlimit(0), _options(po) {

	// preprocessing
	preprocess(NULL);

	_jglp = mex::jglp(copyFactors());  // copyFactors()
	//_jglp = mex::jglp(_mplp.factors());
	_jglp.setProperties("DoMatch=1,DoJG=1");
	_jglp.setIBound(ib);

	// set the var types (excluding the dummy var)
	mex::vector<bool> varTypes(p->getN());
	for (int var = 0; var < p->getN(); ++var) {
		if (p->isMap(var)) varTypes[var] = true;
		else varTypes[var] = false;
	}
	_jglp.setVarTypes(varTypes);

	if (!_pt) {
		return;
	}  // no pseudo tree -> skip remainder

	std::vector<int> dfsOrd;
	findDfsOrder(dfsOrd);
	// copy elimination order information (excluding dummy variable)
	//mex::VarOrder ord(pt->getElimOrder().begin(), --pt->getElimOrder().end());
	mex::VarOrder ord(dfsOrd.begin(), --dfsOrd.end());
	_jglp.setOrder(ord);                   // copy elimination order information

	mex::VarOrder parents(_jglp.gmOrig().nvar());  // copy pseudotree information
	for (size_t i = 0; i < _jglp.gmOrig().nvar(); ++i) {
		int par = _pt->getNode(i)->getParent()->getVar();
		parents[i] = (par == _pt->getRoot()->getVar()) ? -1 : par;
	}
	_jglp.setPseudotree(parents);

}

