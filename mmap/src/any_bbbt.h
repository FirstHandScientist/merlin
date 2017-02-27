/*
 * any_bbbt.h
 *
 *  Created on: Oct 16, 2014
 *      Author: radu
 */

#ifndef ANY_BBBT_H_
#define ANY_BBBT_H_


#include "search.h"

#include "aobb_search_node.h"
#include "aobb_search_space.h"
#include "bound_propagator.h"


/**
 * Anytime BBBT search over MAP variables
 * - static vo and incremental mbte heuristic (Yuan like)
 * - dynamic vo and full mbte heuristic (Park like)
 * - conditioned SUM subproblem is approximated by a lower bound on its
 *   corresponding partition function
 *     -- WMB scheme
 *     -- negative tree reweighted BP
 *     -- ABP
 *     -- probabilistic lower bounds (Markov inequality)
 */

class AnyBBBT : public Search {

protected:

	// number of expansions (over MAP variables)
	size_t m_expansions;

	// number of cache hits (over SUM variables)
	size_t m_cacheHits;
	size_t m_deadendsCP;
	size_t m_deadendsUB;

	// flag indicating static or dynamic VO (default true)
	bool m_useStaticVo;

	// last MAP variable in the chain pseudo tree (roots the SUM subproblem)
	int m_lastMapVar;

	std::list<int> m_mapVars;

	std::vector<int> m_assumptions;

	// Conditioned SUM problem order
	mex::VarOrder _order;
	size_t _inducedWidth;
	bool _saved;

#ifdef UAI_COMPETITION
	void writeSolution(std::vector<int>& assignment);
#endif

private:

	// lower bound the partition function of the conditioned SUM subproblem
	double _lowerBound(const std::vector<int>& assignment, double& exact);

//	// AO search
//	std::stack<AobbSearchNode*> ao_stack;
//	std::vector<AobbSearchNode*> ao_expand;
//	scoped_ptr<AobbSearchSpace> ao_space;
//	scoped_ptr<BoundPropagator> ao_propagator;
//
//	void _initAoSearch();
//	void _heuristicOR(AobbSearchNode*);
//	bool _doProcess(AobbSearchNode*);
//	void _addCacheContext(AobbSearchNode*, const set<int>&) const;
//	bool _doCaching(AobbSearchNode* node);
//	AobbSearchNode* _nextLeaf();
//	bool _generateChildrenAND(AobbSearchNode*, vector<AobbSearchNode*>&);
//	bool _generateChildrenOR(AobbSearchNode* n, vector<AobbSearchNode*>& chi);
//	AobbSearchNode* _nextNode();
//	bool _doExpand(AobbSearchNode* n);
//	bool _doPruning(AobbSearchNode* n);
//	bool _canBePruned(AobbSearchNode* n);
//
//	// Solve the conditioned SUM subproblem by AO search
//	double _ao(const std::vector<int>& assignment);

protected:

	// get the path assignment (as SAT variables)
//	void getPathAssignment(AobbSearchNode* node, std::vector<int>& vars);

	// get the value of a full MAP assignment
	double getMapAssignmentValue(const std::vector<int>& assignment);

	// depth first recursive branch-and-bound (static variable ordering)
	void bbStatic(std::list<int> M, std::vector<int> assignment);

	// depth first recursive branch-and-bound (dynamic variable ordering)
	void bbDynamic(std::set<int> M, std::vector<int> assignment, std::set<int> F);

	// select next unassigned MAP variable (dynamic variable ordering)
	int selectNextVar(const std::set<int>& M,
			std::vector<std::pair<int, double> >& values,
			const std::vector<int>& assignment);

	// custom compare function
	struct heurGreater {
		bool operator()(const std::pair<int, double>& a,
				const std::pair<int, double>& b) const {
			return a.second > b.second;
		}
	};

	// custom compare function
	struct heurLess {
		bool operator()(const std::pair<int, double>& a,
				const std::pair<int, double>& b) const {
			return a.second < b.second;
		}
	};

public:

	// initialize
	void init();

	// main solve routine
	int solve();

public:

	// Constructor
	AnyBBBT(ProgramOptions* o);
	virtual ~AnyBBBT();

protected:

	// Time and node checkpoints
	double m_prevTimeCheckpoint;
	double m_prevNodeCheckpoint;
};

/* inline definitions */

inline AnyBBBT::AnyBBBT(ProgramOptions* o) : Search(o) {
	m_expansions = 0;
	m_cacheHits = 0;
	m_deadendsCP = 0;
	m_deadendsUB = 0;
	m_useStaticVo = true;
	m_lastMapVar = NONE;
	m_prevNodeCheckpoint = 0;
	m_prevTimeCheckpoint = 0;
	_saved = false;
	_inducedWidth = 0;
}

inline AnyBBBT::~AnyBBBT() {

}

#endif /* ANY_BBBT_H_ */
