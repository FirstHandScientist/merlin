/*
 * any_hbfs.h
 *
 *  Created on: 23 Sep 2015
 *      Author: radu
 */

#ifndef SRC_ANY_WRBFS_H_
#define SRC_ANY_WRBFS_H_

#include "search.h"

#include "aobb_search_node.h"
#include "aobb_search_space.h"
#include "bound_propagator.h"

/**
 * Implementation of Hansen's anytime weighted RBFS search for Marginal MAP.
 * It assumes an OR search space over the MAP variables, and it's linear space
 * (the summation subproblems are solved exactly by AND/OR search with caching).
 *
 */

class AnyWRBFS : public Search {

	class SearchNode {
	protected:
		int m_var;
		int m_val;
		double m_h;
		double m_g;
		double m_f;
		double m_F;
		bool m_goal;
		double m_cost;
	public:

		SearchNode(int var = UNKNOWN, int val = UNKNOWN) : m_var(var), m_val(val), m_h(ELEM_NAN),
			m_g(ELEM_NAN), m_f(ELEM_NAN), m_F(ELEM_NAN), m_goal(false), m_cost(ELEM_NAN) {};
		~SearchNode() {};
		SearchNode(const SearchNode& n) {
			m_var = n.m_var;
			m_val = n.m_val;
			m_h = n.m_h;
			m_g = n.m_g;
			m_f = n.m_f;
			m_F = n.m_F;
			m_goal = n.m_goal;
			m_cost = n.m_cost;
		};
		inline const int& var() const { return m_var; }
		inline const int& val() const { return m_val; }
		inline double& h() { return m_h; };
		inline double& g() { return m_g; };
		inline double& f() { return m_f; };
		inline double& F() { return m_F; };
		inline bool& goal() { return m_goal; };
		inline double& cost() { return m_cost; };
		inline const double& h() const { return m_h; };
		inline const double& g() const { return m_g; };
		inline const double& f() const { return m_f; };
		inline const double& F() const { return m_F; };
		inline const bool& goal() const { return m_goal; };
		inline const double& cost() const { return m_cost; };
		friend ostream& operator <<(ostream& out, const SearchNode& n) {
			out << "node x" << n.var() << "=" << n.val() << " : g=" << n.g()
					<< " h=" << n.h() << " f=" << n.f() << " F=" << n.F()
					<< " goal " << n.goal();
			return out;
		}
	};

protected:

	// Overestimation
	double m_overestimation;

	// Node expansions internal counters
	size_t m_numExpanded;
	size_t m_numBacktracks;

	// Checkpoints
	double m_prevTimeCheckpoint;
	double m_prevNodeCheckpoint;

	// Weight (default 1.0) to inflate the heuristic
	double m_weight;

	// Number of SUM evaluations
	size_t m_numSumEvals;

	// Global upper bound (best solution found so far)
	double m_upperBound;
	double m_lowerBound;

	// Last MAP variable in the order
	int m_lastMapVar;

protected:

	double eval(SearchNode& n);
	double* heuristic(int var);
	std::vector<SearchNode> expand(SearchNode& n);
	double rbfs(SearchNode& n, double f, double threshold, int depth);
	size_t getBestFValueIndex(std::vector<SearchNode>& nodes);
	size_t getNextBestFValueIndex(std::vector<SearchNode>& nodes, size_t bestIndex);

public:

	// initialize
	void init();

	// solver
	virtual int solve();

protected:

	// AO search (for the conditioned SUM problem)
	std::stack<AobbSearchNode*> ao_stack;
	std::vector<AobbSearchNode*> ao_expand;
	scoped_ptr<AobbSearchSpace> ao_space;
	scoped_ptr<BoundPropagator> ao_propagator;

	// Solve the conditioned SUM subproblem by AO search
	double _ao(int var, const std::vector<int>& assignment);

	void _initAoSearch(int var);
	void _heuristicOR(AobbSearchNode*);
	bool _doProcess(AobbSearchNode*);
	void _addCacheContext(AobbSearchNode*, const set<int>&) const;
	bool _doCaching(AobbSearchNode* node);
	AobbSearchNode* _nextLeaf();
	bool _generateChildrenAND(AobbSearchNode*, vector<AobbSearchNode*>&);
	bool _generateChildrenOR(AobbSearchNode* n, vector<AobbSearchNode*>& chi);
	AobbSearchNode* _nextNode();
	bool _doExpand(AobbSearchNode* n);

public:

	AnyWRBFS(ProgramOptions* o);
	virtual ~AnyWRBFS() {};

private:
	struct fnAsc {
		bool operator() (const std::pair<int, double>& a,
				const std::pair<int, double>& b) {
			return (a.second < b.second);
		}
	};
};


#endif /* SRC_ANY_WRBFS_H_ */
