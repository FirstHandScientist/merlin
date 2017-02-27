/*
 * astar.h
 *
 *  Created on: Mar 31, 2014
 *      Author: radu
 */

#ifndef IBM_MAP_ASTAR_H_
#define IBM_MAP_ASTAR_H_

#include "search.h"

#include "aobb_search_node.h"
#include "aobb_search_space.h"
#include "bound_propagator.h"

/*
 * Standard A* Search solver (Astar)
 *
 * - best-first search over the MAP variables
 * - conditioned SUM subproblem solved by AOC
 * - dynamic variable ordering (for MAP variables)
 * - dynamic heuristic computation (MCTE(i))
 *
 */

class Astar : public Search {

protected:

	// a search node
	class AstarSearchNode {
	protected:
		int m_var;
		int m_val;
		double m_f;
		double m_g;
		double m_h;
		bool m_terminal;
		int m_depth;
		list<pair<int, int> > m_assignment; // partial MAP assignment

	public:
		inline int var() const {
			return m_var;
		}
		inline int& var() {
			return m_var;
		}
		inline int val() const {
			return m_val;
		}
		inline int& val() {
			return m_val;
		}
		inline double f() const {
			return m_f;
		}
		inline double& f() {
			return m_f;
		}
		inline double g() const {
			return m_g;
		}
		inline double& g() {
			return m_g;
		}
		inline double h() const {
			return m_h;
		}
		inline double& h() {
			return m_h;
		}
		inline bool terminal() const {
			return m_terminal;
		}
		inline bool& terminal() {
			return m_terminal;
		}
		inline int depth() const {
			return m_depth;
		}
		inline int& depth() {
			return m_depth;
		}
		inline const list<pair<int, int> >& assignment() const {
			return m_assignment;
		}
		inline list<pair<int, int> >& assignment() {
			return m_assignment;
		}
		// comparison operator (needed for std::less<T>)
		inline friend bool operator<(const AstarSearchNode& a, const AstarSearchNode& n) {
			return a.m_f < n.m_f;
		}

		friend std::ostream& operator<<(std::ostream& out, const AstarSearchNode& n) {
			out << "node  (x" << n.m_var << "=" << n.m_val << ") f=" << n.m_f;
			out << " [ ";
			for (list<pair<int, int> >::const_iterator li = n.m_assignment.begin();
					li != n.m_assignment.end(); ++li) {
				out << "x" << li->first << "=" << li->second << ", ";
			}
			out << "]";
			return out;
		}
		;

	public:

		// copy constructor
		AstarSearchNode(const AstarSearchNode& n) {
			m_var = n.m_var;
			m_val = n.m_val;
			m_f = n.m_f;
			m_g = n.m_g;
			m_h = n.m_h;
			m_terminal = n.m_terminal;
			m_depth = n.m_depth;
			m_assignment = n.m_assignment;
		}

		// assignment operator =
		AstarSearchNode& operator=(const AstarSearchNode& n) {
			m_var = n.m_var;
			m_val = n.m_val;
			m_f = n.m_f;
			m_g = n.m_g;
			m_h = n.m_h;
			m_terminal = n.m_terminal;
			m_depth = n.m_depth;
			m_assignment = n.m_assignment;
			return *this;
		}

		// constructor
		AstarSearchNode(int var, int val) : m_var(var), m_val(val),
			m_f(ELEM_NAN), m_g(ELEM_NAN), m_h(ELEM_NAN), m_terminal(false), m_depth(0) {};

		// destructor
		~AstarSearchNode() {};
	};

	class Heap {
	protected:

		// the elements
		std::vector<AstarSearchNode> m_heap;

		// compare class (for the OPEN list)
		class comp_t {
		public:
			bool operator()(const AstarSearchNode& x, const AstarSearchNode& y) const {
				double fx = x.f(), fy = y.f();
				if (fx > fy) {
					return true;
				} else if (fx == fy) {
					return (x.depth() < y.depth());
				} else {
					return false;
				}
			}
		};

	public:

		// add an element to the heap
		void push(AstarSearchNode& n) {
			m_heap.push_back(n);
			std::push_heap(m_heap.begin(), m_heap.end(), comp_t());
		}

		// get the highest priority node
		AstarSearchNode top() {
			assert( m_heap.empty() == false);
			return m_heap.front();
		}

		// remove the highest priority node
		void pop() {
			pop_heap(m_heap.begin(), m_heap.end(), comp_t());
			m_heap.pop_back();
		}

		// heapify the container
		void heapify() {
			std::make_heap(m_heap.begin(), m_heap.end(), comp_t());
		}

		// get the size of the heap
		inline size_t size() {
			return m_heap.size();
		}

		// check if empty
		inline bool empty() {
			return m_heap.empty();
		}

		// get the elements
		inline std::vector<AstarSearchNode>& elements() {
			return m_heap;
		}

	public:
		// constructor
		Heap() {};

		// destructor
		~Heap() {};
	};

protected:

	// number of expansions (over MAP variables)
	size_t m_expansions;

	// list of MAP variables
	std::set<int> m_mapVars;

	// the OPEN list
	//std::priority_queue<AstarSearchNode> m_open;
	Heap m_open;

	// last MAP variable in the pseudo tree
	int m_lastMapVar;
	int m_firstMapVar;

	// weight
	double m_epsilon;

	// number of SUM evaluations
	size_t m_sumEvals;

protected:

	// AO search (for the conditioned SUM problem)
	std::stack<AobbSearchNode*> ao_stack;
	std::vector<AobbSearchNode*> ao_expand;
	scoped_ptr<AobbSearchSpace> ao_space;
	scoped_ptr<BoundPropagator> ao_propagator;

	// Solve the conditioned SUM subproblem by AO search
	double _ao(const std::vector<int>& assignment);

	void _initAoSearch();
	void _heuristicOR(AobbSearchNode*);
	bool _doProcess(AobbSearchNode*);
	void _addCacheContext(AobbSearchNode*, const set<int>&) const;
	bool _doCaching(AobbSearchNode* node);
	AobbSearchNode* _nextLeaf();
	bool _generateChildrenAND(AobbSearchNode*, vector<AobbSearchNode*>&);
	bool _generateChildrenOR(AobbSearchNode* n, vector<AobbSearchNode*>& chi);
	AobbSearchNode* _nextNode();
	bool _doExpand(AobbSearchNode* n);
	bool _doPruning(AobbSearchNode* n);
	bool _canBePruned(AobbSearchNode* n);

protected:

	// get the value of a full MAP assignment
	double getMapAssignmentValue(const std::vector<int>& assignment);

	// get the value of a partial MAP assignment
	double getPartialMapAssignmentValue(const std::vector<int>& assignment);

	// select next unassigned MAP variable (static variable ordering)
	int selectNextVar(int currVar,
			std::vector<std::pair<int, double> >& values,
			const std::vector<int>& assignment);

	// select next unassigned MAP variable (dynamic variable ordering)
	int selectNextVar(const set<int>& M,
			std::vector<std::pair<int, double> >& values,
			const std::vector<int>& assignment);

	// A* search over the MAP variables
	virtual void astar();
	void expand(const AstarSearchNode& n);
	void getMapAssignment(const AstarSearchNode& n,
			std::vector<int>& assignment);
public:

	// initialize
	void init();

	// main solve routine
	int solve();

public:

	// Constructor
	Astar(ProgramOptions* o);
	virtual ~Astar();

protected:

	// Time and node checkpoints
	double m_prevTimeCheckpoint;
	double m_prevNodeCheckpoint;

};


/* inline definitions */

inline Astar::Astar(ProgramOptions* o) : Search(o) {
	m_expansions = 0;
	m_prevNodeCheckpoint = 0;
	m_prevTimeCheckpoint = 0;
	ao_expand.reserve(128);
	m_lastMapVar = UNKNOWN;
	m_firstMapVar = UNKNOWN;
	m_epsilon = 1.0; //o->weight;
	m_sumEvals = 0;
}

inline Astar::~Astar() {

}


#endif /* ASTAR_H_ */
