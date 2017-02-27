/*
 * base.h
 *
 * This is the base file for the anytime heuristic search library
 *
 *  Created on: Mar 22, 2013
 *      Author: radu
 */

#ifndef IBM_MAP_BASE_H_
#define IBM_MAP_BASE_H_

//#define UAI_COMPETITION // for UAI competitions
#define DETERMINISM_THRESHOLD 40
//#define NO_ASSIGNMENT

// C/C++ kernel
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <memory.h>
#include <sys/types.h>
#include <sys/timeb.h>


// STL kernel
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iterator>
#include <vector>
#include <map>
#include <functional>
#include <algorithm>
#include <deque>
#include <list>
#include <queue>
#include <set>
#include <stack>
#include <exception>
#include <string>
#include <limits>
#include <climits>

#include <boost/scoped_ptr.hpp>
using boost::scoped_ptr;

/* type for storing contexts in binary */
typedef std::vector<int> signature_t;

// Constants
#define UNKNOWN -1
#define NOID -1
#define NONE -1
#define TYPE_BAYES 0
#define TYPE_COST 1
#define TYPE_MARKOV 2
#define SEARCH_SUCCESS 0
#define SEARCH_TIMEOUT 1
#define SEARCH_FAILURE 2
#define SEARCH_OPTIMAL 3
#define SEARCH_OUT_OF_MEMORY 4
#define ASCENDING 1
#define DESCENDING 2
#define SORT_HEURISTIC 0
#define SORT_DEPTH 1
#define SORT_SCORE 2
#define CP_NONE 0
#define CP_UNIT 1
#define CP_FULL 2

/* SEARCH NODE TYPE */
#define NODE_AND 1
#define NODE_OR 2

#define TRUE 1
#define FALSE 0

#define VERSIONINFO "2016.03 STANDALONE MAP (March, 2016)"

#define ELEM_NAN std::numeric_limits<double>::quiet_NaN()
#define ELEM_INF ( std::numeric_limits<double>::infinity() )
#define ELEM_ZERO 0.0
#define ELEM_ONE 1.0
#define ISNAN(x) ( x!=x )

#define OP_TIMES *
#define OP_TIMESEQ *=
#define OP_DIVIDE /
#define OP_DIVIDEEQ /=


//#define ELEM_ENCODE(X) log10( X )
//#define ELEM_DECODE(X) pow(10.0, X )

#define ELEM_ENCODE(X) log( X )
#define ELEM_DECODE(X) exp( X )

#define SCALE_LOG(X) ( X )
//#define SCALE_NORM(X) pow(10.0, X )
#define SCALE_NORM(X) exp( X )


#define DOUBLE_PRECISION 1e-10
#define EPSILON 1e-6

const int SUBPROB_WIDTH_INC = 0; // default
const int SUBPROB_WIDTH_DEC = 1;
const int SUBPROB_HEUR_INC = 2;
const int SUBPROB_HEUR_DEC = 3;
const std::string subprob_order[4]
  = {"width-inc","width-dec","heur-inc","heur-dec"};
const std::string solver_status[5]
  = {"success","timeout","failure","optimal","out-of-memory"};

#define SEARCH_CHECKPOINT_TIME 10 // 10 minutes
#ifndef SEARCH_CHECKPOINT_TIME
#define SEARCH_CHECKPOINT_NODE 10000
#endif

// MAP evaluation
typedef enum {
	EVAL_VE,
	EVAL_AO
} evaluation_t;

// Various hashing functions (for AOBF class)
//#define HASH_DEFAULT
#define HASH_MURMUR

// Algorithms
typedef enum {
	ALGO_GENERATOR,				// no algorithm, prepare MAP problems
	ALGO_EVALUATOR,				// no algorithm, evaluates a MAP assignment
	ALGO_SUMMATION,				// no algorithm, evaluates the induced width of the summation problem
	ALGO_BB,					// BB-MAP with mini-bucket-tree heur
	ALGO_BBBT,					// BBBT on MAP vars and AO on SUM vars
	ALGO_AO,					// AO-MAP
	ALGO_AOBB,					// AOBB-MAP
	ALGO_BRAOBB,				// BRAOBB-MAP
	ALGO_PARK,					// BB with join-tree heur by James Park
	ALGO_YUAN,					// BB with join-tree heur by Chengue Yuan
	ALGO_ASTAR,					// A* search with dyanmic MCTE(i)
	ALGO_AOBF,					// AOBF (AND contexts)
	ALGO_AOBF2,					// AOBF (OR contexts)
	ALGO_RBFAOO,				// RBFAOO
	ALGO_SPRBFAOO,				// SHARED MEMORY PARALLEL RBFAOO
	ALGO_ANY_BBBT,				// Anytime BBBT
	ALGO_ANY_BBMB,				// Anytime OR BnB (over MAP vars) with caching
	ALGO_ANY_WAOBF,				// Anytime Weighted AOBF with Restarts
	ALGO_ANY_WRAOBF,			// Anytime Weighted AOBF with Repairs
	ALGO_ANY_WAOBF2,			// Anytime Weighted AOBF2 with Restarts
	ALGO_ANY_WRAOBF2,			// Anytime Weighted AOBF2 with Repairs
	ALGO_ANY_WRBFAOO,			// Anytime Weighted RBFAOO with Restarts
	ALGO_ANY_ARASTAR,			// Anytime Reparing A* (with dynamic MCTE(i))
	ALGO_ANY_LAOBF,				// Anytime AOBF with DFBnB/DFS probing
	ALGO_ANY_SAOBF,				// Anytime stochastic AOBF
	ALGO_ANY_WRBFS,				// Anytime weighted RBFS (based on [Hansen and Zhou, AIJ2007])
	ALGO_ANY_WRBFS2,			// Anytime weighted RBFS with caching
	ALGO_ANY_AAOBF,				// Anytime AAOBF (alternating depth-first and best-first AO search)
	ALGO_ANY_LDFS,				// Anytime learning depth-first AO search
	ALGO_AFSE,					// Anytime Factor Set Elimination (Maua and de Campos, 2012)
	ALGO_ANY_LDFS_BOUNDED,		// Anytime Learning Depth-First AND/OR Search with bounded summation
	ALGO_ANY_BBMB_BOUNDED		// Anytime OR BnB (over MAP vars) with WMB and bounded summation

} algorithm_t;

// Heuristics
typedef enum {
	HEUR_MBE,					// mini-buckets
	HEUR_WMB_MM,				// mini-buckets with weighted moment matching
	HEUR_WMB_JG,				// mini-buckets with JGLP propagation
	HEUR_WMB_JG2,				// experimental WMB-JG (only for UAI-14 competition)
	HEUR_JOIN_TREE,				// join tree
	HEUR_MINI_CLUSTER_TREE		// mini-cluster-tree
} heuristic_t;

// Variable types
typedef enum {
	VAR_SUM,
	VAR_MAX
} var_t;

// Input file format
typedef enum {
	FORMAT_UAI,
	FORMAT_ERGO
} format_t;

// Domain values values
typedef int val_t;

/* type for storing contexts in binary (RBFAOO, BRAOBB) */
typedef std::vector<val_t> context_t;

// Weight schedules
typedef enum {
	WEIGHT_SCHED_SUBTRACT,
	WEIGHT_SCHED_DIVIDE,
	WEIGHT_SCHED_INVERSE,
	WEIGHT_SCHED_PIECEWISE,
	WEIGHT_SCHED_SQRT
} weight_schedule_t;

const std::string schedule_names[5]
  = {"subtract","divide","inverse","piecewise","sqrt"};

// Probing policy
typedef enum {
	PROBE_DFS,		// plain depth-first search probing (with caching, no pruning)
	PROBE_BNB		// branch and bound probing (with caching, with pruning)
} probing_t;

const std::string probing_names[5]
  = {"dfs","bnb"};

/*
 * Minor extension of STL stack class, used by AOBB for the 'stack tree'
 */
class AobbSearchNode;
class MyStack: public std::stack<AobbSearchNode*> {
protected:
	size_t m_children;
	MyStack* m_parent;
public:
	size_t getChildren() const {
		return m_children;
	}
	void addChild() {
		m_children += 1;
	}
	void delChild() {
		m_children -= 1;
	}
	MyStack* getParent() const {
		return m_parent;
	}
public:
	MyStack(MyStack* p) :
			m_children(0), m_parent(p) {
	}
};



/*//////////////////////////////////////////////////////////////*/
/* Boost random number library */
#include <boost/random/linear_congruential.hpp>
/* Boost lexical cast (for version string) */
#include <boost/lexical_cast.hpp>

/* static random number generator */
class rand {
private:
	static int state;
	static boost::minstd_rand _r;
public:

	static void seed(const int& s) {
		_r.seed(s);
	}
	static int seed() {
		return state;
	}
	static int max() {
		return _r.max();
	}
	static int next() {
		return state = _r();
	}
	static int next(const int& hi) {
		return static_cast<int>((state = _r()) / (_r.max() + 1.0) * hi);
	}
	static double probability() {
		return ( static_cast<double>(state = _r()) / static_cast<double> (_r.max() ) );
	}
};

/*//////////////////////////////*/

size_t strHasher(const std::string& str);

// UAI competition only
int uai(int argc, char** argv);

#endif /* BASE_H_ */
