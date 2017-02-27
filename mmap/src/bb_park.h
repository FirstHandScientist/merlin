/*
 * bb_park.h
 *
 *  Created on: Sep 26, 2013
 *      Author: radu
 */

#ifndef IBM_MAP_BB_PARK_H_
#define IBM_MAP_BB_PARK_H_


#include "search.h"


// OR Branch-and-Bound with join/bucket-tree heuristics [Park and Darwiche, 2003]
class BBPark : public Search {

protected:

	// number of expansions
	size_t m_expansions;

protected:

	// depth first recursive branch-and-bound
	void bb(std::set<int> M, std::vector<int> assignment);

	// select next unassigned MAP variable
	int selectNextVar(const set<int>& M,
			std::vector<std::pair<int, double> >& values,
			const std::vector<int>& assignment);

	// custom compare function
	struct descHeurComp {
		bool operator()(const std::pair<int, double>& a,
				const std::pair<int, double>& b) const {
			return a.second > b.second;
		}
	};

public:

	// initialize
	void init();

	// main solver routine
	int solve();

public:

	// Constructor
	BBPark(ProgramOptions* opt);
	virtual ~BBPark();

protected:

	// Time and node checkpoints
	double m_prevTimeCheckpoint;
	double m_prevNodeCheckpoint;
};


/* inline definitions */

inline BBPark::BBPark(ProgramOptions* opt) : Search(opt) {
	m_expansions = 0;
	m_prevTimeCheckpoint = 0;
	m_prevNodeCheckpoint = 0;
}

inline BBPark::~BBPark() {

}



#endif /* BB_PARK_H_ */
