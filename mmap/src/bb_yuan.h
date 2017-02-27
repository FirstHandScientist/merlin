/*
 * bb_yaun.h
 *
 *  Created on: Sep 30, 2013
 *      Author: radu
 */

#ifndef IBM_MAP_BB_YAUN_H_
#define IBM_MAP_BB_YAUN_H_

#include "search.h"

/**
 * OR Branch-and-Bound with Incremental Join/Bucket-tree Bounds
 *
 */

class BBYuan : public Search {

protected:

	// number of expansions
	size_t m_expansions;

protected:

	// depth first recursive branch-and-bound (static variable ordering)
	void bb(std::list<int> M, std::vector<int> assignment);

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
	BBYuan(ProgramOptions* o);
	virtual ~BBYuan();

protected:

	// Time and node checkpoints
	double m_prevTimeCheckpoint;
	double m_prevNodeCheckpoint;

};

/* inline definitions */

inline BBYuan::BBYuan(ProgramOptions* o) : Search(o) {
	m_expansions = 0;
	m_prevNodeCheckpoint = 0;
	m_prevTimeCheckpoint = 0;
}

inline BBYuan::~BBYuan() {

}

#endif /* BB_YAUN_H_ */
