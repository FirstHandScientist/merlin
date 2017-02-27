/*
 * evaluator.h
 *
 *  Created on: Oct 22, 2014
 *      Author: radu
 */

#ifndef MAP_SUMMATION_H_
#define MAP_SUMMATION_H_


#include "solver.h"

/**
 * Computes the induced width of the summation problem
 */

class MapSummation : public Solver {

public:

	// init
	void init();

	// main solve routine
	int solve();

public:

	// Constructor
	MapSummation(ProgramOptions* o);
	virtual ~MapSummation();
};

/* inline definitions */

inline MapSummation::MapSummation(ProgramOptions* o) : Solver(o) {

}

inline MapSummation::~MapSummation() {

}


#endif /* EVALUATOR_H_ */
