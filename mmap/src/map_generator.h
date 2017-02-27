/*
 * map_generator.h
 *
 *  Created on: Sep 27, 2013
 *      Author: radu
 */

#ifndef IBM_MAP_MAP_GENERATOR_H_
#define IBM_MAP_MAP_GENERATOR_H_

#include "solver.h"

class MapGenerator : public Solver {

public:

	// initialize
	void init();

	// Main solver routine
	int solve();

public:

	MapGenerator(ProgramOptions *o);
	virtual ~MapGenerator();

};

/* inline definitions */
inline MapGenerator::MapGenerator(ProgramOptions* o) : Solver(o) {

}

inline MapGenerator::~MapGenerator() {

}


#endif /* MAP_GENERATOR_H_ */
