/*
 * evaluator.h
 *
 *  Created on: Oct 22, 2014
 *      Author: radu
 */

#ifndef EVALUATOR_H_
#define EVALUATOR_H_


#include "ao.h"

/**
 * Evaluates a MAP assignment by solving the conditioned SUM problem
 * - by Variable Elimination (BE)
 * - by AND/OR search with adaptive caching
 */

class MapEvaluator : public AO {

protected:

	evaluation_t m_eval;

public:

	// main solve routine
	int solve();

public:

	// Constructor
	MapEvaluator(ProgramOptions* o);
	virtual ~MapEvaluator();
};

/* inline definitions */

inline MapEvaluator::MapEvaluator(ProgramOptions* o) : AO(o) {
	m_eval = EVAL_VE; // default
}

inline MapEvaluator::~MapEvaluator() {

}


#endif /* EVALUATOR_H_ */
