/*
 * ve.h
 *
 *  Created on: Apr 17, 2013
 *      Author: radu
 */

#ifndef IBM_ANYTIME_VE_H_
#define IBM_ANYTIME_VE_H_


#include "base.h"
#include "solver.h"

/* Variable Elimination solver */

class VariableElimination : public Solver {

public:
	VariableElimination(ProgramOptions *o);
	virtual ~VariableElimination();
};

#endif /* VE_H_ */
