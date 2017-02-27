/*
 * any_wrbfaoo.h
 *
 *  Created on: Oct 29, 2014
 *      Author: radu
 */

#ifndef ANY_WRBFAOO_H_
#define ANY_WRBFAOO_H_


#include "rbfaoo.h"
#include "weight_schedule.h"

/**
 * Anytime Restarting Weighted Recursive Best-First AND/OR Search (aRBFAOO-W)
 */

class AnyWRBFAOO : public RBFAOO {

protected:

	// weight update schedule
	scoped_ptr<WeightSchedule> m_schedule;

public:

	int solve();

public:

	// Constructor
	AnyWRBFAOO(ProgramOptions* opt);

	// Destructor
	~AnyWRBFAOO();

};

/* inline definitions */
inline AnyWRBFAOO::AnyWRBFAOO(ProgramOptions* o) : RBFAOO(o) {
	m_schedule.reset(new WeightSchedule(o->weightSched, o->weight, o->weightParamK, o->weightParamD));
}

inline AnyWRBFAOO::~AnyWRBFAOO() {

}



#endif /* ANY_WRBFAOO_H_ */
