/*
 * any_waobf.h
 *
 *  Created on: Oct 29, 2014
 *      Author: radu
 */

#ifndef ANY_WAOBF2_H_
#define ANY_WAOBF2_H_


#include "aobf2.h"
#include "weight_schedule.h"

class AnyWAOBF2 : public AOBF2 {

protected:

	// weight update schedule
	scoped_ptr<WeightSchedule> m_schedule;

public:

	// solver
	virtual int solve();

public:
	AnyWAOBF2(ProgramOptions* o);
	~AnyWAOBF2();
};

/* inline functions */
inline AnyWAOBF2::AnyWAOBF2(ProgramOptions* o) : AOBF2(o) {
	m_schedule.reset(new WeightSchedule(o->weightSched, o->weight, o->weightParamK, o->weightParamD));
}

inline AnyWAOBF2::~AnyWAOBF2() {

}

#endif /* ANY_WAOBF2_H_ */
