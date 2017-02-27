/*
 * any_waobf.h
 *
 *  Created on: Oct 29, 2014
 *      Author: radu
 */

#ifndef ANY_WAOBF_H_
#define ANY_WAOBF_H_


#include "aobf.h"
#include "weight_schedule.h"

class AnyWAOBF : public AOBF {

protected:

	// weight update schedule
	scoped_ptr<WeightSchedule> m_schedule;

public:

	// solver
	virtual int solve();

public:
	AnyWAOBF(ProgramOptions* o);
	~AnyWAOBF();
};

/* inline functions */
inline AnyWAOBF::AnyWAOBF(ProgramOptions* o) : AOBF(o) {
	m_schedule.reset(new WeightSchedule(o->weightSched, o->weight, o->weightParamK, o->weightParamD));
}

inline AnyWAOBF::~AnyWAOBF() {

}

#endif /* ANY_WAOBF_H_ */
