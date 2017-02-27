/*
 * any_arastar.h
 *
 *  Created on: Oct 29, 2014
 *      Author: radu
 */

#ifndef ANY_ARASTAR_H_
#define ANY_ARASTAR_H_


#include "astar.h"
#include "weight_schedule.h"

class AnyARAstar : public Astar {

protected:

	// weight update schedule
	scoped_ptr<WeightSchedule> m_schedule;

	void astar();

	void repair();

public:

	AnyARAstar(ProgramOptions* o);
	~AnyARAstar();
};

/* inline functions */
inline AnyARAstar::AnyARAstar(ProgramOptions* o) : Astar(o) {
	m_schedule.reset(new WeightSchedule(o->weightSched, o->weight, o->weightParamK, o->weightParamD));
}

inline AnyARAstar::~AnyARAstar() {

}


#endif /* ANY_ARASTAR_H_ */
