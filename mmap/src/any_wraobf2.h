/*
 * any_wraobf.h
 *
 *  Created on: Oct 29, 2014
 *      Author: radu
 */

#ifndef ANY_WRAOBF2_H_
#define ANY_WRAOBF2_H_

#include "aobf2.h"
#include "weight_schedule.h"

/* Anytime Repairing Weighted AOBF solver*/
class AnyWRAOBF2 : public AOBF2 {

protected:

	// weight update schedule
	scoped_ptr<WeightSchedule> m_schedule;

	// compare class (for repairing the search space)
	struct CompNodeDepthAsc : public std::binary_function<AobfSearchNode*, AobfSearchNode*, bool> {
		bool operator()(const AobfSearchNode *x, const AobfSearchNode *y) const {
			return (x->getDepth() < y->getDepth()); // deeper nodes first
		}
	};

	// AOstart search with repair
	int aostar(bool verbose = false);

	// repair the search space
	void repair();

public:
	AnyWRAOBF2(ProgramOptions* o);
	virtual ~AnyWRAOBF2();

};

/* inline definitions */
inline AnyWRAOBF2::AnyWRAOBF2(ProgramOptions* o) : AOBF2(o) {
	m_schedule.reset(new WeightSchedule(o->weightSched, o->weight, o->weightParamK, o->weightParamD));
	//m_useLayers = true;
}

inline AnyWRAOBF2::~AnyWRAOBF2() {

}



#endif /* ANY_WRAOBF2_H_ */
