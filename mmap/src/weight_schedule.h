/*
 * weight_schedule.h
 *
 *  Created on: Apr 23, 2013
 *      Author: radu
 */

#ifndef IBM_ANYTIME_WEIGHT_SCHEDULE_H_
#define IBM_ANYTIME_WEIGHT_SCHEDULE_H_

#include "base.h"

/* A weighting schedule for anytime best-first search */
class WeightSchedule {
protected:
	int m_scheduleType;
	double m_weight;
	double m_paramK;
	double m_paramD;
	double m_iteration;
	double m_initialWeight;

public:

	// returns the next weight (updates m_weight)
	double getNextWeight();

	// returns the initial weight
	double getInitialWeight();

	// returns the current weight
	double getCurrentWeight();

public:
	 WeightSchedule(int schedType, double w, double delta, double denom);
	 ~WeightSchedule();
};

/* inline declarations */
inline WeightSchedule::WeightSchedule(int schedType, double w, double k, double d) :
		m_scheduleType(schedType), m_weight(w), m_paramK(k),
		m_paramD(d), m_iteration(1), m_initialWeight(w) {

	assert(w >= 1);
	assert(k > 0);
	assert(d > 0);
}

inline WeightSchedule::~WeightSchedule() {

}

inline double WeightSchedule::getInitialWeight() {
	return m_initialWeight;
}

inline double WeightSchedule::getCurrentWeight() {
	return m_weight;
}

inline double WeightSchedule::getNextWeight() {

	switch (m_scheduleType) {
	case WEIGHT_SCHED_SUBTRACT:
		m_weight -= m_paramK;
		if (m_weight < 1.001) {
			m_weight = 1.0;
		}
		break;
	case WEIGHT_SCHED_DIVIDE:
		m_weight /= m_paramK;
		if (m_weight < 1.001) {
			m_weight = 1.0;
		}
		break;
	case WEIGHT_SCHED_INVERSE:
		m_iteration++;
		m_weight = m_initialWeight / m_iteration;
		if (m_weight < 1.001) {
			m_weight = 1.0;
		}
		break;
	case WEIGHT_SCHED_PIECEWISE:
		m_iteration++;
		if (m_weight >= m_paramD) {
			m_weight = m_initialWeight / m_iteration;
		} else {
			m_weight /= m_paramK;
		}
		if (m_weight < 1.001) {
			m_weight = 1.0;
		}
		break;
	case WEIGHT_SCHED_SQRT:
		m_weight = sqrt(m_weight) / m_paramK;
		if (m_weight < 1.001) {
			m_weight = 1.0;
		}
		break;
	default:
		assert(false);
		break;
	}

	return m_weight;
}

#endif /* WEIGHT_SCHEDULE_H_ */
