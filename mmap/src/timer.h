// timer.h -- A timer

/*
 * Author: Radu Marinescu
 *
 * Copyright (c) IBM Research, 2012
 *
 * This program is under IBM license. Do not distribute.
 *
 */

/*
 * NOTE: This is an internal header file.
 * You should not attempt to use it directly.
 */

#ifndef IBM_ANYTIME_TIMER_H_
#define IBM_ANYTIME_TIMER_H_


#include "base.h"

#include <unistd.h>
#include <sys/time.h>
#include <sys/times.h>

class Timer {

protected:

	// start time
	double m_start;

	// stop time
	double m_stop;

	// time limit (default NONE(-1), no limit)
	double m_limit;

	// time elapsed in seconds
	double m_elapsed;

	// timeout flag
	bool m_timeout;

	// running flag
	bool m_running;

protected:

	// get the CPU time in seconds
	double cpuTime();

public:

	// start the timer
	void start();

	// stop the timer
	void stop();

	// reset the timer
	void reset(double limit = NONE);

	// check for timeout
	bool timeout();

	// get the elapsed time from start
	double elapsed();

public:

	// default constructor
	Timer(double limit = NONE);

	// destructor
	~Timer();

private:
	Timer(const Timer&);
};

/* inline definitions */

// constructor
inline Timer::Timer(double limit) :
		m_start(0), m_stop(0), m_limit(limit), m_elapsed(0),
		m_timeout(false), m_running(false) {
	/* empty */
}

// destructor
inline Timer::~Timer() {

}

// get the CPU time (in seconds)
inline double Timer::cpuTime() {

	static struct tms buf;

	times(&buf);
	return ((double) (buf.tms_utime+buf.tms_stime+buf.tms_cutime+buf.tms_cstime)) / ((double) sysconf(_SC_CLK_TCK));
}

// reset the timer
inline void Timer::reset(double limit) {
	m_start = m_stop = m_elapsed = 0;
	m_limit = limit;
	m_timeout = m_running = false;
}

// start the timer
inline void Timer::start() {
	if (!m_running) {
		m_start = cpuTime();
		m_running = true;
	}
}

// stop the timer
inline void Timer::stop() {
	if (m_running) {
		m_stop = cpuTime();
		m_elapsed = (m_stop - m_start);
		m_running = false;
	}
}

// check for timeout
inline bool Timer::timeout() {
	if (m_limit < 0) {
		return false;
	}

	if (!m_timeout) {
		double e = (cpuTime() - m_start);
		if (e >= m_limit) {
			m_timeout = true;
		}
	}

	return m_timeout;
}

// get the elapsed time
inline double Timer::elapsed() {
	if (m_running) {
		m_elapsed = cpuTime() - m_start;
	}

	return m_elapsed;
}

#endif /* TIMER_H_ */
