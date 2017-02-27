/*
 * CacheTable.h
 *
 *  Copyright (C) 2008-2012 Lars Otten
 *  This file is part of DAOOPT.
 *
 *  DAOOPT is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  DAOOPT is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with DAOOPT.  If not, see <http://www.gnu.org/licenses/>.
 *  
 *  Created on: Dec 4, 2008
 *      Author: Lars Otten <lotten@ics.uci.edu>
 */

#ifndef IBM_ANYTIME_BRAOBB_CACHETABLE_H_
#define IBM_ANYTIME_BRAOBB_CACHETABLE_H_

#include "base.h"
#include "utils.h"

#include <vector>
#include <string>
#include <malloc.h>

#include "boost/unordered_map.hpp"
#include "boost/unordered_set.hpp"

//#define hash_set boost::unordered_set
//#define hash_map boost::unordered_map

//typedef boost::unordered_map<context_t, double> context_hash_map;

#ifndef NO_ASSIGNMENT
typedef boost::unordered_map<context_t, pair<double, vector<val_t> > > context_hash_map;
#else
typedef boost::unordered_map<context_t, double> context_hash_map;
#endif

class AobbCacheTable {

private:
	bool m_full;
	int m_size;
	std::vector<context_hash_map*> m_tables;

public:

#ifndef NO_ASSIGNMENT
	virtual void write(int n, const context_t& ctxt, double v, vector<val_t>& sol) throw (int);
	virtual pair<double, vector<val_t> > read(int n, const context_t& ctxt) const throw (int);
#else
	virtual void write(int n, const context_t& ctxt, double v) throw (int);
	virtual double read(int n, const context_t& ctxt) const throw (int);
#endif

//	virtual void write(int n, const context_t& ctxt, double v)
//			throw (int);
//	virtual double read(int n, const context_t& ctxt) const
//			throw (int);

	virtual void reset(int n);

private:
	int memused() const;

public:
	virtual void printStats() const;

public:
	AobbCacheTable(int size);
	virtual ~AobbCacheTable();
};

class UnCacheTable: public AobbCacheTable {
public:

//	void write(int n, const context_t& ctxt, double v) throw (int) {
//		throw UNKNOWN;
//	}
//	double read(int n, const context_t& ctxt) const throw (int) {
//		throw UNKNOWN;
//	}
#ifndef NO_ASSIGNMENT
	void write(int n, const context_t& ctxt, double v, const vector<val_t>& sol) throw (int) {
		throw UNKNOWN;
	}
	pair<double, vector<val_t> > read(int n, const context_t& ctxt) const throw (int) {
		throw UNKNOWN;
	}
#else
	void write(int n, const context_t& ctxt, double v) throw (int) {
		throw UNKNOWN;
	}
	double read(int n, const context_t& ctxt) const throw (int) {
		throw UNKNOWN;
	}
#endif

	void reset(int n) {
	}

public:
	void printStats() const {
	}

public:
	UnCacheTable() :
			AobbCacheTable(0) {
	}
	~UnCacheTable() {
	}

};

/* Inline definitions */

/* inserts a value into the respective cache table, throws an int if
 * insert non successful (memory limit or index out of bounds) */
#ifndef NO_ASSIGNMENT
inline void AobbCacheTable::write(int n, const context_t& ctxt,
		double v, vector<val_t>& sol) throw (int) {
#else
inline void AobbCacheTable::write(int n, const context_t& ctxt,
		double v) throw (int) {
#endif
	assert(n < m_size);

	// create hash table if needed
	if (!m_tables[n]) {
		m_tables[n] = new context_hash_map;
	}

	// this will write only if entry not present yet
	//m_tables[n]->insert( make_pair(ctxt, v) );
#ifndef NO_ASSIGNMENT
	  m_tables[n]->insert( make_pair(ctxt, make_pair(v,sol) ) );
#else
	  m_tables[n]->insert( make_pair(ctxt, v) );
#endif
}

/* tries to read a value from a table, throws an int (UNKNOWN) if not found */
#ifndef NO_ASSIGNMENT
inline pair<double, vector<val_t> > AobbCacheTable::read(int n, const context_t& ctxt) const throw (int) {
#else
inline double AobbCacheTable::read(int n, const context_t& ctxt) const throw (int) {
#endif

	assert(n < m_size);
	// does cache table exist?
	if (!m_tables[n])
		throw UNKNOWN;
	// look for actual entry
	context_hash_map::const_iterator it = m_tables[n]->find(ctxt);
	if (it == m_tables[n]->end())
		throw UNKNOWN;
	return it->second;
}

inline void AobbCacheTable::reset(int n) {
	assert(n<m_size);
	if (m_tables[n]) {
//    m_tables[n]->clear();
		delete m_tables[n];
		m_tables[n] = NULL;
		m_full = false;
//    cout << "Reset cache table " << n << endl;
	}
}

inline AobbCacheTable::AobbCacheTable(int size) :
		m_full(false), m_size(size), m_tables(size) {
	for (int i = 0; i < size; ++i)
		m_tables[i] = NULL;
}

inline void AobbCacheTable::printStats() const {
	ostringstream ss;
	ss << "Cache statistics:";
	for (vector<context_hash_map *>::const_iterator it = m_tables.begin();
			it != m_tables.end(); ++it) {
		if (*it)
			ss << " " << (*it)->size();
		else
			ss << " .";
	}
	ss << endl;

	std::cout << (ss.str());
}

inline AobbCacheTable::~AobbCacheTable() {
	for (vector<context_hash_map *>::iterator it = m_tables.begin();
			it != m_tables.end(); ++it)
		if (*it)
			delete *it;
}

inline int AobbCacheTable::memused() const {

	// read mem stats
	struct mallinfo info;
	info = mallinfo();

	int m = info.hblkhd + info.uordblks;
	m /= 1024 * 1024; // to MByte

	return m;
}

#endif /* CACHETABLE_H_ */
