/*
 * SearchSpace.h
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
 *  Created on: Nov 4, 2008
 *      Author: Lars Otten <lotten@ics.uci.edu>
 */

#ifndef IBM_MAP_BRAOBB_SEARCHSPACE_H_
#define IBM_MAP_BRAOBB_SEARCHSPACE_H_

#include "aobb_search_node.h"
#include "aobb_cache_table.h"

/* main search space structure for worker nodes */
struct AobbSearchSpace {

	size_t nodesOR;         // total number of OR nodes expanded
	size_t nodesAND;        // total number of AND nodes expanded
	size_t nodesORmap;      // number of OR nodes expanded (MAP only)
	size_t nodesANDmap;     // number of AND nodes expanded	(MAP only)
	size_t numSumEvals;
	AobbSearchNode* root; // true root of the search space, always a dummy OR node
	AobbCacheTable* cache; // Cache table

	AobbSearchNode* getTrueRoot() const;
	AobbSearchSpace();
	~AobbSearchSpace();
};

inline AobbSearchSpace::AobbSearchSpace() :
		nodesOR(0), nodesAND(0), nodesORmap(0), nodesANDmap(0), numSumEvals(0),
		root(NULL), cache(NULL) {
	/* intentionally empty at this point */
}

inline AobbSearchSpace::~AobbSearchSpace() {
	if (root) {
		delete root;
	}
	if (cache) {
		delete cache;
	}
}

/* returns the relevant root node in both conditioned and unconditioned cases */
inline AobbSearchNode* AobbSearchSpace::getTrueRoot() const {
	assert(root);
	return root;
}

#endif /* SEARCHSPACE_H_ */
