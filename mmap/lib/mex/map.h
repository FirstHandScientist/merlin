#ifndef __MEX_MAP_H
#define __MEX_MAP_H

#include <assert.h>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <stdlib.h>
#include <stdint.h>

#include <map>
#include "mxObject.h"

/*
 * mex::map<Key,Value>      : templated class for map (unique mappings) with matlab containers
 * mex::multimap<Key,Value> : templated class for multimap (non-unique mappings) with matlab containers
*/

namespace mex {

template <class Key, class T>
class map : public virtual mxObject, public std::map<Key,T> { 
public:
	bool checkValid() { return true; };
};

template <class Key, class T>
class multimap : public virtual mxObject, public std::multimap<Key,T> { 
public:
	bool checkValid() { return true; };
};

}       // namespace mex
#endif  // re-include

