/*
 * random.cpp
 *
 *  Created on: Apr 2, 2013
 *      Author: radu
 */



#include "base.h"

/*
 * Random number generator is static, implemented in base.h.
 * Here only the static member variables are initialized.
 */
int rand::state = 1;
boost::minstd_rand rand::_r;

