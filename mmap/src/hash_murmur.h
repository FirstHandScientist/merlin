/*
 * murmur_hash.h
 *
 *  Created on: Apr 18, 2013
 *      Author: radu
 */

#ifndef IBM_ANYTIME_MURMUR_HASH_H_
#define IBM_ANYTIME_MURMUR_HASH_H_

#include "base.h"

#include <stdint.h>

//-----------------------------------------------------------------------------

void MurmurHash3_x86_32  ( const void * key, int len, size_t seed, void * out );

void MurmurHash3_x64_64  ( const void * key, int len, size_t seed, void * out );

void MurmurHash3_x64_128 ( const void * key, int len, size_t seed, void * out );

#endif /* MURMUR_HASH_H_ */
