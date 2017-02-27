/*
 * zobrist_hash.cpp
 *
 *  Created on: Apr 18, 2013
 *      Author: radu
 */

#include <assert.h>
#include <stdlib.h>

#include "problem.h"
#include "hash_zobrist.h"

typedef uint64_t * zobrist_c0_array;
typedef uint32_t * zobrist_c1_array;

uint64_t **Zobrist::c0_random_table = NULL;
uint32_t **Zobrist::c1_random_table = NULL;
vector<val_t> Zobrist::random_table_sizes;
size_t Zobrist::total_size = 0;
bool Zobrist::initialization_flag = false;

void Zobrist::initOnce(Problem &problem) {

  assert (!initialization_flag);
  initialization_flag = true;

  size_t size = problem.getN();
  total_size = size;
  c0_random_table = new zobrist_c0_array[size];
  c1_random_table = new zobrist_c1_array[size];

  srandom(1); // kishi: check if this is fine with the orignal source of Random.[cpp,h] (uses boost random for multi-threads)
  for (size_t i = 0; i < size; i ++) {
    size_t n = problem.getDomainSize(i);
    c0_random_table[i] = new uint64_t[n];
    c1_random_table[i] = new uint32_t[n];
    random_table_sizes.push_back(n);
  }
  ; assert(size == random_table_sizes.size());
  initialization_flag = true;

  for(size_t var = 0; var< size; var++) {
    int n = random_table_sizes[var];
    for (int val = 0; val < n; val ++) {
      c0_random_table[var][val] = (static_cast<uint64_t>(random()) << static_cast<uint64_t>(31)) |
	static_cast<uint64_t>(random());
      c1_random_table[var][val] = static_cast<uint32_t>(random());
    }
  }
}

void Zobrist::finishOnce() {
	; assert(initialization_flag);
	if (c0_random_table) {
		for (size_t i = 0; i < total_size; ++i) {
			delete[] c0_random_table[i];
		}
		delete[] c0_random_table;
	}
	if (c1_random_table) {
		for (size_t i = 0; i < total_size; ++i) {
			delete[] c1_random_table[i];
		}
		delete[] c1_random_table;
	}
}

void Zobrist::encodeOR(const std::set<int> &vars, const std::vector<val_t> &assignment) {
  m_c0 = 0;
  m_c1 = 0;

  for (set<int>::const_iterator itC=vars.begin(); itC!=vars.end(); ++itC) {
    int var = *itC;
    ; assert(var >= 0 && static_cast<size_t>(var) < assignment.size());
    int val = assignment[var];
    ; assert(val >= 0 && val < random_table_sizes[var]);

    m_c0 ^= c0_random_table[var][val];
    m_c1 ^= c1_random_table[var][val];
  }
}

void Zobrist::encodeAND(const Zobrist &z, val_t value) {
  *this = z;
  this->m_c0 += value + 1;
}



