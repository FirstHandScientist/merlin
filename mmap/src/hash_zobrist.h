/*
 * zobrist_hash.h
 *
 *  Created on: Apr 18, 2013
 *      Author: radu
 */

#ifndef IBM_MAP_ZOBRIST_HASH_H_
#define IBM_MAP_ZOBRIST_HASH_H_


// Zobrist Hash Key implementation
#include <stdint.h>

#include "base.h"

class Problem;

class Zobrist {

 private:

  uint64_t m_c0;
  uint32_t m_c1;
  static uint64_t **c0_random_table;
  static uint32_t **c1_random_table;
  static std::vector<val_t> random_table_sizes;
  static size_t total_size;
  static bool initialization_flag;

 public:

  uint64_t getC0() const {
    return m_c0;
  }
  uint32_t getC1() const {
    return m_c1;
  }
  bool operator==(const struct Zobrist &other) const {
    return m_c0 == other.m_c0 && m_c1 == other.m_c1;
  }
  Zobrist &operator=(const Zobrist &other) {
    m_c0 = other.m_c0;
    m_c1 = other.m_c1;
    return *this;
  }

  void encodeOR(const std::set<int> &vars, const std::vector<val_t> &assignment);

  void encodeAND(const Zobrist &z, val_t value);

  static void initOnce(Problem &problem);

  static void finishOnce();
};



#endif /* ZOBRIST_HASH_H_ */
