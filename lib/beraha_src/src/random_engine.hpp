//
// Created by mario on 14/05/19.
//

#ifndef CPPMODEL_RANDOM_ENGINE_HPP
#define CPPMODEL_RANDOM_ENGINE_HPP

#include <random>
#include "options.hpp"

class RandomEngine {
 public:

  static RandomEngine& Instance() {
    static RandomEngine s;
    return s;
  }
  
  stats::rand_engine_t& get();

 private:
  RandomEngine(stats::uint_t seed_val=std::random_device{}());
  ~RandomEngine() {}

  RandomEngine(RandomEngine const&) = delete;
  RandomEngine& operator= (RandomEngine const&) = delete;

  stats::rand_engine_t mt;

};

#endif //CPPMODEL_RANDOM_ENGINE_HPP
