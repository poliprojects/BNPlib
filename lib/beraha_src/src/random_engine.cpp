//
// Created by mario on 14/05/19.
//

#include "random_engine.hpp"


RandomEngine::RandomEngine(stats::uint_t seed_val) {
  mt.seed(seed_val);
}

stats::rand_engine_t& RandomEngine::get() {
  return mt;
}
