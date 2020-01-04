/*################################################################################
  ##
  ##   Copyright (C) 2011-2018 Keith O'Hara
  ##
  ##   This file is part of the StatsLib C++ library.
  ##
  ##   Licensed under the Apache License, Version 2.0 (the "License");
  ##   you may not use this file except in compliance with the License.
  ##   You may obtain a copy of the License at
  ##
  ##       http://www.apache.org/licenses/LICENSE-2.0
  ##
  ##   Unless required by applicable law or agreed to in writing, software
  ##   distributed under the License is distributed on an "AS IS" BASIS,
  ##   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  ##   See the License for the specific language governing permissions and
  ##   limitations under the License.
  ##
  ################################################################################*/

/*
 * Sample from a multinomial distribution
 */

#ifndef _statslib_rmultinom_HPP
#define _statslib_rmultinom_HPP

#ifdef STATS_WITH_MATRIX_LIB

template<typename mT, typename eT = double>
statslib_inline
mT rmultinom(const mT& prob);

template<typename mT, typename eT = double>
statslib_inline
uint_t rdiscrete(const mT& prob, uint_t seed_val=std::random_device{}());

template<typename mT, typename eT = double>
statslib_inline
uint_t rdiscrete(const mT& prob, rand_engine_t& engine);

statslib_inline
uint_t rdiscreteunif(uint_t a, uint_t b, uint_t seed_val=std::random_device{}());

statslib_inline
uint_t rdiscreteunif(uint_t a, uint_t b, rand_engine_t& engine);


#include "rmultinom.ipp"

#endif

#endif
