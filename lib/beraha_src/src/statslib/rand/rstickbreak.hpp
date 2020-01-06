//
// Created by mario on 14/05/19.
//

#ifndef _statslib_rstickbreak_HPP
#define _statslib_rstickbreak_HPP

#ifdef STATS_WITH_MATRIX_LIB

statslib_inline
arma::vec rstickbreak(double a, double b, uint_t numComponents,
                      rand_engine_t& engine);

statslib_inline
arma::vec rstickbreak(double a, double b, uint_t numComponents,
                      uint_t seed_val=std::random_device{}());


statslib_inline
arma::vec rstickbreak(const arma::vec& avec, const arma::vec& bvec,
                      rand_engine_t& engine);

statslib_inline
arma::vec rstickbreak(const arma::vec& avec, const arma::vec& bvec,
                      uint_t seed_val=std::random_device{}());


#include "rstickbreak.ipp"

#endif

#endif