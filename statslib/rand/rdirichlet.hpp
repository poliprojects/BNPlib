#ifndef _statslib_rdirichlet_HPP
#define _statslib_rdirichlet_HPP

// specializations
statslib_inline
arma::vec rdirichlet(const arma::vec& alphas, rand_engine_t& engine);

statslib_inline
arma::vec rdirichlet(const arma::vec& alphas, uint_t seed_val = std::random_device{}());

#include "rdirichlet.ipp"

#endif
