//
// Created by mario on 14/05/19.
//

#ifndef _statslib_rstickbreak_IPP
#define _statslib_rstickbreak_IPP

statslib_inline
arma::vec
rstickbreak(double a, double b, uint_t numComponents,
            rand_engine_t& engine) {

  arma::vec avec(numComponents - 1);
  avec.fill(a);
  arma::vec bvec(numComponents - 1);
  bvec.fill(b);
  return rstickbreak(avec, bvec, engine);
}

statslib_inline
arma::vec
rstickbreak(double a, double b, uint_t numComponents,
            uint_t seed_val) {
  rand_engine_t engine(seed_val);
  return rstickbreak(a, b, numComponents, engine);
}


statslib_inline
arma::vec
rstickbreak(const arma::vec& avec, const arma::vec& bvec, rand_engine_t& engine) {
  uint_t numComponents = mat_ops::n_elem(avec) + 1;
  arma::vec out(numComponents);
  arma::vec weights(numComponents);

  for (int i = 0; i < numComponents - 1; i++)
    weights(i) = stats::rbeta(avec(i), bvec(i), engine);

  weights(numComponents - 1) = 1.0;
  out(0) = weights(0);

  arma::vec one_minus_weights(numComponents - 1);
  one_minus_weights.fill(1);
  one_minus_weights -= weights.head(numComponents - 1);

  arma::vec cumprods = arma::cumprod(one_minus_weights);

  for (int i = 0; i < numComponents - 1; i++)
    out(i + 1) = weights(i + 1) * cumprods(i);

  return out;
}



statslib_inline
arma::vec rstickbreak(const arma::vec& avec, const arma::vec& bvec,
                      uint_t seed_val) {
  rand_engine_t engine(seed_val);
  return rstickbreak(avec, bvec, engine);
}
#endif
