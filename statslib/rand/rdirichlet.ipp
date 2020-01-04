statslib_inline
arma::vec rdirichlet(const arma::vec& alphas, rand_engine_t& engine) {
  const uint_t K = alphas.n_elem;
  arma::vec ret(K);
  for (uint_t i = 0; i < K; i++)
    ret(i) = stats::rgamma(alphas(i), 1.0, engine);

  ret /= arma::sum(ret);
  return ret;
}

statslib_inline
arma::vec rdirichlet(const arma::vec& alphas, uint_t seed_val) {
  rand_engine_t engine(seed_val);
  return rdirichlet(alphas, engine);
}
