#ifndef STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_EXP_MOD_NORMAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/erfc.hpp>
#include <stan/math/prim/scal/fun/size_zero.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

template <bool propto, typename T_y, typename T_loc, typename T_scale,
          typename T_inv_scale>
return_type_t<T_y, T_loc, T_scale, T_inv_scale> exp_mod_normal_lpdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma,
    const T_inv_scale& lambda) {
  static const char* function = "exp_mod_normal_lpdf";
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale, T_inv_scale>;

  if (size_zero(y, mu, sigma, lambda)) {
    return 0.0;
  }

  T_partials_return logp(0.0);

  check_not_nan(function, "Random variable", y);
  check_finite(function, "Location parameter", mu);
  check_positive_finite(function, "Inv_scale parameter", lambda);
  check_positive_finite(function, "Scale parameter", sigma);
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma, "Inv_scale paramter",
                         lambda);

  if (!include_summand<propto, T_y, T_loc, T_scale, T_inv_scale>::value) {
    return 0.0;
  }

  using std::exp;
  using std::log;
  using std::sqrt;

  operands_and_partials<T_y, T_loc, T_scale, T_inv_scale> ops_partials(
      y, mu, sigma, lambda);

  scalar_seq_view<T_y> y_vec(y);
  scalar_seq_view<T_loc> mu_vec(mu);
  scalar_seq_view<T_scale> sigma_vec(sigma);
  scalar_seq_view<T_inv_scale> lambda_vec(lambda);
  size_t N = max_size(y, mu, sigma, lambda);

  for (size_t n = 0; n < N; n++) {
    const T_partials_return y_dbl = value_of(y_vec[n]);
    const T_partials_return mu_dbl = value_of(mu_vec[n]);
    const T_partials_return sigma_dbl = value_of(sigma_vec[n]);
    const T_partials_return lambda_dbl = value_of(lambda_vec[n]);

    if (include_summand<propto>::value) {
      logp -= LOG_TWO;
    }
    if (include_summand<propto, T_inv_scale>::value) {
      logp += log(lambda_dbl);
    }
    logp += lambda_dbl
                * (mu_dbl + 0.5 * lambda_dbl * sigma_dbl * sigma_dbl - y_dbl)
            + log(erfc((mu_dbl + lambda_dbl * sigma_dbl * sigma_dbl - y_dbl)
                       / (SQRT_TWO * sigma_dbl)));

    const T_partials_return deriv_logerfc
        = NEG_TWO_OVER_SQRT_PI
          * exp(-(mu_dbl + lambda_dbl * sigma_dbl * sigma_dbl - y_dbl)
                / (SQRT_TWO * sigma_dbl)
                * (mu_dbl + lambda_dbl * sigma_dbl * sigma_dbl - y_dbl)
                / (sigma_dbl * SQRT_TWO))
          / erfc((mu_dbl + lambda_dbl * sigma_dbl * sigma_dbl - y_dbl)
                 / (sigma_dbl * SQRT_TWO));

    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_[n]
          += -lambda_dbl + deriv_logerfc * -1.0 / (sigma_dbl * SQRT_TWO);
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_[n]
          += lambda_dbl + deriv_logerfc / (sigma_dbl * SQRT_TWO);
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_[n]
          += sigma_dbl * lambda_dbl * lambda_dbl
             + deriv_logerfc
                   * (-mu_dbl / (sigma_dbl * sigma_dbl * SQRT_TWO)
                      + lambda_dbl / SQRT_TWO
                      + y_dbl / (sigma_dbl * sigma_dbl * SQRT_TWO));
    }
    if (!is_constant_all<T_inv_scale>::value) {
      ops_partials.edge4_.partials_[n]
          += 1 / lambda_dbl + lambda_dbl * sigma_dbl * sigma_dbl + mu_dbl
             - y_dbl + deriv_logerfc * sigma_dbl / SQRT_TWO;
    }
  }
  return ops_partials.build(logp);
}

template <typename T_y, typename T_loc, typename T_scale, typename T_inv_scale>
inline return_type_t<T_y, T_loc, T_scale, T_inv_scale> exp_mod_normal_lpdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma,
    const T_inv_scale& lambda) {
  return exp_mod_normal_lpdf<false>(y, mu, sigma, lambda);
}

}  // namespace math
}  // namespace stan
#endif
