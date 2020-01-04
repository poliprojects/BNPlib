#ifndef STAN_MATH_PRIM_PROB_BETA_LOG_HPP
#define STAN_MATH_PRIM_PROB_BETA_LOG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/prob/beta_lpdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * The log of the beta density for the specified scalar(s) given the specified
 * sample size(s). y, alpha, or beta can each either be scalar or a vector.
 * Any vector inputs must be the same length.
 *
 * <p> The result log probability is defined to be the sum of
 * the log probabilities for each observation/alpha/beta triple.
 *
 * Prior sample sizes, alpha and beta, must be greater than 0.
 *
 * @deprecated use <code>beta_lpdf</code>
 *
 * @param y (Sequence of) scalar(s).
 * @param alpha (Sequence of) prior sample size(s).
 * @param beta (Sequence of) prior sample size(s).
 * @return The log of the product of densities.
 * @tparam T_y Type of scalar outcome.
 * @tparam T_scale_succ Type of prior scale for successes.
 * @tparam T_scale_fail Type of prior scale for failures.
 */
template <bool propto, typename T_y, typename T_scale_succ,
          typename T_scale_fail>
return_type_t<T_y, T_scale_succ, T_scale_fail> beta_log(
    const T_y& y, const T_scale_succ& alpha, const T_scale_fail& beta) {
  return beta_lpdf<propto, T_y, T_scale_succ, T_scale_fail>(y, alpha, beta);
}

/** \ingroup prob_dists
 * @deprecated use <code>beta_lpdf</code>
 */
template <typename T_y, typename T_scale_succ, typename T_scale_fail>
inline return_type_t<T_y, T_scale_succ, T_scale_fail> beta_log(
    const T_y& y, const T_scale_succ& alpha, const T_scale_fail& beta) {
  return beta_lpdf<T_y, T_scale_succ, T_scale_fail>(y, alpha, beta);
}

}  // namespace math
}  // namespace stan
#endif
