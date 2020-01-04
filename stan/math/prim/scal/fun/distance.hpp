#ifndef STAN_MATH_PRIM_SCAL_FUN_DISTANCE_HPP
#define STAN_MATH_PRIM_SCAL_FUN_DISTANCE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/scal/fun/abs.hpp>

namespace stan {
namespace math {

/**
 * Returns the distance between two scalars.
 *
 * @param x1 First scalar.
 * @param x2 Second scalar.
 * @return Distance between two scalars
 * @throw std::domain_error If the arguments are not finite.
 */
template <typename T1, typename T2>
inline return_type_t<T1, T2> distance(const T1& x1, const T2& x2) {
  check_finite("distance", "x1", x1);
  check_finite("distance", "x2", x2);
  return abs(x1 - x2);
}
}  // namespace math
}  // namespace stan
#endif
