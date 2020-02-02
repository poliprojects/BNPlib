#ifndef STAN_MATH_PRIM_MAT_FUN_DISTANCE_HPP
#define STAN_MATH_PRIM_MAT_FUN_DISTANCE_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/squared_distance.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the distance between the specified vectors.
 *
 * @tparam T1 type of elements in first vector
 * @tparam T2 type of elements in second vector
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return Dot product of the vectors.
 * @throw std::domain_error If the vectors are not the same
 * size or if they are both not vector dimensioned.
 */
template <typename T1, int R1, int C1, typename T2, int R2, int C2>
inline return_type_t<T1, T2> distance(const Eigen::Matrix<T1, R1, C1>& v1,
                                      const Eigen::Matrix<T2, R2, C2>& v2) {
  using std::sqrt;
  check_vector("distance", "v1", v1);
  check_vector("distance", "v2", v2);
  check_matching_sizes("distance", "v1", v1, "v2", v2);
  return sqrt(squared_distance(v1, v2));
}

}  // namespace math
}  // namespace stan

#endif
