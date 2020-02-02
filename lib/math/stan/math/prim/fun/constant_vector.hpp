#ifndef STAN_MATH_PRIM_FUN_CONSTANT_VECTOR_HPP
#define STAN_MATH_PRIM_FUN_CONSTANT_VECTOR_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return a vector with elements set to the same constant
 *
 * @param K size of the vector
 * @param c constant value
 * @return A vector of size K with all elements initialised to
 * the same constant.
 * @throw std::domain_error if K is negative.
 */
inline Eigen::VectorXd constant_vector(int K, double c) {
  check_nonnegative("constant_vector", "size", K);
  return Eigen::VectorXd::Constant(K, c);
}

}  // namespace math
}  // namespace stan

#endif
