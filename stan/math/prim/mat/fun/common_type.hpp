#ifndef STAN_MATH_PRIM_MAT_FUN_COMMON_TYPE_HPP
#define STAN_MATH_PRIM_MAT_FUN_COMMON_TYPE_HPP

#include <stan/math/prim/arr/fun/common_type.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Struct which calculates type promotion over two types.
 *
 * <p>This specialization is for matrix types.
 *
 * @tparam T1 type of elements in the first matrix
 * @tparam T2 type of elements in the second matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 */
template <typename T1, typename T2, int R, int C>
struct common_type<Eigen::Matrix<T1, R, C>, Eigen::Matrix<T2, R, C> > {
  using type = Eigen::Matrix<typename common_type<T1, T2>::type, R, C>;
};

}  // namespace math
}  // namespace stan

#endif
